# get l2pfba solutions for dropouts 


import pandas as pd
import os
from itertools import chain
from optlang.symbolics import Zero
from cobra.util.context import get_context
from functools import partial
from micom.solution import solve
from micom.util import check_modification
from micom.problems import regularize_l2_norm
from collections.abc import Sized
import numpy as np
import re
import argparse
from micom import load_pickle

# functions
def load_medium_dict(medium_path):
    """Loads a two-column CSV medium file into a dictionary."""
    medium_df = pd.read_csv(medium_path, names=["reaction", "flux"])
    return dict(zip(medium_df.iloc[:, 0], medium_df.iloc[:, 1]))

def get_biomass_objectives_and_define_as_constraint(community):
    biomass_obj = []
    coefficients = dict()
    n_taxa = len(community.taxonomy)
    
    # get all biomass reactions
    for rec in community.reactions:
        if rec.id.startswith("Growth"):
            biomass_obj.append(rec)
            
    # coefficient scaled to abundance -> for equal abundance = 1/community_size
    for rxn in biomass_obj:
        coefficients[rxn.forward_variable] = 1/n_taxa
        coefficients[rxn.reverse_variable] = -1/n_taxa

    constraint = community.problem.Constraint(0, lb = 0, ub = None)
    community.add_cons_vars(constraint)
    community.solver.update()
    constraint.set_linear_coefficients(coefficients = coefficients)
    return constraint

def add_pfba_objective_totalcom(community, minimal_growth, constraint, atol=1e-4, rtol=1e-4):
    """Add pFBA objective.

    Add objective to minimize the summed flux of all reactions to the
    current objective. This one will work with any objective (even non-linear
    ones).

    See Also
    --------
    pfba

    Parameters
    ----------
    community : micom.Community
        The community to add the objective to.
    """
    constraint.lb = (1-rtol) * minimal_growth - atol
    if community.solver.objective.name == "_pfba_objective":
        raise ValueError("model already has pfba objective")
    reaction_variables = (
        (rxn.forward_variable, rxn.reverse_variable) for rxn in community.reactions
    )
    variables = chain(*reaction_variables)
    community.objective = Zero
    community.objective_direction = "min"
    community.objective.set_linear_coefficients(dict.fromkeys(variables, 1.0))
    if community.modification is None:
        community.modification = "pFBA"
    else:
        community.modification += " and pFBA :)"
    community.solver.update()

def fix_individual_growth_rates(community, growth_rates, atol=1e-6, rtol=1e-6):
    """
    Add one constraint per taxon, fixing its individual biomass (Growth) flux
    close to the value found in `growth_rates` (e.g. from an L2 solution).

    growth_rates : dict {taxon_id: growth_rate}
    Mirrors the same lb/ub bounding pattern used in add_pfba_objective_totalcom.
    """
    constraints = []
    for rxn in community.reactions:
        if not rxn.id.startswith("Growth"):
            continue

        # MICOM tags each per-taxon reaction with the owning taxon id;
        # fall back to string-splitting the reaction id if that's unavailable
        taxon = getattr(rxn, "community_id", None)
        if taxon is None:
            # e.g. "Growth__<taxon>" or "Growth_<taxon>"
            taxon = rxn.id.split("Growth", 1)[1].lstrip("_")

        if taxon not in growth_rates:
            continue

        g = growth_rates[taxon]
        if pd.isna(g):
            continue

        lb = (1 - rtol) * g - atol
        ub = (1 + rtol) * g + atol
        # keep lb <= ub even if g is ~0 or slightly negative due to solver noise
        lb, ub = min(lb, ub), max(lb, ub)

        constraint = community.problem.Constraint(
            rxn.forward_variable - rxn.reverse_variable,
            lb=lb,
            ub=ub,
            name=f"fix_growth_{taxon}",
        )
        community.add_cons_vars(constraint)
        constraints.append(constraint)

    community.solver.update()
    return constraints

def run_l2_then_pfba(community, fractions, fluxes=True, atol=1e-4, rtol=1e-4):
    """
    For each tradeoff fraction:
      1. Fix community-level (weighted) growth at `fraction * max_growth`.
      2. Solve with an L2 objective on individual growth rates -> picks the
         "fair"/non-degenerate distribution of growth across taxa.
      3. Fix each taxon's growth at the value found in step 2.
      4. Switch to a pFBA objective (minimize total flux) and re-solve
         -> resolves remaining flux-level degeneracy without disturbing
         the growth distribution chosen by L2.
    """
    results = []
    opt_sol = community.optimize()
    print("optimal solution:", opt_sol.growth_rate)
    for fr in fractions:
        with community:
            minimal_growth = opt_sol.growth_rate * fr

            # --- Stage 1: L2-regularized growth distribution ---
            regularize_l2_norm(community, 0.0)

            biomass_constraint = get_biomass_objectives_and_define_as_constraint(community)
            biomass_constraint.lb = (1 - rtol) * minimal_growth - atol
            community.solver.update()

            community.solver.problem.parameters.advance.set(0)
            community.solver.problem.cleanup(1e-10)

            sol_l2 = solve(community, fluxes=False, pfba=False, atol=atol, rtol=rtol)

            if sol_l2 is None or getattr(sol_l2, "members", None) is None:
                print(f"Tradeoff {fr}: L2 stage failed, skipping")
                results.append((fr, None))
                continue

            members = sol_l2.members
            if "growth_rate" not in members.columns:
                growth_col = [c for c in members.columns if "growth" in c.lower()]
                if not growth_col:
                    print(f"Tradeoff {fr}: no growth_rate column in L2 solution.members, skipping")
                    results.append((fr, None))
                    continue
                members = members.rename(columns={growth_col[0]: "growth_rate"})
        
            growth_rates = members["growth_rate"].dropna().to_dict()
            growth_rates.pop("medium", None)
            # --- Stage 2: fix individual growth rates, then switch to pFBA ---
            fix_individual_growth_rates(community, growth_rates, atol=atol, rtol=rtol)

            add_pfba_objective_totalcom(
                community, minimal_growth, constraint=biomass_constraint, atol=atol, rtol=rtol
            )

            community.solver.problem.parameters.advance.set(0)
            #community.solver.problem.cleanup(1e-10)

            sol_pfba = community.optimize(fluxes=fluxes, raise_error=False)
            results.append((fr, sol_pfba))

            print("---------")
            print(f"Tradeoff-value: {fr}")
            if sol_pfba is not None:
                print(f"Community objective: {community.solver.variables.community_objective.primal}")
                print(f"pFBA solution (total flux): {sol_pfba.objective_value}")
            else:
                print(f"Solver status: {sol_pfba}")
                print(community.solver.status)
            print("---------")

    results_df = pd.DataFrame.from_records(results, columns=["tradeoff", "solution"])
    return results_df

def test_dropouts(dropout_dir, medium, fractions):
    com_growth_dict = {}
    indi_growth_dict = {}
    exchange_flux_dict = {}

    for file_id in os.listdir(dropout_dir):
        if not file_id.endswith(".pkl"):
            continue

        fp = os.path.join(dropout_dir, file_id)

        try:
            com = load_pickle(fp)
            com.medium = medium
            results_df = run_l2_then_pfba(com, fractions = fractions, fluxes=True)
            results = results_df.iloc[0]["solution"]

            # 1. Community Growth
            com_growth_dict[file_id] = results.growth_rate

            # 2. Individual Growth
            members_df = results.members.dropna(subset=["growth_rate"])
            indi_growth_dict[file_id] = members_df

            # 3. Exchange Fluxes
            fluxes_df = results.fluxes
            ex_cols = [
                col
                for col in fluxes_df.columns
                if isinstance(col, str)
                and (col.startswith("EX_") or col.endswith("_e"))
            ]
            ex_df = fluxes_df[ex_cols]
            exchange_flux_dict[file_id] = ex_df
        except Exception as exc:
            print(f"Skipping {file_id} due to {type(exc).__name__}: {exc}")
            continue

    # 1. Community Growth
    com_growth = pd.DataFrame.from_dict(
        com_growth_dict, orient="index", columns=["community_growth"]
    ) if com_growth_dict else pd.DataFrame(columns=["community_growth"])

    # 2. Individual Growth
    if indi_growth_dict:
        indi_growth = pd.concat(
            indi_growth_dict.values(),
            axis=1,
            keys=indi_growth_dict.keys(),
        )
        indi_growth.index.name = "taxa"
    else:
        indi_growth = pd.DataFrame()

    # 3. Exchange Fluxes
    if exchange_flux_dict:
        ex_df = pd.concat(
            exchange_flux_dict.values(),
            axis=1,
            keys=exchange_flux_dict.keys(),
        )
        ex_df.index.name = "metabolites"
    else:
        ex_df = pd.DataFrame()

    return com_growth, indi_growth, ex_df


def process_and_save(dropout_dir, medium_dict, fractions, save_dir=None):
    print(f"Processing directory: {dropout_dir}")
    com_gr, indi_gr, fluxes = test_dropouts(dropout_dir, medium_dict, fractions)

    if save_dir is None:
        save_dir = dropout_dir
    else:
        os.makedirs(save_dir, exist_ok=True)

    # Define output file paths in the requested save directory
    com_gr_path = os.path.join(save_dir, "community_growth.csv")
    indi_gr_path = os.path.join(save_dir, "individual_growth.csv")
    fluxes_path = os.path.join(save_dir, "exchange_fluxes.csv")

    # Save to CSV
    com_gr.to_csv(com_gr_path)
    indi_gr.to_csv(indi_gr_path)
    fluxes.to_csv(fluxes_path)

    print(f"Saved results to:\n - {com_gr_path}\n - {indi_gr_path}\n - {fluxes_path}\n")

# MAIN 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Optimize community dropouts."
    )
    parser.add_argument(
        "--dropout_dirs",
        nargs="+",
        required=True,
        help="One or more directories containing .pkl dropout communities.",
    )
    parser.add_argument(
        "--medium",
        required=True,
        help="Path to the medium CSV file (reaction, flux).",
    )
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=[1.0],
        help="Tradeoff fractions to evaluate for the L2->pFBA workflow. Defaults to [1.0].",
    )
    parser.add_argument(
        "--save_dir",
        default=None,
        help="Optional directory where output CSV files should be written. Defaults to the dropout directory when omitted.",
    )

    args = parser.parse_args()

    # Load medium dictionary once
    medium_dict = load_medium_dict(args.medium)
    fractions = args.fractions

    # Loop through all directories passed via CLI
    for dropout_dir in args.dropout_dirs:
        if os.path.exists(dropout_dir):
            process_and_save(dropout_dir, medium_dict, fractions, args.save_dir)
        else:
            print(f"Warning: Directory not found -> {dropout_dir}")