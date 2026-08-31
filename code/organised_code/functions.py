# functions used in several notebooks

# imports
import os
import pandas as pd
from collections import defaultdict
from cobra import Reaction, Metabolite
from itertools import chain, combinations
from optlang.symbolics import Zero
from micom.solution import solve
from micom.problems import regularize_l2_norm

##### media generation and saving


#create media dict from dataframe
def create_medium(media_df):
    media_dict = dict(zip(media_df.iloc[:,0], media_df.iloc[:,1]))
    return media_dict

# compute medium series from csv (needed for complete_medium)
def compute_medium_series(csv_medium):
    new_id = csv_medium.copy()
    new_id[0] = new_id[0].astype(str).str.strip().str.replace(r'\d+$', '', regex=True)
    new_id[0] = new_id[0].str.replace('_e', '_m')
    flux_values = pd.to_numeric(new_id[1], errors='coerce').fillna(1000.0)
    med_series = pd.Series(flux_values.values, index=new_id[0].astype(str).str.strip())
    med_series = med_series[med_series.index.notna() & (med_series.index != 'nan')]
    return med_series

# save medium locally
def save_media(medium_dict, medianame, media_save_dir):
    media_df = pd.DataFrame.from_dict(
        medium_dict, orient="index")
    filename = f"{medianame}.csv"
    filepath = os.path.join(media_save_dir, filename)
    media_df.to_csv(filepath)    
    return media_df

# assign new medium in cobra
def get_exchanges(model):
    model_exchanges = []
    for ex in model.exchanges:
        ex_id = str(ex).split(' ',1)[0].replace(':', '')
        model_exchanges.append(ex_id)
    return model_exchanges

def get_model_specific_medium(model, medium_dict):
    exchanges = get_exchanges(model)
    temp_li = []
    for k, v in medium_dict.items():
        if k in exchanges:
            temp_li.append(k)
    specific_medium = {i: medium_dict[i] for i in temp_li}
    return specific_medium

# get metabolite classes per media
def assign_metaboliteclass(rxn_df, met_class_dict):
    medium_met_class = defaultdict(list)
    count_dict = defaultdict(int)

    if isinstance(rxn_df, pd.DataFrame):
        if "reaction" in rxn_df.columns:
            reactions = rxn_df["reaction"].astype(str)
        else:
            reactions = rxn_df.index.astype(str)
    else:
        reactions = pd.Index(rxn_df).astype(str)

    for ex_r in reactions:
        if ex_r.endswith("_e"):
            lookup_r = ex_r
        else:
            lookup_r = ex_r[:-2] + "_e"

        if lookup_r in met_class_dict:
            mc = met_class_dict[lookup_r]
            medium_met_class[mc].append(lookup_r)
            count_dict[mc] += 1
        else:
            print(f"{lookup_r} could not be found in metabolite lookup")

    return medium_met_class, count_dict

#compare two media 
def compare_media(source_dict, target_dict):
    only_in_source = set()
    in_both = set()
    all_source_converted = set()
    for ex_r, flux in source_dict.items():
        com_rxn_id = ex_r[:-2] + "_m"
        all_source_converted.add(com_rxn_id)
        
        if com_rxn_id in target_dict:
            in_both.add(com_rxn_id)
        else:
            only_in_source.add(com_rxn_id)
                
    only_in_target = set(target_dict.keys()) - all_source_converted
    
    return list(only_in_source), list(in_both), list(only_in_target)


def tidy_exchange_fluxes(fluxes_df: pd.DataFrame) -> pd.DataFrame:
    """Converts a MICOM reaction flux matrix into a long-format DataFrame with

    columns: [reaction, taxon, flux, direction].
    """
    df = fluxes_df.copy()

    # 1. Drop non-taxon rows (e.g. "medium") directly from index
    df = df[~df.index.isin(["medium"])].copy()

    # 2. Extract taxon/compartment name from index into a column
    df["taxon"] = df.index

    # 3. Filter for exchange/sink reactions (columns starting with EX_ and ending with _e)
    ex_cols = [
        col
        for col in df.columns
        if isinstance(col, str) and col.startswith("EX_") and col.endswith("_e")
    ]
    df = df[["taxon"] + ex_cols]

    # 4. Melt from wide matrix to long tidy format
    tidy = df.melt(id_vars=["taxon"], var_name="reaction", value_name="flux")

    # 5. Drop NaN and zero fluxes (inactive reactions)
    tidy = tidy.dropna(subset=["flux"])
    tidy = tidy[tidy["flux"] != 0.0].copy()

    # 6. Assign direction based on standard COBRA/MICOM conventions:
    #    Negative flux (< 0) -> Uptake / Import into system
    #    Positive flux (> 0) -> Secretion / Export out of system
    tidy["direction"] = tidy["flux"].apply(
        lambda x: "import" if x < 0 else "export"
    )

    # 7. Reorder columns as requested
    tidy = tidy[["reaction", "taxon", "flux", "direction"]]

    return tidy.reset_index(drop=True)
def classify_interaction(
    met_df: pd.DataFrame, taxon: str, partner: str, abundance: float #hardcode abundance as all strains are equally abundant in my case 
) -> pd.DataFrame:
    """Checks if and how two taxa interact for a specific metabolite."""
    if met_df.empty:
        return None

    tol = 1e-06
    
    # Keep active fluxes exceeding tolerance threshold
    f = met_df[(met_df.flux.abs() * abundance) > tol]
    
    # Must involve at least 2 active taxa and cannot both be exporting without consumer
    if (f.shape[0] < 2) or (f.direction == "export").all():
        return None
        
    if (f["direction"] == "import").sum() == 2:
        int_type = "co-consumed"
    else:
        focal_dirs = f.loc[f["taxon"] == taxon, "direction"]
        if not focal_dirs.empty and (focal_dirs == "export").all():
            int_type = "provided"
        else:
            int_type = "received"

    return pd.DataFrame(
        {
            "focal": taxon,
            "partner": partner,
            "reaction": met_df["reaction"].iloc[0],
            "class": int_type,
            "flux": (f.flux.abs() * abundance).min(),
        },
        index=[0],
    )
def interactions_per_taxon(clean_flux_df: pd.DataFrame, taxon: str) -> pd.DataFrame:
    """Quantifies interactions for a single focal taxon against all partner taxa."""
    # Find metabolites where the focal taxon has active exchange
    focal_mets = clean_flux_df.loc[
        clean_flux_df["taxon"] == taxon, "reaction"
    ].unique()
    
    # Get all potential partner taxa
    partners = [t for t in clean_flux_df["taxon"].unique() if t not in (taxon, "medium")]
    abundance = 1/len(partners)
    ints = []
    for p in partners:
        # Subset fluxes for focal taxon and current partner
        pair_df = clean_flux_df[
            clean_flux_df["taxon"].isin([taxon, p]) & 
            clean_flux_df["reaction"].isin(focal_mets)
        ]
        
        if pair_df.empty:
            continue

        res = (
            pair_df.groupby("reaction", group_keys=False)
            .apply(lambda df: classify_interaction(df, taxon, p, abundance))
        )
        ints.append(res)
        
    # Drop empty or None results
    valid_ints = [i for i in ints if i is not None and not i.empty]
    return pd.concat(valid_ints, ignore_index=True) if valid_ints else pd.DataFrame()

def get_interactions(clean_flux_df: pd.DataFrame) -> pd.DataFrame:
    """Runs interaction pipeline across all unique taxa."""
    # Ensure taxon values are treated consistently as strings
    df = clean_flux_df.copy()
    df["taxon"] = df["taxon"].astype(str)
    
    taxons = [t for t in df["taxon"].unique() if t != "medium"]
    
    all_ints = []
    for taxon in taxons:
        res = interactions_per_taxon(df, taxon)
        if not res.empty:
            all_ints.append(res)
            
    return pd.concat(all_ints, ignore_index=True) if all_ints else pd.DataFrame()

    
##### add metabolites and reactions to models
# taken from Lisa 
def add_new_rxn(model, id, name, lb, ub, stoich):
    if id not in model.reactions:
        new_rxn = Reaction(id=id, name=name, lower_bound=lb, upper_bound=ub)

        met_objs = {}
        for met_id, coeff in stoich.items():
            met_obj = model.metabolites.get_by_id(met_id) if met_id in model.metabolites else None
            if met_obj is None:
                return f"{met_id} not in {model.id}; {id} was not added"
            met_objs[met_obj] = coeff

        new_rxn.add_metabolites(met_objs)
        model.add_reactions([new_rxn])
        return None
    else:
        return f"{id} already in {model.id}"
    
# taken from Lisa    
def add_new_met(model, id, name, formula, charge, compartment):
    if id not in model.metabolites:
        new_met = Metabolite(id, name=name, formula=formula, charge=charge, compartment=compartment)
        model.add_metabolites([new_met])
        return None
    else:
        return f"{id} already in {model.id}"
    
def overwrite_reaction(model, rxn_id, new_rxn_dict):
    if rxn_id in model.reactions:
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            old_metabolites = {met.id: coeff for met, coeff in rxn.metabolites.items()}
            rxn.subtract_metabolites(rxn.metabolites)
            rxn.add_metabolites(new_rxn_dict)
            new_metabolites = {met.id: coeff for met, coeff in rxn.metabolites.items()}
        except KeyError as e:
            print(f"{model} does not contain one of the metabolites, {e}")

def delete_reaction(model, rxn_id):
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        old_metabolites = {met.id: coeff for met, coeff in rxn.metabolites.items()}
        model.remove_reactions([rxn])



###################################################################################################################
##### Growth ########
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
def run_fba(community, fractions, fluxes):
    results = []
    opt_sol = community.optimize()
    for fr in fractions:
        with community:
            minimal_growth = opt_sol.growth_rate * fr
            biomass_constraint = get_biomass_objectives_and_define_as_constraint(community)
            add_pfba_objective_totalcom(community, minimal_growth, constraint = biomass_constraint, atol=1e-6, rtol=1e-6)

            community.solver.problem.parameters.advance.set(0)
            community.solver.problem.cleanup(1e-10)

            sol_pfba2 = community.optimize(fluxes=fluxes, raise_error = False)
            results.append((fr, sol_pfba2))

            print("---------")
            print(f"Tradeoff-value: {fr}")
            if sol_pfba2 != None:
                print(f"Community {community.solver.variables.community_objective.primal}")
                print(f"pFBA solution: {sol_pfba2.objective_value}")
            else:
                print(f"Solver status: {sol_pfba2}")
            print("---------")

    results_df = pd.DataFrame.from_records(results, columns=["tradeoff", "solution"])
    return results_df
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
    opt_sol = community.optimize() #FBA 
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

def save_tradeoff_results(results_df, save_dir, prefix):
    """Saves community growth rates, individual taxon growth rates, and flux tables

    from tradeoff results to CSV files in `save_dir`.
    """
    save_path = save_dir

    summary_records = []
    taxon_growth_list = []

    for _, row in results_df.iterrows():
        tradeoff = row["tradeoff"]
        solution = row["solution"]

        # Skip failed/empty solutions
        if solution is None:
            continue

        # 1. Capture community growth rate & objective
        com_growth = getattr(solution, "growth_rate", None)
        obj_val = getattr(solution, "objective_value", None)

        summary_records.append(
            {
                "tradeoff": tradeoff,
                "community_growth": com_growth,
                "objective_value": obj_val,
                "status": getattr(solution, "status", "unknown"),
            }
        )

        # 2. Capture individual taxon growth rates (from solution.members)
        if hasattr(solution, "members") and solution.members is not None:
            members_df = solution.members.copy()
            if isinstance(members_df, pd.DataFrame):
                members_df["tradeoff"] = tradeoff
                taxon_growth_list.append(members_df)
                tradeoff_str = f"{int(tradeoff * 10):02d}"
                members_df.to_csv(os.path.join(save_path, f"{prefix}_{tradeoff_str}_individualtaxon_growth.csv"))


        # 3. Capture reaction flux matrix (from solution.fluxes)
        if hasattr(solution, "fluxes") and solution.fluxes is not None:
            fluxes_df = solution.fluxes.copy()
            if isinstance(fluxes_df, pd.DataFrame):
                fluxes_df["tradeoff"] = tradeoff
                tradeoff_str = f"{int(tradeoff * 10):02d}"
                fluxes_df.to_csv(os.path.join(save_path, f"{prefix}_{tradeoff_str}_fluxes.csv"))

    # --- Save Summary Table ---
    if summary_records:
        summary_df = pd.DataFrame(summary_records)
        summary_df.to_csv(os.path.join(save_path, f"{prefix}_community_growth_summary.csv"), index=False)

    # --- Save Taxon Growth Rates ---
    if taxon_growth_list:
        all_taxon_growth = pd.concat(taxon_growth_list, ignore_index=False)
        all_taxon_growth.to_csv(os.path.join(save_path, f"{prefix}_all_taxon_growth.csv"), index=False)

    print(f"Results successfully saved to: {save_path}")



###########################################################################################################################################
################### INTERACTION ANALYSIS ########################################################################################################################
def tidy_exchange_fluxes(fluxes_df: pd.DataFrame) -> pd.DataFrame:
    """Converts a MICOM reaction flux matrix into a long-format DataFrame with

    columns: [reaction, taxon, flux, direction].
    """
    df = fluxes_df.copy()

    # 1. Drop non-taxon rows (e.g. "medium") directly from index
    df = df[~df.index.isin(["medium"])].copy()

    # 2. Extract taxon/compartment name from index into a column
    df["taxon"] = df.index

    # 3. Filter for exchange/sink reactions (columns starting with EX_ and ending with _e)
    ex_cols = [
        col
        for col in df.columns
        if isinstance(col, str) and col.startswith("EX_") and col.endswith("_e")
    ]
    df = df[["taxon"] + ex_cols]

    # 4. Melt from wide matrix to long tidy format
    tidy = df.melt(id_vars=["taxon"], var_name="reaction", value_name="flux")

    # 5. Drop NaN and zero fluxes (inactive reactions)
    tidy = tidy.dropna(subset=["flux"])
    tidy = tidy[tidy["flux"] != 0.0].copy()

    # 6. Assign direction based on standard COBRA/MICOM conventions:
    #    Negative flux (< 0) -> Uptake / Import into system
    #    Positive flux (> 0) -> Secretion / Export out of system
    tidy["direction"] = tidy["flux"].apply(
        lambda x: "import" if x < 0 else "export"
    )

    # 7. Reorder columns as requested
    tidy = tidy[["reaction", "taxon", "flux", "direction"]]

    return tidy.reset_index(drop=True)
def classify_interaction(
    met_df: pd.DataFrame, taxon: str, partner: str, abundance: float #hardcode abundance as all strains are equally abundant in my case 
) -> pd.DataFrame:
    """Checks if and how two taxa interact for a specific metabolite."""
    if met_df.empty:
        return None

    tol = 1e-06
    
    # Keep active fluxes exceeding tolerance threshold
    f = met_df[(met_df.flux.abs() * abundance) > tol]
    
    # Must involve at least 2 active taxa and cannot both be exporting without consumer
    if (f.shape[0] < 2) or (f.direction == "export").all():
        return None
        
    if (f["direction"] == "import").sum() == 2:
        int_type = "co-consumed"
    else:
        focal_dirs = f.loc[f["taxon"] == taxon, "direction"]
        if not focal_dirs.empty and (focal_dirs == "export").all():
            int_type = "provided"
        else:
            int_type = "received"

    return pd.DataFrame(
        {
            "focal": taxon,
            "partner": partner,
            "reaction": met_df["reaction"].iloc[0],
            "class": int_type,
            "flux": (f.flux.abs() * abundance).min(),
        },
        index=[0],
    )
def interactions_per_taxon(clean_flux_df: pd.DataFrame, taxon: str) -> pd.DataFrame:
    """Quantifies interactions for a single focal taxon against all partner taxa."""
    # Find metabolites where the focal taxon has active exchange
    focal_mets = clean_flux_df.loc[
        clean_flux_df["taxon"] == taxon, "reaction"
    ].unique()
    
    # Get all potential partner taxa
    partners = [t for t in clean_flux_df["taxon"].unique() if t not in (taxon, "medium")]
    abundance = 1/len(partners)
    ints = []
    for p in partners:
        # Subset fluxes for focal taxon and current partner
        pair_df = clean_flux_df[
            clean_flux_df["taxon"].isin([taxon, p]) & 
            clean_flux_df["reaction"].isin(focal_mets)
        ]
        
        if pair_df.empty:
            continue

        res = (
            pair_df.groupby("reaction", group_keys=False)
            .apply(lambda df: classify_interaction(df, taxon, p, abundance))
        )
        ints.append(res)
        
    # Drop empty or None results
    valid_ints = [i for i in ints if i is not None and not i.empty]
    return pd.concat(valid_ints, ignore_index=True) if valid_ints else pd.DataFrame()

def get_interactions(clean_flux_df: pd.DataFrame) -> pd.DataFrame:
    """Runs interaction pipeline across all unique taxa."""
    # Ensure taxon values are treated consistently as strings
    df = clean_flux_df.copy()
    df["taxon"] = df["taxon"].astype(str)
    
    taxons = [t for t in df["taxon"].unique() if t != "medium"]
    
    all_ints = []
    for taxon in taxons:
        res = interactions_per_taxon(df, taxon)
        if not res.empty:
            all_ints.append(res)
            
    return pd.concat(all_ints, ignore_index=True) if all_ints else pd.DataFrame()

# MRO = metabolic resource overlap
def mro_score(taxa_uptakes):
    """Calculates the Metabolic Resource Overlap (MRO) score for a community.

    Parameters:
    -----------
    uptake_df : pd.DataFrame
        DataFrame containing uptake data with columns for taxon/species and reactions/metabolites.
    taxon_col : str
        Column name identifying the taxa (default: "taxon").
    reaction_col : str
        Column name identifying the exchange reactions or imported metabolites.

    Returns:
    --------
    float
        The calculated MRO score between 0.0 (no overlap) and 1.0 (identical resource usage).
    """
    # Group uptake reactions (as sets) by taxon

    taxa = list(taxa_uptakes.index)
    num_taxa = len(taxa)

    # Total sum of uptakes across all individual species 
    total_uptakes_sum = sum(len(uptakes) for uptakes in taxa_uptakes)
    if total_uptakes_sum == 0:
        return 0.0

    # pairwise intersection sum
    intersection_sum = 0
    for strain_a, strain_b in combinations(taxa, 2):
        shared_uptakes = taxa_uptakes[strain_a].intersection(taxa_uptakes[strain_b])
        intersection_sum += len(shared_uptakes)

    # Number of pairwise combinations
    num_combinations = (num_taxa * (num_taxa - 1)) / 2

    # Final MRO Score
    mro = num_taxa * intersection_sum / (num_combinations * total_uptakes_sum)
    return mro


# compare uptakes from medium to uptakes from community 
def uptake_mediumvscommunity(uptake_df, medium_dict):
    count_medium_up = []
    all_uptakes = set(uptake_df["reaction"].dropna())
    medium_rxns = {rxn[:-2] + "_e" if isinstance(rxn, str) and rxn.endswith("_m") else rxn for rxn in medium_dict.keys()}
    for up_rxn in all_uptakes:
        if up_rxn not in count_medium_up and up_rxn in medium_rxns:
            count_medium_up.append(up_rxn)
    count_community_exchanges = len(all_uptakes) - len(count_medium_up)
    rel_com = count_community_exchanges
    rel_med = len(count_medium_up)
    mip = (len(all_uptakes) - rel_med)/len(all_uptakes)
    return rel_med, count_medium_up, rel_com, mip

def fMIP(uptake_df, medium_dict):
    fluxes_all = uptake_df["flux"].sum()
    fluxes_medium = 0
    medium_rxns = {rxn[:-2] + "_e" if isinstance(rxn, str) and rxn.endswith("_m") else rxn for rxn in medium_dict.keys()}
    in_medium_mask = uptake_df["reaction"].isin(medium_rxns)
    fluxes_medium = uptake_df.loc[in_medium_mask, "flux"].sum()
    fmip_score = (abs(fluxes_all) - abs(fluxes_medium)) / abs(fluxes_all)
    return fluxes_all, fluxes_medium, fmip_score