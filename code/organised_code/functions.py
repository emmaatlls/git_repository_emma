# functions used in several notebooks

# imports
import os
import pandas as pd
from collections import defaultdict
from cobra import Reaction, Metabolite

##### media generation and saving


#create media dict from dataframe
def create_medium(media_df):
    media_dict = dict(zip(media_df.iloc[:,0], media_df.iloc[:,1]))
    return media_dict

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