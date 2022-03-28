import pandas as pd

#
# Leaving only C atoms
#
atom_mapping = pd.read_csv("rdt_analysis/atom_mapping.csv", sep=',')
print(atom_mapping)
atom_mapping = atom_mapping[["reaction", "reactant", "num", "product"]]

cond = atom_mapping[atom_mapping['reactant'].str.contains('_C_') & atom_mapping['product'].str.contains('_C_')]
cond.to_csv(r'data/atom_mapping_C_atoms.csv', index=False)


#
# From LR reactions to master reactions
#

LR_atom_mapping = pd.read_csv("data/atom_mapping_C_atoms.csv", sep=",")
LR_atom_mapping.columns = ["LR_Rhea_ID", "reactant", "num", "product"]
LR_atom_mapping = LR_atom_mapping.dropna()

Rhea_RDT_mapping = pd.read_csv("data/undirected_to_LR.tsv", sep="\t")
Rhea_RDT_mapping.columns = ["reaction", "LR_Rhea_ID"]

atom_mapping = pd.merge(LR_atom_mapping, Rhea_RDT_mapping, how="left")
atom_mapping = atom_mapping[["reaction", "reactant", "num", "product"]]

atom_mapping = atom_mapping.dropna()
atom_mapping = atom_mapping.reset_index(drop=True)


#
# Replacing polymers IDs with ChEBI IDs
#

polymers_x = pd.read_csv("pre_data/polymers_to_chebi.tsv", sep="\t")
polymers_x = polymers_x[["Polymer_ID", "ChEBI_ID"]]
polymers_x.columns = ["metabolite.x", "metabolite.x.chebi"]

polymers_y = pd.read_csv("pre_data/polymers_to_chebi.tsv", sep="\t")
polymers_y = polymers_y[["Polymer_ID", "ChEBI_ID"]]
polymers_y.columns = ["metabolite.y", "metabolite.y.chebi"]


#
# Making supplementary files
#

tmp = atom_mapping["reactant"].str.split("_", n=2, expand=True)
tmp.columns = ["metabolite.x", "element.x", "coordinates.x"]

tmp = pd.merge(tmp, polymers_x, how="left")
tmp.loc[tmp["metabolite.x"].str.contains("POLYMER:"), "metabolite.x"] = tmp["metabolite.x.chebi"]
tmp = tmp[["metabolite.x", "element.x", "coordinates.x"]]

tmp["atom.x"] = tmp[["metabolite.x", "coordinates.x"]].agg("_".join, axis=1)

atom_mapping["metabolite.x"] = tmp["metabolite.x"]
atom_mapping["element.x"] = tmp["element.x"]
atom_mapping["atom.x"] = tmp["atom.x"]


tmp = atom_mapping["product"].str.split("_", n=2, expand=True)
tmp.columns = ["metabolite.y", "element.y", "coordinates.y"]

tmp = pd.merge(tmp, polymers_y, how="left")
tmp.loc[tmp["metabolite.y"].str.contains("POLYMER:"), "metabolite.y"] = tmp["metabolite.y.chebi"]
tmp = tmp[["metabolite.y", "element.y", "coordinates.y"]]

tmp["atom.y"] = tmp[["metabolite.y", "coordinates.y"]].agg("_".join, axis=1)

atom_mapping["metabolite.y"] = tmp["metabolite.y"]
atom_mapping["element.y"] = tmp["element.y"]
atom_mapping["coordinates.y"] = tmp["coordinates.y"]
atom_mapping["atom.y"] = tmp["atom.y"]


atom_mapping = atom_mapping[["reaction", "metabolite.x", "element.x", "atom.x", "metabolite.y", "element.y", "atom.y"]]
atom_mapping = atom_mapping.drop_duplicates()
atom_mapping = atom_mapping.reset_index(drop=True)

full_atom_mapping = atom_mapping
atom_mapping = atom_mapping[["reaction", "metabolite.x", "metabolite.y"]]

atom_mapping = atom_mapping.drop_duplicates()
atom_mapping = atom_mapping.reset_index(drop=True)

atom_mapping["index"] = atom_mapping.index
atom_mapping["index"] = atom_mapping["index"] + 1
atom_mapping["index"] = atom_mapping["index"].astype(str)

atom_mapping["length"] = 5 - atom_mapping["index"].str.len()

atom_mapping["rpair"] = atom_mapping.length.apply(lambda x: ("RP" + "0" * x))
atom_mapping["rpair"] = atom_mapping["rpair"] + atom_mapping["index"]

atom_mapping = atom_mapping[["reaction", "metabolite.x", "metabolite.y", "rpair"]]
atom_mapping.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair"]

atom_mapping = pd.DataFrame(atom_mapping)
atom_mapping.to_csv(r"data/rpairs_preprocessing.csv", index=False)

atom_mapping.columns = ["reaction", "metabolite.x", "metabolite.y", "rpair"]
full_atom_mapping_m = pd.merge(full_atom_mapping, atom_mapping,
                               how="left",
                               left_on=["reaction", "metabolite.x", "metabolite.y"],
                               right_on=["reaction", "metabolite.x", "metabolite.y"])

# rpair2align = full_atom_mapping_m[["rpair", "atom.x", "atom.y"]]
# rpair2align = rpair2align.drop_duplicates()
# rpair2align = rpair2align.reset_index(drop=True)
# rpair2align.to_csv(r"network_files/rpair2align.csv", index=False)
# 
# reaction2pair = full_atom_mapping_m[["rpair", "reaction"]]
# reaction2pair = reaction2pair.drop_duplicates()
# reaction2pair = reaction2pair.reset_index(drop=True)
# reaction2pair.to_csv(r"network_files/reaction2pair.csv", index=False)

reaction2align = full_atom_mapping_m[["reaction", "atom.x", "atom.y"]]
reaction2align = reaction2align.drop_duplicates()
reaction2align = reaction2align.reset_index(drop=True)
reaction2align.to_csv(r"network_files/reaction2align.csv", index=False)

atoms_x = full_atom_mapping_m[["atom.x", "metabolite.x", "element.x"]]
atoms_x.columns = ["atom", "metabolite", "element"]
atoms_y = full_atom_mapping_m[["atom.y", "metabolite.y", "element.y"]]
atoms_y.columns = ["atom", "metabolite", "element"]

atoms = pd.concat([atoms_x, atoms_y])
atoms = atoms.drop_duplicates()
atoms = atoms.reset_index(drop=True)
atoms.to_csv(r"network_files/atoms.csv", index=False)

metabolite2atom_x = full_atom_mapping_m[["atom.x", "metabolite.x"]]
metabolite2atom_x.columns = ["atom", "metabolite"]
metabolite2atom_y = full_atom_mapping_m[["atom.y", "metabolite.y"]]
metabolite2atom_y.columns = ["atom", "metabolite"]

metabolite2atom = pd.concat([metabolite2atom_x, metabolite2atom_y])
metabolite2atom = metabolite2atom.drop_duplicates()
metabolite2atom = metabolite2atom.reset_index(drop=True)
metabolite2atom.to_csv(r"network_files/metabolite2atom.csv", index=False)


#
# Making RPAIRS
#

full_atom_mapping_m.to_csv(r"network_files/full_atom_mapping_1.csv", index=False)

full_atom_mapping_m["el.check.x"] = 1
full_atom_mapping_m["el.check.y"] = 1

full_atom_mapping_m.loc[full_atom_mapping_m["element.y"] != "C", "el.check.y"] = 0

full_atom_mapping = full_atom_mapping_m[["reaction", "metabolite.x", "metabolite.y", "rpair", "el.check.x", "el.check.y"]]

full_atom_mapping = full_atom_mapping.groupby(['reaction', 'metabolite.x', 'metabolite.y', 'rpair'],
                                              as_index=False).agg({'el.check.x': "sum",
                                                                   'el.check.y': "sum"})

full_atom_mapping.loc[full_atom_mapping['el.check.x'] >= 1, 'check'] = 1
full_atom_mapping.loc[full_atom_mapping['el.check.y'] >= 1, 'check'] = 1
full_atom_mapping.loc[full_atom_mapping['el.check.x'] < 1, 'check'] = 0
full_atom_mapping.loc[full_atom_mapping['el.check.y'] < 1, 'check'] = 0

full_atom_mapping['rtype'] = "not_main"
full_atom_mapping.loc[full_atom_mapping['check'] == 1, 'rtype'] = "main"

full_atom_mapping = full_atom_mapping[["reaction", "metabolite.x", "el.check.x", "metabolite.y", "el.check.y", "rpair", "rtype"]]
full_atom_mapping.columns = ["rxn", "metabolite.x", "el.check.x", "metabolite.y", "el.check.y", "rpair", "rtype"]

atom_mapping.columns = ['rxn', 'metabolite.x', 'metabolite.y', "rpair"]

merged = pd.merge(atom_mapping, full_atom_mapping, how='left',
                  left_on=['rxn', 'metabolite.x', 'metabolite.y', "rpair"],
                  right_on=['rxn', 'metabolite.x', 'metabolite.y', "rpair"])

merged = merged[["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]]
merged.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]

# # Deleting CO2
# merged.loc[merged['metabolite.x'] == "CHEBI:16526", 'rtype'] = "not_main"
# merged.loc[merged['metabolite.y'] == "CHEBI:16526", 'rtype'] = "not_main"
#
# # Deleting Coenzyme A (CHEBI:57287)
# merged.loc[merged['metabolite.x'] == "CHEBI:57287", 'rtype'] = "not_main"
# merged.loc[merged['metabolite.y'] == "CHEBI:57287", 'rtype'] = "not_main"

merged.to_csv(r'network_files/rpairs_preprocessing_2.csv', index=False)


#
# Fixing RPAIRS file
#

rpairs = pd.read_csv('network_files/rpairs_preprocessing_2.csv', sep=",")
rpairs = rpairs[rpairs.rtype == "main"]
rpairs = rpairs.reset_index(drop=True)

Rhea_RDT_mapping = pd.read_csv("data/undirected_to_LR.tsv", sep="\t")
Rhea_RDT_mapping = Rhea_RDT_mapping[["Rhea_ID", "LR_Rhea_ID"]]
Rhea_RDT_mapping.columns = ["Rhea_ID", "rxn"]

rpairs = pd.merge(rpairs, Rhea_RDT_mapping, how="left")
rpairs = rpairs[["Rhea_ID", "metabolite.x", "metabolite.y", "rpair", "rtype"]]
rpairs.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]


# Fixing order of reactants
rpairs["bad_order"] = 0
rpairs.loc[rpairs['metabolite.x'] > rpairs['metabolite.y'], 'bad_order'] = 1

good = rpairs[rpairs.bad_order == 0]
good = good[["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]]

bad = rpairs[rpairs.bad_order == 1]
bad = bad[["rxn", "metabolite.y", "metabolite.x", "rpair", "rtype"]]
bad.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]

rpairs = pd.concat([good, bad])

rpairs.to_csv(r'network_files/rpairs.csv', index=False)
