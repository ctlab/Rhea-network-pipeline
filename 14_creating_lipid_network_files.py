import pandas as pd


# Reading metabolites data
metabolites = pd.read_csv("network_files/annotation_of_nodes_master_graph.csv", sep=";")
metabolites = metabolites[metabolites.SwissLipids_ID != "-"]
metabolites = metabolites.reset_index(drop=True)
lipids = list(set(metabolites.metabolite.values))

# Reading atom_mapping data
atom_mapping = pd.read_csv("network_files/full_atom_mapping_1.csv", sep=",")
am = atom_mapping[(atom_mapping["metabolite.x"].isin(lipids)) | (atom_mapping["metabolite.y"].isin(lipids))]
am = am.reset_index(drop=True)
am.to_csv(r'network_files/atom_mapping_lipids.tsv', sep="\t", index=False, )

# Leaving only lipids
ChEBI = list(set(am["metabolite.x"].values).union(set(am["metabolite.y"].values)))
metabolites = metabolites[metabolites.metabolite.isin(ChEBI)]
metabolites = metabolites.reset_index(drop=True)

# Get lipid reaction list
lipid_reactions = list(set(am.reaction.values))

# Creating url annotation for lipid reactions
urls = pd.read_csv("network_files/reaction_urls.tsv", sep="\t")
urls = urls[urls.Rhea_ID.isin(lipid_reactions)]
urls = urls.reset_index(drop=True)
urls.to_csv(r'network_files/reaction_urls_lipids.tsv', sep="\t", index=False, )

# Creating enzyme2reaction for lipid reactions
enzyme2reaction = pd.read_csv("network_files/proteins_genes_unmerged.tsv", sep="\t")
enzyme2reaction = enzyme2reaction[["Rhea_ID", "Enzyme_ID"]]
enzyme2reaction.columns = ["reaction", "enzyme"]
enzyme2reaction = enzyme2reaction.drop_duplicates()
enzyme2reaction = enzyme2reaction.reset_index(drop=True)
enzyme2reaction = enzyme2reaction[enzyme2reaction.reaction.isin(lipid_reactions)]
enzyme2reaction = enzyme2reaction.reset_index(drop=True)
enzyme2reaction.to_csv(r'network_files/enzyme2reaction_lipids.tsv',
                       sep="\t", index=False, )

# # Creating genes annotation for lipid reactions
# gene2reaction = pd.read_csv("network_files/gene2reaction.tsv", sep="\t")
# gene2reaction = gene2reaction[gene2reaction.reaction.isin(lipid_reactions)]
# gene2reaction = gene2reaction.reset_index(drop=True)
# gene2reaction.to_csv(r'network_files/gene2reaction_lipids.tsv',
#                      sep="\t", index=False, )

# Creating reaction2align for lipid reactions
reaction2align = pd.read_csv("network_files/reaction2align.csv", sep=",")
reaction2align = reaction2align[reaction2align.reaction.isin(lipid_reactions)]
reaction2align = reaction2align.reset_index(drop=True)
reaction2align.to_csv(r'network_files/reaction2align_lipids.tsv', sep="\t", index=False)

# Creating full rhea metabolites & atoms files for lipidomic reactions
atoms1 = reaction2align[["atom.x"]]
atoms1.columns = ["atom"]
atoms2 = reaction2align[["atom.y"]]
atoms2.columns = ["atom"]
atoms_full = pd.concat([atoms1, atoms2])
atoms_full["metabolite"], atoms_full["coords"] = atoms_full["atom"].str.split("_", 1).str
lipid_mets = atoms_full[["metabolite"]]
lipid_mets = lipid_mets.drop_duplicates()
lipid_mets = lipid_mets.reset_index(drop=True)

atoms = pd.read_csv("network_files/atoms.csv", sep=",")
atoms = atoms[atoms.metabolite.isin(lipid_mets.metabolite.values)]
atoms = atoms.reset_index(drop=True)
atoms.to_csv(r'network_files/atoms_lipids.tsv', sep="\t", index=False)

metabolites = pd.read_csv("network_files/annotation_of_nodes_master_graph.csv", sep=";")
metabolites = metabolites[["metabolite", "metabolite_name", "metabolite_url"]]
metabolites = metabolites[metabolites.metabolite.isin(lipid_mets.metabolite.values)]
metabolites = metabolites.reset_index(drop=True)
metabolites.to_csv(r'network_files/metabolites_lipids.tsv', sep="\t", index=False)

# Creating mapFrom for LipidMaps db
mapFromLipidMaps = metabolites[["metabolite", "LipidMaps_ID"]]
mapFromLipidMaps = mapFromLipidMaps[mapFromLipidMaps.LipidMaps_ID != "-"]
mapFromLipidMaps = mapFromLipidMaps.drop_duplicates()
mapFromLipidMaps = mapFromLipidMaps.reset_index(drop=True)
mapFromLipidMaps.columns = ["metabolite", "LipidMaps"]
mapFromLipidMaps.to_csv(r'network_files/mapFromLipidMaps.tsv', sep="\t", index=False)

# Creating mapFrom for SwissLipids db
mapFromSwissLipids = metabolites[["metabolite", "SwissLipids_ID"]]
mapFromSwissLipids = mapFromSwissLipids[mapFromSwissLipids.SwissLipids_ID != "-"]
mapFromSwissLipids = mapFromSwissLipids.drop_duplicates()
mapFromSwissLipids = mapFromSwissLipids.reset_index(drop=True)
mapFromSwissLipids.columns = ["metabolite", "SwissLipids"]
mapFromSwissLipids.to_csv(r'network_files/mapFromSwissLipids_ver1.tsv', sep="\t", index=False)

