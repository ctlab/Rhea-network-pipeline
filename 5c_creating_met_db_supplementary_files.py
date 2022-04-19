import pandas as pd


# Creating hmdb.csv
hmdb = pd.DataFrame(columns=["HMDB", "metabolite"])
hmdb.to_csv(r'network_files/hmdb.csv', sep=";", index=False, )

# Creating anomers_suppl.csv
anomers = pd.read_csv("chebi_pH7_3_mapping.tsv", sep="\t")
anomers = anomers[["CHEBI", "CHEBI_PH7_3"]]
anomers.columns = ["metabolite", "base_metabolite"]
anomers.metabolite = "CHEBI:" + anomers.metabolite.astype(int).astype(str)
anomers.base_metabolite = "CHEBI:" + anomers.base_metabolite.astype(int).astype(str)

anomers_1 = anomers[["base_metabolite", "metabolite"]]
anomers_2 = anomers[["metabolite", "base_metabolite"]]
anomers_2.columns = ["metabolite", "base_metabolite"]
anomers = pd.concat([anomers_1, anomers_2])
anomers = anomers.drop_duplicates()
anomers.to_csv(r'network_files/anomers_suppl.csv', sep=";", index=False)

# Creating metabolite2atom
mets2atom = pd.read_csv("network_files/metabolite2atom.csv", sep=",")
mets2atom = mets2atom[["metabolite"]].drop_duplicates()
mets2atom = mets2atom.reset_index(drop=True)


# Reading data
chemical_data = pd.read_csv("data/chemical_data.tsv", sep="\t")
names = pd.read_csv("data/names.tsv", sep="\t")

lm_che_swl = pd.read_csv("data/lipidmaps_to_chebi_to_swisslipids.tsv", sep="\t")
lipid_classes_data = pd.read_csv("data/lipidmaps_classes.tsv", sep="\t")
swisslipids_classes = pd.read_csv("data/swisslipids_classes.tsv", sep="\t")


#
# Creating annotation for graph nodes
#

# Preprocessing chemical data
chemical_data.loc[chemical_data["SOURCE"] == "ChEBI", "SOURCE"] = "!ChEBI"
chemical_data = chemical_data.sort_values(['COMPOUND_ID', 'SOURCE'], ascending=[True, True])
chemical_data = chemical_data.drop_duplicates(subset=['COMPOUND_ID', 'TYPE'], keep="first")

chemical_data = chemical_data[["COMPOUND_ID", "TYPE", "CHEMICAL_DATA"]]
chemical_data.columns = ["ChEBI_ID", "Type", "Chemical_Data"]
chemical_data = chemical_data.drop_duplicates()
chemical_data = chemical_data.sort_values(by=['ChEBI_ID'])
chemical_data = chemical_data.reset_index(drop=True)
chemical_data["ChEBI_ID"] = "CHEBI:" + chemical_data["ChEBI_ID"].astype(str)

chemical_data = chemical_data.pivot(index="ChEBI_ID", columns="Type", values="Chemical_Data")
chemical_data = chemical_data.reset_index()
chemical_data.columns = ["ChEBI_ID", "Charge", "Formula", "Mass", "Monoisotopic_mass"]
chemical_data.loc[chemical_data["Charge"].isnull(), "Charge"] = "-"
chemical_data.loc[chemical_data["Formula"].isnull(), "Formula"] = "-"
chemical_data.loc[chemical_data["Mass"].isnull(), "Mass"] = "-"
chemical_data.loc[chemical_data["Monoisotopic_mass"].isnull(), "Monoisotopic_mass"] = "-"


# Preprocessing standard names
names["COMPOUND_ID"] = "CHEBI:" + names["COMPOUND_ID"].astype(str)

s = names.NAME.str.len().sort_values().index
names = names.reindex(s)
names = names.reset_index(drop=True)

names.loc[names["TYPE"] == "BRAND NAME", "TYPE"] = "!BRAND NAME"
names.loc[names["TYPE"] == "NAME", "TYPE"] = "!NAME"
names.loc[names["TYPE"] == "IUPAC NAME", "TYPE"] = "ZIUPAC NAME"

names = names.sort_values(['COMPOUND_ID', 'TYPE', 'SOURCE', 'ADAPTED'], ascending=[True, True, True, False])
names = names.drop_duplicates(subset=['COMPOUND_ID'], keep="first")

names = names[["COMPOUND_ID", "NAME"]]
names.columns = ["ChEBI_ID", "Name"]
names = names.drop_duplicates()
names = names.sort_values(by=['ChEBI_ID'])
names = names.reset_index(drop=True)


# Merging ChemData with Standard names
chebi_data = pd.merge(names, chemical_data, how="outer",
                      left_on="ChEBI_ID", right_on="ChEBI_ID")

chebi_data.loc[chebi_data["Name"].isnull(), "Name"] = "-"
chebi_data.loc[chebi_data["Charge"].isnull(), "Charge"] = "-"
chebi_data.loc[chebi_data["Formula"].isnull(), "Formula"] = "-"
chebi_data.loc[chebi_data["Mass"].isnull(), "Mass"] = "-"
chebi_data.loc[chebi_data["Monoisotopic_mass"].isnull(), "Monoisotopic_mass"] = "-"
chebi_data["Link"] = "http://www.ebi.ac.uk/chebi/searchId.do?chebiId=" + chebi_data["ChEBI_ID"].astype(str)


chebi_data = chebi_data[["ChEBI_ID", "Name", "Link", "Charge", "Mass", "Formula"]]
chebi_data.columns = ["metabolite", "metabolite_name", "metabolite_url",
                      "metabolite_charge", "metabolite_mass", "metabolite_formula"]


# Preprocessing lm_che_pub_swlp
lm_che_swl.columns = ["metabolite", "LipidMaps_ID", "SwissLipids_ID"]


# Merging
gna = pd.merge(mets2atom, chebi_data, how="left")
gna = pd.merge(gna, lm_che_swl, how="left")
gna.loc[gna["LipidMaps_ID"].isnull(), "LipidMaps_ID"] = "-"
gna.loc[gna["SwissLipids_ID"].isnull(), "SwissLipids_ID"] = "-"

gna = pd.merge(gna, lipid_classes_data, how="left")
gna = pd.merge(gna, swisslipids_classes, how="left")

gna.loc[gna["Class_LipidMaps"].isnull(), "Class_LipidMaps"] = "-"
gna.loc[gna["Abbreviation_LipidMaps"].isnull(), "Abbreviation_LipidMaps"] = "-"
gna.loc[gna["Synonyms_LipidMaps"].isnull(), "Synonyms_LipidMaps"] = "-"

gna.loc[gna["Class_SwissLipids"].isnull(), "Class_SwissLipids"] = "-"
gna.loc[gna["Abbreviation_SwissLipids"].isnull(), "Abbreviation_SwissLipids"] = "-"

gna = gna[["metabolite", "metabolite_name", "metabolite_url", "metabolite_charge",
           "metabolite_mass", "metabolite_formula",
           "LipidMaps_ID", "SwissLipids_ID",
           "Class_LipidMaps", "Abbreviation_LipidMaps", "Synonyms_LipidMaps",
           "Class_SwissLipids", "Abbreviation_SwissLipids"]]


'''As some of the metabolites have same ChEBI ID for different LipidMaps IDs we need to fix repetitive annotation.'''


def fixing(x):
    y = list(set(x))
    if ("-" in y) & (len(y) > 1):
        y = list(set(y) - {"-"})
    x = " // ".join(y)
    return x


gna_1 = gna.groupby(['metabolite', 'metabolite_name', 'metabolite_url', 'metabolite_charge', 'metabolite_mass',
                     'metabolite_formula']).agg({'LipidMaps_ID': lambda x: fixing(x),
                                                 'SwissLipids_ID': lambda x: fixing(x),
                                                 'Class_LipidMaps': lambda x: fixing(x),
                                                 'Abbreviation_LipidMaps': lambda x: fixing(x),
                                                 'Synonyms_LipidMaps': lambda x: fixing(x),
                                                 'Class_SwissLipids': lambda x: fixing(x),
                                                 'Abbreviation_SwissLipids': lambda x: fixing(x)})


gna_1 = gna_1.reset_index()

# Saving result dataframe
gna_1.to_csv(r'network_files/annotation_of_nodes_master_graph.csv', sep=";", index=False, )

# Creating mapFrom for SwissLipids
mapFromSwissLipids = gna_1[["metabolite", "SwissLipids_ID"]]
mapFromSwissLipids = mapFromSwissLipids.drop_duplicates()
mapFromSwissLipids = mapFromSwissLipids.reset_index(drop=True)
mapFromSwissLipids.to_csv(r'network_files/mapFromSwissLipids.csv', sep=";", index=False, )

