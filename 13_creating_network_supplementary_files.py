import re
import time
import requests
import pandas as pd
from bioservices import *
ch = ChEBI()

#
# Getting genes from enzymes
#

rhea2ec = pd.read_csv("data/rhea2ec.tsv", sep="\t")
rhea2ec = rhea2ec[["MASTER_ID", "ID"]]
rhea2ec.columns = ["Rhea_ID", "Enzyme_ID"]

master = pd.read_csv("data/rhea-directions.tsv", sep="\t")
master = master[["RHEA_ID_MASTER"]]
master.columns = ["Rhea_ID"]

rhea_enz = pd.merge(master, rhea2ec, how="left")
rhea_enz.loc[rhea_enz['Enzyme_ID'].isnull(), 'Enzyme_ID'] = "-"
rhea_enz.columns = ["reaction", "enzyme"]
rhea_enz.to_csv(r'network_files/pre_enzyme2reaction.csv', sep=";", index=False, )

# Adding KEGG enzymes to Rhea reactions
rhea2enz = pd.read_csv("network_files/pre_enzyme2reaction.csv", sep=";")
rhea2enz.columns = ["rhea", "enzyme"]

kegg2enz = pd.read_csv("data/react2enz_kegg.csv")
kegg2enz.columns = ["kegg", "enzyme"]

kegg2rhea = pd.read_csv("data/rhea2kegg_reaction.tsv", sep="\t")
kegg2rhea = kegg2rhea[["MASTER_ID", "ID"]]
kegg2rhea.columns = ["rhea", "kegg"]

rhea2enz_kegg = pd.merge(kegg2rhea, kegg2enz, how="left")
rhea2enz_kegg = rhea2enz_kegg[["rhea", "enzyme"]]

rhea2enz_prep = pd.concat([rhea2enz, rhea2enz_kegg])
rhea2enz_prep = rhea2enz_prep.drop_duplicates()
rhea2enz_prep = rhea2enz_prep[rhea2enz_prep.enzyme != "-"]
rhea2enz_prep = rhea2enz_prep.reset_index(drop=True)

reactions = rhea2enz[["rhea"]]
reactions = reactions.drop_duplicates()
reactions = reactions.reset_index(drop=True)

rhea2enz_full = pd.merge(reactions, rhea2enz_prep, how="left")
rhea2enz_full.loc[rhea2enz_full.enzyme.isnull(), "enzyme"] = "-"
rhea2enz_full.columns = ["reaction", "enzyme"]

rhea2enz_full.to_csv(r'network_files/enzyme2reaction.csv', sep=";", index=False, )

rhea_enz = pd.read_csv("network_files/enzyme2reaction.csv", sep=";")
rhea_enz.columns = ["Rhea_ID", "Enzyme_ID"]

kegg_hsa = pd.read_csv("genes_data/hsa.tsv", sep="\t", header=None)
kegg_hsa.columns = ["Gene_ID", "Enzyme_ID"]

kegg_hsa['Gene_ID'] = kegg_hsa['Gene_ID'].str.replace('hsa:', '')
kegg_hsa['Enzyme_ID'] = kegg_hsa['Enzyme_ID'].str.replace('ec:', '')


# Adding genes

annotation = pd.merge(rhea_enz, kegg_hsa, how="left")
annotation.loc[annotation["Gene_ID"].isnull(), "Gene_ID"] = "-"
annotation["Rhea_link"] = "http://www.ebi.ac.uk/rhea/reaction.xhtml?id=" + annotation["Rhea_ID"].astype(str)
annotation = annotation[["Rhea_ID", "Rhea_link", "Enzyme_ID", "Gene_ID"]]


annotation = annotation.drop_duplicates()
annotation = annotation.reset_index(drop=True)

annotation.to_csv(r'data/reactions_annotation.tsv',
                  sep="\t", index=False, )

annotation2 = annotation[["Rhea_ID", "Rhea_link"]]
annotation2.columns = ["reaction", "reaction_url"]
annotation2 = annotation2.drop_duplicates()
annotation2 = annotation2.reset_index(drop=True)
annotation2.to_csv(r'network_files/reaction_urls.tsv', sep="\t", index=False, )

annotation = annotation[["Rhea_ID", "Enzyme_ID", "Gene_ID"]]
annotation = annotation.drop_duplicates()
annotation = annotation.reset_index(drop=True)
annotation.to_csv(r'data/reactions_enzymes_genes_rhea_network_unmerged.tsv',
                  sep="\t", index=False)


#
# Getting human proteins
#

human_proteins = pd.read_csv("genes_data/hsa.proteins.tsv", sep="\t")
human_proteins = human_proteins[["Entry", "Cross-reference (GeneID)"]]
human_proteins.columns = ["Protein_ID", "Gene_ID_protein"]

human_proteins.loc[human_proteins.Gene_ID_protein.isnull(), "Gene_ID_protein"] = "-"
human_proteins.Gene_ID_protein = human_proteins.Gene_ID_protein.str.replace(r'(;)$', '')
human_proteins = human_proteins.reset_index(drop=True)

prot_ids_list = list(human_proteins.Protein_ID)
i = -1
gene_id_l = []
indexes = []

for gene_id in human_proteins.Gene_ID_protein:
    i += 1

    gene_id = str(gene_id)
    gene_id_l = list(gene_id)

    if ";" in gene_id_l:
        indexes.append(i)

count = 0
if indexes != []:
    for i in indexes:
        count += 1
        if count == 10:
            break

if indexes != []:
    for i in indexes:
        prot_id = human_proteins.Protein_ID[i]
        gene_ids_list = human_proteins.Gene_ID_protein[i]
        gene_ids_list = gene_ids_list.split(';')

        for gene_id in gene_ids_list:
            human_proteins = human_proteins.append({"Protein_ID": prot_id,
                                                    "Gene_ID_protein": gene_id},
                                                   ignore_index=True)

human_proteins = human_proteins.drop(indexes)
human_proteins = human_proteins.reset_index(drop=True)

print("Dataframe length =", len(human_proteins))
human_proteins.to_csv(r'data/proteins_genes_unmerged.csv', index=False)

#
#
# R GENES ANNOTATION INSERT HERE
#
#

#
# Merging enzyme genes with protein genes
#

human_proteins = pd.read_csv("data/proteins_genes_unmerged.csv", sep=",")

rhea2uniprot = pd.read_csv("data/rhea2uniprot_sprot.tsv", sep="\t")
rhea2uniprot = rhea2uniprot[["MASTER_ID", "ID"]]
rhea2uniprot.columns = ["Rhea_ID", "Protein_ID"]


proteins = pd.merge(rhea2uniprot, human_proteins, how="inner")
proteins = proteins[["Rhea_ID", "Protein_ID", "Gene_ID_protein"]]
proteins.columns = ["Rhea_ID", "Protein_ID", "Gene_ID"]
proteins = proteins[proteins.Gene_ID != "-"]

proteins = proteins.drop_duplicates()
proteins = proteins.reset_index(drop=True)


enzymes_genes = pd.read_csv("data/reactions_enzymes_genes_rhea_network_unmerged.tsv", sep="\t")
enzymes_genes = enzymes_genes[["Rhea_ID", "Enzyme_ID", "Gene_ID"]]
enzymes_genes = enzymes_genes.drop_duplicates()
enzymes_genes = enzymes_genes.reset_index(drop=True)


proteins = proteins[["Rhea_ID", "Gene_ID"]]
proteins.Gene_ID = proteins.Gene_ID.astype(int).astype(str)


proteins = pd.merge(proteins, enzymes_genes,
                    how="left",
                    left_on=["Rhea_ID", "Gene_ID"],
                    right_on=["Rhea_ID", "Gene_ID"])

proteins.loc[proteins["Enzyme_ID"].isnull(), "Enzyme_ID"] = "-"
proteins = proteins[["Rhea_ID", "Enzyme_ID", "Gene_ID"]]


final_genes = pd.concat([proteins, enzymes_genes])
final_genes = final_genes[(final_genes.Enzyme_ID != "-") | (final_genes.Gene_ID != "-")]
final_genes = final_genes.drop_duplicates()
final_genes = final_genes.reset_index(drop=True)


enzymes_genes = enzymes_genes[["Rhea_ID"]]
enzymes_genes = enzymes_genes.drop_duplicates()
enzymes_genes = enzymes_genes.reset_index(drop=True)


final_genes = pd.merge(enzymes_genes, final_genes, how="left")
final_genes.loc[final_genes["Enzyme_ID"].isnull(), "Enzyme_ID"] = "-"
final_genes.loc[final_genes["Gene_ID"].isnull(), "Gene_ID"] = "-"


final_genes = final_genes[["Rhea_ID", "Gene_ID"]]
final_genes.columns = ["reaction", "gene"]
final_genes = final_genes.drop_duplicates()
final_genes = final_genes.reset_index(drop=True)
final_genes.to_csv(r'network_files/gene2reaction.tsv', sep="\t", index=False, )
