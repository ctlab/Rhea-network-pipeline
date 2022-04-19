import pandas as pd


undirected_to_LR = pd.read_csv("data/rhea-directions.tsv", sep="\t")


# For test run
# undirected_to_LR = undirected_to_LR[0:5]


LR_Rhea_IDs = list(undirected_to_LR.RHEA_ID_LR.values)

for i in range(len(LR_Rhea_IDs)):
    LR_Rhea_IDs[i] = str(LR_Rhea_IDs[i]) + ".rxn"


with open("data/LR_Rhea_mapping_rxn.txt", "x") as file:
    print(*LR_Rhea_IDs, file=file, sep="\n")


undirected_to_LR = undirected_to_LR[["RHEA_ID_MASTER", "RHEA_ID_LR"]]
undirected_to_LR.columns = ["Rhea_ID", "LR_Rhea_ID"]

undirected_to_LR.to_csv(r'data/undirected_to_LR.tsv', sep="\t", index=False, )
