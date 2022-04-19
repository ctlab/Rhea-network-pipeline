import pandas as pd
import numpy as np


def translation(db_from, db_to):
    from_ID = ""
    to_ID = ""

    if db_from == "LipidMaps":
        file = "data/databases/structures.txt"
        line_from_ID_start = "> <LM_ID>\n"

        if db_to == "ChEBI":
            line_to_ID_start = '> <CHEBI_ID>\n'

        elif db_to == "SwissLipids":
            line_to_ID_start = '> <SWISSLIPIDS_ID>\n'

    elif db_from == "ChEBI":
        file = "data/databases/ChEBI_complete_3star.txt"
        line_from_ID_start = "> <ChEBI ID>\n"

        if db_to == "LipidMaps":
            line_to_ID_start = '> <LIPID MAPS instance Database Links>\n'

        elif db_to == "SwissLipids":
            line_to_ID_start = '> <SwissLipids Database Links>\n'

    flag_1 = False
    flag_2 = False
    flag_3 = False
    flag_4 = False

    annotation = pd.DataFrame(data=[], columns=[(db_from + "_ID"), (db_to + "_ID")])

    handle = open(file, "r")

    for line in handle:

        if flag_1:
            from_ID = line[:-1]
            flag_1 = False
            flag_3 = True

        if flag_2 and flag_3:
            if db_to == "ChEBI":
                to_ID = "CHEBI:" + line[:-1]
            else:
                to_ID = line[:-1]
            flag_2 = False
            flag_4 = True

        if to_ID != "" and from_ID != "" and flag_4:
            annotation = annotation.append({db_from + "_ID": from_ID,
                                            db_to + "_ID": to_ID},
                                           ignore_index=True)
            flag_1 = False
            flag_2 = False
            flag_3 = False
            flag_4 = False

        if line == line_from_ID_start:
            flag_1 = True

        if line == "$$$$\n":
            flag_1 = False
            flag_2 = False
            flag_3 = False

        if flag_3 & (line == line_to_ID_start):
            flag_2 = True

    handle.close()

    return annotation


lipidmaps_to_chebi = translation(db_from="LipidMaps", db_to="ChEBI")
lipidmaps_to_swisslipids = translation(db_from="LipidMaps", db_to="SwissLipids")
chebi_to_lipidmaps = translation(db_from="ChEBI", db_to="LipidMaps")
chebi_to_swisslipids = translation(db_from="ChEBI", db_to="SwissLipids")

lipids = pd.read_csv("data/databases/lipids.tsv", sep="\t", encoding="latin-1")
swisslipids_to_chebi = lipids[["Lipid ID", "CHEBI"]]
swisslipids_to_chebi.columns = ["SwissLipids_ID", "ChEBI_ID"]
swisslipids_to_chebi = swisslipids_to_chebi[swisslipids_to_chebi.ChEBI_ID.notnull()]
swisslipids_to_chebi = swisslipids_to_chebi.drop_duplicates()
swisslipids_to_chebi = swisslipids_to_chebi.reset_index(drop=True)
swisslipids_to_chebi.ChEBI_ID = swisslipids_to_chebi.ChEBI_ID.astype(str)
bad = swisslipids_to_chebi[swisslipids_to_chebi.ChEBI_ID.str.contains(' | ')]
bad = bad.reset_index(drop=True)
for i in range(len(bad.SwissLipids_ID)):
    lipid_id = bad.SwissLipids_ID[i]
    chebi_ids_list = bad.ChEBI_ID[i]
    chebi_ids_list = chebi_ids_list.split(' | ')
    for chebi_id in chebi_ids_list:
        swisslipids_to_chebi = swisslipids_to_chebi.append({"SwissLipids_ID": lipid_id,
                                                            "ChEBI_ID": chebi_id}, ignore_index=True)
swisslipids_to_chebi = swisslipids_to_chebi[~swisslipids_to_chebi.ChEBI_ID.str.contains(' | ')]
swisslipids_to_chebi.ChEBI_ID = "CHEBI:" + swisslipids_to_chebi.ChEBI_ID.astype(str)
swisslipids_to_chebi = swisslipids_to_chebi.drop_duplicates()
swisslipids_to_chebi = swisslipids_to_chebi.reset_index(drop=True)

swisslipids_to_lipidmaps = lipids[["Lipid ID", "LIPID MAPS"]]
swisslipids_to_lipidmaps.columns = ["SwissLipids_ID", "LipidMaps_ID"]
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps[swisslipids_to_lipidmaps.LipidMaps_ID.notnull()]
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps.drop_duplicates()
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps.reset_index(drop=True)
swisslipids_to_lipidmaps.LipidMaps_ID = swisslipids_to_lipidmaps.LipidMaps_ID.astype(str)
bad = swisslipids_to_lipidmaps[swisslipids_to_lipidmaps.LipidMaps_ID.str.contains(' | ')]
bad = bad.reset_index(drop=True)
for i in range(len(bad.SwissLipids_ID)):
    lipid_id = bad.SwissLipids_ID[i]
    lm_ids_list = bad.LipidMaps_ID[i]
    lm_ids_list = lm_ids_list.split(' | ')
    for lm_id in lm_ids_list:
        swisslipids_to_lipidmaps = swisslipids_to_lipidmaps.append({"SwissLipids_ID": lipid_id,
                                                                    "LipidMaps_ID": lm_id}, ignore_index=True)
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps[~swisslipids_to_lipidmaps.LipidMaps_ID.str.contains(' | ')]
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps.drop_duplicates()
swisslipids_to_lipidmaps = swisslipids_to_lipidmaps.reset_index(drop=True)
swisslipids_to_lipidmaps.columns = ["SwissLipids_ID", "LipidMaps_ID"]


# Merging
chebi_and_lipidmaps = pd.concat([chebi_to_lipidmaps[["ChEBI_ID", "LipidMaps_ID"]],
                                 lipidmaps_to_chebi[["ChEBI_ID", "LipidMaps_ID"]]])
chebi_and_lipidmaps = chebi_and_lipidmaps.drop_duplicates()
chebi_and_lipidmaps = chebi_and_lipidmaps.reset_index(drop=True)

chebi_and_swisslipids = pd.concat([chebi_to_swisslipids[["ChEBI_ID", "SwissLipids_ID"]],
                                   swisslipids_to_chebi[["ChEBI_ID", "SwissLipids_ID"]]])
chebi_and_swisslipids = chebi_and_swisslipids.drop_duplicates()
chebi_and_swisslipids = chebi_and_swisslipids.reset_index(drop=True)

swisslipids_and_lipidmaps = pd.concat([swisslipids_to_lipidmaps[["LipidMaps_ID", "SwissLipids_ID"]],
                                       lipidmaps_to_swisslipids[["LipidMaps_ID", "SwissLipids_ID"]]])
swisslipids_and_lipidmaps = swisslipids_and_lipidmaps.drop_duplicates()
swisslipids_and_lipidmaps = swisslipids_and_lipidmaps.reset_index(drop=True)

ch_lm_swl = pd.merge(left=chebi_and_lipidmaps, right=chebi_and_swisslipids,
                     left_on='ChEBI_ID', right_on='ChEBI_ID', how='outer')
ch_lm_swl = ch_lm_swl.replace(np.nan, '-', regex=True)

for i in range(len(ch_lm_swl)):
    swl_id = ch_lm_swl.SwissLipids_ID[i]
    if swl_id == "-":
        try:
            swl_id = swisslipids_and_lipidmaps[swisslipids_and_lipidmaps["LipidMaps_ID"] ==
                                               ch_lm_swl.LipidMaps_ID[i]].SwissLipids_ID.iloc[0]
        except IndexError:
            swl_id = "-"
    ch_lm_swl.SwissLipids_ID[i] = swl_id

for i in range(len(ch_lm_swl)):
    lm_id = ch_lm_swl.LipidMaps_ID[i]
    if lm_id == "-":
        try:
            lm_id = swisslipids_and_lipidmaps[swisslipids_and_lipidmaps["SwissLipids_ID"] ==
                                              ch_lm_swl.SwissLipids_ID[i]].LipidMaps_ID.iloc[0]
        except IndexError:
            lm_id = "-"
    ch_lm_swl.LipidMaps_ID[i] = lm_id

ch_lm_swl = ch_lm_swl.drop_duplicates()
ch_lm_swl = ch_lm_swl.reset_index(drop=True)

ch_lm_swl.to_csv(r'data/lipidmaps_to_chebi_to_swisslipids.tsv', index=False, sep="\t")

#
# # Dealing with classes
#

# LipidMaps

def lipidmaps_ids_with_classes():
    lipidmaps_ID = ""
    lm_name = ""
    lm_class = ""
    abbrev = ""
    synonyms = ""
    flag_lm_id_line = False
    flag_found_lm_id = False
    flag_name_line = False
    flag_found_name = False
    flag_class_line = False
    flag_found_class = False
    flag_abbrv_line = False
    flag_found_abbrv = False
    flag_syn_line = False
    flag_syn_found = False

    annotation = pd.DataFrame(data=[],
                              columns=["LipidMaps_ID", "Name_LipidMaps", "Class_LipidMaps",
                                       "Abbreviation_LipidMaps", "Synonyms_LipidMaps"])
    file = open("data/databases/structures.txt", "r")

    for line in file:

        # Adding LipidMaps ID
        if flag_lm_id_line:
            lipidmaps_ID = line[:-1]
            flag_lm_id_line = False
            flag_found_lm_id = True

        # Adding LipidMaps name
        if flag_name_line and flag_found_lm_id:
            lm_name = line[:-1]
            flag_name_line = False
            flag_found_name = True

        # Adding lipid class
        if flag_class_line and flag_found_lm_id:
            lm_class = line[:-1]
            flag_class_line = False
            flag_found_class = True

        # Adding abbreviation
        if (flag_class_line or flag_found_class) and \
                (flag_abbrv_line) and (flag_found_lm_id):
            abbrev = line[:-1]
            flag_abbrv_line = False
            flag_found_abbrv = True

        # Adding synonyms
        if (flag_class_line or flag_found_class) and \
                (flag_syn_line) and (flag_found_abbrv):
            synonyms = line[:-1]
            flag_syn_line = False
            flag_syn_found = True

        if line == "> <LM_ID>\n":
            flag_lm_id_line = True

        if line == "$$$$\n":

            if lipidmaps_ID != "" and \
                    (flag_found_name or flag_found_class or \
                     flag_found_abbrv or flag_syn_found):

                if lm_name == "":
                    lm_name = "-"
                if lm_class == "":
                    lm_class = "-"
                if abbrev == "":
                    abbrev = "-"
                if synonyms == "":
                    synonyms = "-"

                annotation = annotation.append({"LipidMaps_ID": lipidmaps_ID,
                                                "Name_LipidMaps": lm_name,
                                                "Class_LipidMaps": lm_class,
                                                "Abbreviation_LipidMaps": abbrev,
                                                "Synonyms_LipidMaps": synonyms},
                                               ignore_index=True)
            # Resetting flags
            flag_lm_id_line = False
            flag_found_lm_id = False
            flag_name_line = False
            flag_found_name = False
            flag_class_line = False
            flag_found_class = False
            flag_abbrv_line = False
            flag_found_abbrv = False
            flag_syn_line = False
            flag_syn_found = False

            lipidmaps_ID = ""
            lm_name = ""
            lm_class = ""
            abbrev = ""
            synonyms = ""

        # Next line will contain necessary data
        if (flag_found_lm_id == True) & (line == '> <NAME>\n'):
            flag_name_line = True
        elif (flag_found_lm_id == True) & (line == '> <CATEGORY>\n'):
            flag_class_line = True
        elif (flag_found_lm_id == True) & (line == '> <ABBREVIATION>\n'):
            flag_abbrv_line = True
        elif (flag_found_lm_id == True) & (line == '> <SYNONYMS>\n'):
            flag_syn_line = True

    file.close()
    return annotation

lipidmaps_annotation = lipidmaps_ids_with_classes()
lipidmaps_annotation.to_csv(r'data/lipidmaps_classes.tsv', index=False, sep="\t")

# SwissLipids

lipids = pd.read_csv("data/databases/lipids.tsv", sep="\t", encoding="latin-1")
swisslipids_classes = lipids[["Lipid ID", "Lipid class*", "Abbreviation*"]]
classes = lipids[["Lipid class*"]]
classes = classes.drop_duplicates()
classes = classes.reset_index(drop=True)
classes["Class_name"] = ""
for i in range(len(classes)):
    try:
        classes.Class_name[i] = lipids.Name[lipids["Lipid ID"] == classes["Lipid class*"][i]].iloc[0]
    except IndexError:
        classes.Class_name[i] = "-"
classes.columns = ["Lipid_class", "Class_name"]
swisslipids_classes = pd.merge(left=swisslipids_classes, right=classes,
                            left_on='Lipid class*', right_on='Lipid_class', how='left')
swisslipids_classes = swisslipids_classes[["Lipid ID", "Class_name", "Abbreviation*"]]
swisslipids_classes.columns = ["Lipid_ID", "Lipid_class", "Abbreviation"]
swisslipids_classes = swisslipids_classes.replace(np.nan, '-', regex=True)
swisslipids_classes = swisslipids_classes.drop_duplicates()
swisslipids_classes = swisslipids_classes.reset_index(drop=True)
swisslipids_classes = pd.DataFrame(swisslipids_classes)
swisslipids_classes.columns = ["SwissLipids_ID", "Class_SwissLipids", "Abbreviation_SwissLipids"]

swisslipids_classes.to_csv(r'data/swisslipids_classes.tsv', index=False, sep="\t")
