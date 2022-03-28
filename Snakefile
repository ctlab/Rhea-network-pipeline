import os


rule all:
    input:
        "network/network.rhea.rda",
        "network/network.rhea.rds",
        "network/met.rhea.db.rda",
        "network/met.rhea.db.rds",
        "network/network.rhea.lipids.rda",
        "network/network.rhea.lipids.rds",
        "network/met.lipids.db.rda",
        "network/met.lipids.db.rds"


## 1. Getting reactions
rule get_reactions_and_enzymes:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    output:
        "data/rhea2ec.tsv"
    message:
        "...1. Downloading reactions' list..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv -O {output}"


## 2. Downloading all rxn files
rule download_rxn_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv"
    output:
        "data/rhea-rxn.tar.gz"
    message:
        "...2. Downloading rxn files..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/ctfiles/rhea%2Drxn.tar.gz -O {output}"


## 3. Downloading UNDIRECTED to LR to RL
rule download_rxn_directions:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea-rxn.tar.gz"
    output:
        "data/rhea-directions.tsv"
    message:
        "...3. Downloading reactions' directions mapping list..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv -O {output}"


## 4. Getting LR rxn list
rule sorting_rxn_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea-directions.tsv",
        "data/rhea2ec.tsv",

        "4_getting_LR_rxns_list.py"
    output:
        "data/LR_Rhea_mapping_rxn.txt",
        "data/undirected_to_LR.tsv"
    message:
        "...4. Creating LR reactions' list..."
    shell:
        "python3.9 4_getting_LR_rxns_list.py"


rule pre_RDT:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/LR_Rhea_mapping_rxn.txt"
    output:
        "data/reports/rdt_test_done.txt"
    message:
        "...5. Preparing for RDT step..."
    shell:
        "cwd=$(pwd);"
        "mkdir -p data/pre_rdt_output;"
        "cd data/pre_rdt_output;"
        "java -jar -Xmx8G /rdt-2.4.1-jar-with-dependencies.jar -Q SMI -q 'CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O' -g -c -j AAM -f TEXT || true;"
        "touch $cwd/data/reports/rdt_test_done.txt"


## 5. Moving LR to another folder
checkpoint moving_rxns:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/LR_Rhea_mapping_rxn.txt",
        "data/undirected_to_LR.tsv",
        "data/reports/rdt_test_done.txt"
    output:
        reactions=directory("data/LR_rxn")
    message:
        "...5. Moving LR reactions..."
    shell:
        'mkdir -p data/LR_rxn;'
        'for file in `cat data/LR_Rhea_mapping_rxn.txt`; do tar -xzf data/rhea-rxn.tar.gz rxn/$file || true; done;'
        "cp ./rxn/*.rxn ./data/LR_rxn/"


## 6. Running RDT
rule performing_RDT:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        lr_rxn="data/LR_rxn/{reaction}.rxn",
        rxn="rxn/{reaction}.rxn"
    output:
        "rdt_output/{reaction}/ECBLAST_{reaction}_AAM.rxn"
    message:
        "...7. Performing RDT analysis..."
    shell:
        "cwd=$(pwd);"
        "mkdir -p rdt_output failed_rdt;"
        "cd rdt_output/{wildcards.reaction};"
        "java -jar -Xmx8G /rdt-2.4.1-jar-with-dependencies.jar -Q RXN -q $cwd/{input.rxn} -g -j AAM -f TEXT || true;"
        "if [[ ! -f $cwd/rdt_output/{wildcards.reaction}/ECBLAST_{wildcards.reaction}_AAM.rxn ]]; then touch $cwd/failed_rdt/ECBLAST_{wildcards.reaction}_AAM.rxn; fi;"
        "touch $cwd/rdt_output/{wildcards.reaction}/ECBLAST_{wildcards.reaction}_AAM.rxn"


def aggregate_rxns(wildcards):
    checkpoint_output = checkpoints.moving_rxns.get(**wildcards).output[0]
    file_names = expand("rdt_output/{reaction}/ECBLAST_{reaction}_AAM.rxn",
        reaction=glob_wildcards(os.path.join(checkpoint_output,"{reaction}.rxn")).reaction)
    return file_names


rule aggregate:
    input:
        aggregate_rxns
    output:
        "data/rdt_done.txt"
    shell:
        "touch {output};"
        'if [ "$(ls -A failed_rdt 2> /dev/null)" != "" ]; then cwd=$(pwd); cd rdt_output;'
        'for file in $cwd/failed_rdt/*; do filename=$(basename "$file");'
        'prefolder=$(basename "$file" "_AAM.rxn"); folder=${{prefolder/ECBLAST_/}};'
        'rm $folder/$filename; done;'
        'fi'


## 7. Analysing RDT
rule performing_RDT_results_analysis:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rdt_done.txt",

        "7_Rdt_Output_Analysis.R"
    output:
        "rdt_analysis/atom_mapping.csv"
    message:
        "...7. Analysing RDT results and creating atom mapping table..."
    shell:
        "mkdir -p rdt_analysis;"
        "Rscript 7_Rdt_Output_Analysis.R"


## 8. Making RPAIRS and supplementary files
rule creating_RPAIRS_and_supplementary_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "rdt_analysis/atom_mapping.csv",
        "data/undirected_to_LR.tsv",

        "8_Making_RPAIRS_and_supplementary_files.py"
    output:
        "data/atom_mapping_C_atoms.csv",
        "data/rpairs_preprocessing.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",
        "network_files/metabolite2atom.csv",
        "network_files/rpairs.csv",
        "network_files/full_atom_mapping_1.csv"
    message:
        "...8. Creating RPAIRS and supplementary files..."
    shell:
        "mkdir -p network_files;"
        "python3.9 8_Making_RPAIRS_and_supplementary_files.py"


## 9.1. Preprocessing databases
rule preprocessing_various_annotations:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "9_1_Preprocessing_various_annotations.py"
    output:
        "data/lipidmaps_to_chebi_to_swisslipids.tsv",
        "data/swisslipids_classes.tsv",
        "data/lipidmaps_classes.tsv"
    message:
        "...9. Creating met.db.rhea supplementary files..."
    shell:
        "mkdir data/databases;"
        "wget 'https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip' -O data/databases/LM.zip;"
        "unzip data/databases/LM.zip; mv data/databases/structures.sdf data/databases/structures.txt;"
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf;"
        "mv data/databases/ChEBI_complete_3star.sdf data/databases/ChEBI_complete_3star.txt;"
        "wget 'https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv' -O data/databases/lipids.tsv;"
        "python3.9 9_1_Preprocessing_various_annotations.py"

## 9. Creating supplementary files for met.db.rhea
rule supp_files_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/rpairs.csv",
        "data/lipidmaps_to_chebi_to_swisslipids.csv",
        "data/swisslipids_classes.tsv",
        "data/lipidmaps_classes",

        "9_Creating_met_db_supplementary_files.py"
    output:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/hmdb.csv",
        "network_files/anomers_suppl.csv",
        "data/kegg2chebi.tsv",

        "9_Creating_met_db_supplementary_files.py"
    message:
        "...9. Creating met.db.rhea supplementary files..."
    shell:
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv -O data/chemical_data.tsv;"
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz -O data/names.tsv.gz;"
        "wget http://rest.kegg.jp/conv/chebi/compound -O data/kegg2chebi.tsv;"
        "wget https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv -O data/chebi_pH7_3_mapping.tsv;"
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv -O data/rhea2kegg_reaction.tsv;"
        "gzip -d data/names.tsv.gz;"
        "python3.9 9_Creating_met_db_supplementary_files.py"


## 10.1. Getting mappings for met.db
rule getting_mapping_tables_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/kegg2chebi.tsv",

        "10_1_Getting_HMDB_and KEGG_mappings.R"
    output:
        "data/HMDB2metabolite.csv",
        "data/react2enz_kegg.csv"
    message:
        "...9. Creating met.db.rhea supplementary files..."
    shell:
        "Rscript 10_1_Getting_HMDB_and KEGG_mappings.R"


## 10.2. Converting mappings for met.db
rule converting_mappings_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/kegg2chebi.tsv",
        "data/HMDB2metabolite.csv",

        "10_2_Converting_db_IDs.py"
    output:
        "network_files/mapFromHMDB.csv",
        "network_files/mapFromKEGG.csv"
    message:
        "...9. Creating met.db.rhea supplementary files..."
    shell:
        "python3.9 10_2_Converting_db_IDs.py"

## 10. Create met.db.rhea
rule creating_met_db_rhea_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/mapFromHMDB.csv",
        "network_files/mapFromKEGG.csv",
        "network_files/anomers_suppl.csv",

        "10_creating_met.db.rhea.R"
    output:
        "network/met.rhea.db.rda",
        "network/met.rhea.db.rds"
    message:
        "...10. Creating met.db.rhea object..."
    shell:
        "Rscript 10_creating_met.db.rhea.R"


## 11. Download human proteome
rule downloading_human_proteome:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network/met.rhea.db.rda"
    output:
        "genes_data/hsa.proteins.tsv"
    message:
        "...11. Downloading human proteome..."
    shell:
        "mkdir genes_data;"
        "wget 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,database(Ensembl),database(GeneID)&compress=yes' -O genes_data/hsa.proteins.tsv.gz;"
        "gzip -d genes_data/hsa.proteins.tsv.gz"


## 12. Downloading rhea to Uniprot mapping
rule download_rhea2uniprot:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network/met.rhea.db.rda"
    output:
        "data/rhea2uniprot_sprot.tsv"
    message:
        "...12. Downloading Rhea to Uniprot mapping..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv -O {output}"


## 13. Supplementary files for network
rule supp_files_for_network:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv",
        "data/rhea-directions.tsv",
        "data/rhea2uniprot_sprot.tsv",
        "genes_data/hsa.proteins.tsv",

        "13_creating_network_supplementary_files.py"
    output:
        "network_files/enzyme2reaction.csv",
        "network_files/gene2reaction.tsv",
        "network_files/reaction_urls.tsv",
        "data/proteins_genes_unmerged.tsv",
        "data/reactions_annotation.tsv",
        "genes_data/hsa.tsv"
    message:
        "...13. Creating network supplementary files..."
    shell:
        "wget http://rest.kegg.jp/link/enzyme/hsa -O genes_data/hsa.tsv;"
        "python3.9 13_creating_network_supplementary_files.py"


## 14. Creating lipid-specific graph files
rule creating_lipid_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/reaction_urls.tsv",
        "network_files/full_atom_mapping_1.csv",
        "data/proteins_genes_unmerged.tsv",
        "data/reactions_annotation.tsv",

        "14_creating_lipid_network_files.py"
    output:
        "network_files/atoms_lipids.tsv",
        "network_files/enzyme2reaction_lipids.tsv",
        "network_files/reaction_urls_lipids.tsv",
        "network_files/reaction2align_lipids.tsv",
        "network_files/mapFromLipidMaps.tsv",
        "network_files/mapFromSwissLipids_ver1.tsv"
    message:
        "...14. Creating supplementary files..."
    shell:
        "python3.9 14_creating_lipid_network_files.py"


## 15. Creating network object
rule creating_network_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/atoms.csv",
        "network_files/enzyme2reaction.csv",
        "network_files/reaction_urls.tsv",
        "network_files/reaction2align.csv",
        "network_files/gene2reaction.tsv",

        "15_Creating_network_object.R"
    output:
        "network/network.rhea.rda",
        "network/network.rhea.rds"
    message:
        "...15. Creating network object..."
    shell:
        "Rscript 15_Creating_network_object.R"


## 16. Creating lipid network & lipid met.db
rule creating_lipid_network_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/atoms_lipids.tsv",
        "network_files/enzyme2reaction_lipids.tsv",
        "network_files/reaction_urls_lipids.tsv",
        "network_files/reaction2align_lipids.tsv",
        "network_files/mapFromLipidMaps.tsv",
        "network_files/mapFromSwissLipids_ver1.tsv",
        "pre_data/mapFromSpecies.csv",

        "16_Creating_lipid_network.R"
    output:
        "network/network.rhea.lipids.rda",
        "network/network.rhea.lipids.rds",
        "network/met.lipids.db.rda",
        "network/met.lipids.db.rds"
    message:
        "...16. Creating network object..."
    shell:
        "Rscript 16_Creating_lipid_network.R"

