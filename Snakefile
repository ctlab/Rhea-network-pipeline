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


## 1. Downloading list of Rhea reactions
rule get_reactions_and_enzymes:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    output:
        "data/rhea2ec.tsv"
    message:
        "...1. Downloading list of Rhea reactions..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2ec.tsv -O {output}"


## 2. Preprocessing RXN files
# 2a. Downloading all RXN files
rule download_rxn_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv"
    output:
        "data/rhea-rxn.tar.gz"
    message:
        "...2. Downloading RXN files..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/ctfiles/rhea%2Drxn.tar.gz -O {output}"


## 2b. Creating LR reactions' list
rule download_rxn_directions:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea-rxn.tar.gz",
        "data/rhea2ec.tsv",

        "2b_getting_LR_rxns_list.py"
    output:
        "data/rhea-directions.tsv",
        "data/LR_Rhea_mapping_rxn.txt",
        "data/undirected_to_LR.tsv"
    message:
        "...2. Creating reactions directions mapping list..."
    shell:
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv -O data/rhea-directions.tsv;"
        "python3.9 2b_getting_LR_rxns_list.py"


rule pre_RDT:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/LR_Rhea_mapping_rxn.txt"
    output:
        "data/reports/rdt_test_done.txt"
    message:
        "...3. Preparing for RDT step..."
    shell:
        "cwd=$(pwd);"
        "mkdir -p data/pre_rdt_output;"
        "cd data/pre_rdt_output;"
        "java -jar -Xmx8G /rdt-2.4.1-jar-with-dependencies.jar -Q SMI -q 'CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O' -g -c -j AAM -f TEXT || true;"
        "touch $cwd/data/reports/rdt_test_done.txt"


## 2c. Moving LR to another folder
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
        "...3. Moving LR reactions..."
    shell:
        'mkdir -p data/LR_rxn;'
        'for file in `cat data/LR_Rhea_mapping_rxn.txt`; do tar -xzf data/rhea-rxn.tar.gz rxn/$file || true; done;'
        "cp ./rxn/*.rxn ./data/LR_rxn/"


## 3. Running RDT
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


## 4. Analyzing RDT
rule performing_RDT_results_analysis:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rdt_done.txt",

        "4_RDT_output_analysis.R"
    output:
        "rdt_analysis/atom_mapping.csv"
    message:
        "...4. Analyzing RDT results and creating atom mapping table..."
    shell:
        "mkdir -p rdt_analysis;"
        "Rscript 4_RDT_output_analysis.R"


## 5. Creating supplementary files for network & metabolites annotation objects for both metabolite and lipid networks
# 5a. Creating supplementary files for metabolite network object
rule creating_network_supplementary_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "rdt_analysis/atom_mapping.csv",
        "data/undirected_to_LR.tsv",

        "5a_creating_network_supplementary_files.py"
    output:
        "data/atom_mapping_C_atoms.csv",
        "data/rpairs_preprocessing.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",
        "network_files/metabolite2atom.csv",
        "network_files/rpairs.csv",
        "network_files/full_atom_mapping_1.csv"
    message:
        "...5. Creating network supplementary files..."
    shell:
        "mkdir -p network_files;"
        "python3.9 5a_creating_network_supplementary_files.py"


# 5b. Preprocessing ChEBI, LipidMaps & SwissLipids databases
rule preprocessing_various_annotations:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv",

        "5b_preprocessing_annotations.py"
    output:
        "data/lipidmaps_to_chebi_to_swisslipids.tsv",
        "data/swisslipids_classes.tsv",
        "data/lipidmaps_classes.tsv"
    message:
        "...5. Creating metabolites annotation supplementary files..."
    shell:
        "mkdir data/databases;"
        "wget 'https://www.lipidmaps.org/files/?file=LMSD&ext=sdf.zip' -O data/databases/LM.zip;"
        "unzip data/databases/LM.zip; mv data/databases/structures.sdf data/databases/structures.txt;"
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf;"
        "mv data/databases/ChEBI_complete_3star.sdf data/databases/ChEBI_complete_3star.txt;"
        "wget 'https://www.swisslipids.org/api/file.php?cas=download_files&file=lipids.tsv' -O data/databases/lipids.tsv;"
        "python3.9 5b_preprocessing_annotations.py"


## 5c. Creating supplementary files for met.db.rhea
rule supp_files_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/rpairs.csv",
        "data/lipidmaps_to_chebi_to_swisslipids.csv",
        "data/swisslipids_classes.tsv",
        "data/lipidmaps_classes",

        "5c_creating_met_db_supplementary_files.py"
    output:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/hmdb.csv",
        "network_files/anomers_suppl.csv",
        "data/kegg2chebi.tsv",
        "data/rhea2kegg_reaction.tsv"
    message:
        "...5. Creating metabolites annotation supplementary files..."
    shell:
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/chemical_data.tsv -O data/chemical_data.tsv;"
        "wget https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/names.tsv.gz -O data/names.tsv.gz;"
        "wget http://rest.kegg.jp/conv/chebi/compound -O data/kegg2chebi.tsv;"
        "wget https://ftp.expasy.org/databases/rhea/tsv/chebi_pH7_3_mapping.tsv -O data/chebi_pH7_3_mapping.tsv;"
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2kegg_reaction.tsv -O data/rhea2kegg_reaction.tsv;"
        "gzip -d data/names.tsv.gz;"
        "python3.9 5c_creating_met_db_supplementary_files.py"


## 5d. Getting mappings for met.db
rule getting_mapping_tables_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/kegg2chebi.tsv",

        "5d_getting_HMDB_KEGG_mappings.R"
    output:
        "data/HMDB2metabolite.csv",
        "data/react2enz_kegg.csv"
    message:
        "...5. Creating metabolites annotation supplementary files..."
    shell:
        "Rscript 5d_getting_HMDB_KEGG_mappings.R"


## 5e. Converting mappings for met.db
rule converting_mappings_for_met_db_rhea:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/kegg2chebi.tsv",
        "data/HMDB2metabolite.csv",

        "5e_creating_HMDB_KEGG_mapping_files.py"
    output:
        "network_files/mapFromHMDB.csv",
        "network_files/mapFromKEGG.csv"
    message:
        "...5. Creating metabolites annotation supplementary files..."
    shell:
        "python3.9 5e_creating_HMDB_KEGG_mapping_files.py"


## 5f. Downloading network supplementary files for genes mapping
rule downloading_supp_gene_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv",
        "data/rhea-directions.tsv"
    output:
        "genes_data/hsa.tsv",
        "genes_data/hsa.proteins.tsv",
        "data/rhea2uniprot_sprot.tsv"
    message:
        "...5. Downloading network supplementary files..."
    shell:
        "mkdir genes_data;"
        "wget 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&force=true&columns=id,entry%20name,reviewed,protein%20names,database(Ensembl),database(GeneID)&compress=yes' -O genes_data/hsa.proteins.tsv.gz;"
        "gzip -d genes_data/hsa.proteins.tsv.gz;"
        "wget https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv -O data/rhea2uniprot_sprot.tsv;"
        "wget http://rest.kegg.jp/link/enzyme/hsa -O genes_data/hsa.tsv"


## 5g. Creating supplementary files for metabolite network object
rule supp_files_for_network:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "data/rhea2ec.tsv",
        "data/rhea-directions.tsv",
        "data/react2enz_kegg.csv",
        "data/rhea2kegg_reaction.tsv",
        "genes_data/hsa.tsv",
        "genes_data/hsa.proteins.tsv",
        "data/rhea2uniprot_sprot.tsv",

        "5g_creating_network_supplementary_files.py"
    output:
        "network_files/enzyme2reaction.csv",
        "network_files/gene2reaction.tsv",
        "network_files/reaction_urls.tsv",
        "data/proteins_genes_unmerged.tsv",
        "data/reactions_annotation.tsv"
    message:
        "...5. Creating network supplementary files..."
    shell:
        "python3.9 5g_creating_network_supplementary_files.py"


## 5h. Creating supplementary files for lipid subnetwork object and corresponding metabolites annotation
rule creating_lipid_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/reaction_urls.tsv",
        "network_files/full_atom_mapping_1.csv",
        "data/proteins_genes_unmerged.tsv",
        "data/reactions_annotation.tsv",

        "5h_creating_lipid_supplementary_files.py"
    output:
        "network_files/atoms_lipids.tsv",
        "network_files/enzyme2reaction_lipids.tsv",
        "network_files/reaction_urls_lipids.tsv",
        "network_files/reaction2align_lipids.tsv",
        "network_files/mapFromLipidMaps.tsv",
        "network_files/mapFromSwissLipids_ver1.tsv"
    message:
        "...5. Creating lipid-specific supplementary files..."
    shell:
        "python3.9 5h_creating_lipid_supplementary_files.py"


## 6a. Create metabolites annotation object
rule creating_met_db_rhea_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/annotation_of_nodes_master_graph.csv",
        "network_files/mapFromHMDB.csv",
        "network_files/mapFromKEGG.csv",
        "network_files/anomers_suppl.csv",

        "6a_creating_met.db.rhea.R"
    output:
        "network/met.rhea.db.rda",
        "network/met.rhea.db.rds"
    message:
        "...6. Creating metabolites annotation object..."
    shell:
        "Rscript 6a_creating_met.db.rhea.R"


## 6b. Creating network object
rule creating_network_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11_v2"
    input:
        "network_files/atoms.csv",
        "network_files/enzyme2reaction.csv",
        "network_files/reaction_urls.tsv",
        "network_files/reaction2align.csv",

        "6b_creating_network_object.R"
    output:
        "network/network.rhea.rda",
        "network/network.rhea.rds"
    message:
        "...6. Creating network object..."
    shell:
        "Rscript 6b_creating_network_object.R"


## 7. Creating lipid network & lipid met.db
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

        "7_creating_lipid_subnetwork_objects.R"
    output:
        "network/network.rhea.lipids.rda",
        "network/network.rhea.lipids.rds",
        "network/met.lipids.db.rda",
        "network/met.lipids.db.rds"
    message:
        "...7. Creating lipid network object and metabolites annotation object..."
    shell:
        "Rscript 7_creating_lipid_subnetwork_objects.R"
