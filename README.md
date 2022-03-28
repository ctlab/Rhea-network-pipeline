## Pipeline for the creation of Rhea metabolic network & Rhea lipidomic subnetwork

This pipeline produces network-objects and corresponding metabolite annotations based on
Rhea database (https://www.rhea-db.org) which are used for analysis by Shiny GATOM
(https://artyomovlab.wustl.edu/shiny/gatom).

## Docker

The pipeline is designed to work in specifically constructed docker container which can
be found here: https://hub.docker.com/repository/docker/mariaembio/python3.9_r4.0_java11_v2.

## Execution of the pipeline 

The code is written in Python and R, and designed to be executed with Snakemake and 
singularity. Execution will need at least 3 GB of empty space and can take from 2 to 10 
days to run.

## Data

Pipeline takes as input files that are stored in `pre_data` folder. These files include 
mappings between ChEBI IDs and Rhea polymer IDs as well as preprocessed SwissLipids database 
data (https://www.swisslipids.org).

## Pipeline description

1. List of Rhea undirected reactions, also called master reactions, is downloaded;
2. All undirected reactions IDs are translated into LR-reactions and then corresponding 
RXN files are downloaded;
3. Then Reaction Decoder Tool (https://github.com/asad/ReactionDecoder) is used for atom 
mapping of the RXN files;
4. Atom mapping tables are created with the use of ChemmineR
(https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html);
5. Various supplementary files for metabolic and lipidomic networks are created, including
Rhea reactions' annotation and ChEBI metabolites' annotation;
6. Network-objects and corresponding metabolites' annotation object for Rhea metabolic 
network and lipidomic subnetworks are created.
