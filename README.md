# Pipeline for the creation of Rhea metabolic network & Rhea lipidomic subnetwork

This pipeline produces network objects and corresponding metabolites annotation for metabolite 
and lipid networks based on [Rhea database](https://www.rhea-db.org) which are used for analysis by 
[Shiny GATOM](https://artyomovlab.wustl.edu/shiny/gatom). The pipeline uses Snakemake and Singularity 
for execution and uses the provided scripts while running.


## Requirements

* Snakemake 
* Singularity v3.5.3
* At least 10 GB of disk space 


## Quick Start

* Clone the repository:

```git clone git@github.com:ctlab/Rhea-network-pipeline.git```

* Activate the environment with snakemake and singularity:

```conda activate snakemake```

* Execute the pipeline:

```snakemake --use-singularity --cores 8```

Execution will need at least 10 GB of empty space and can take from 2 to 10 
days to run depending on the setup.


## Pipeline structure

The pipeline is designed to work in specifically constructed docker image which can
be found [here](https://hub.docker.com/repository/docker/mariaembio/python3.9_r4.0_java11_v2).
It contains all necessary dependencies for Python and R code as well as Reaction Decoder Tool,
and Snakemake will use it automatically when executing.

Some steps of the pipeline will execute included scripts that are written in Python or R, 
while simpler steps will only run bash commands from Snakefile. 

Pipeline takes as input files that are stored in `pre_data` folder:
* `polymers_to_chebi.tsv` file contains mappings between ChEBI IDs and Rhea polymer IDs;
* `mapFromSpecies.csv` file contains preprocessed mapping tables between lipid species and 
ChEBI IDs.

Network `rds` object and metabolites annotation `rds` object for two kinds of networks 
-- metabolite network and lipid subnetwork -- are considered to be the output of the pipeline. 
The files will be stored in `network` folder.


## Snakemake pipeline steps

1. List of undirected Rhea reactions is downloaded;
2. In order to construct atom network, we need to perform atom mapping with Reaction Decoder Tool.
The Reaction Decoder Tool takes as input RXN files, however, undirected reactions do not have RXN files. 
Thus, the following steps are done:

    a. All RXN-files for Rhea reactions are downloaded;

    b. Mapping table between undirected, left-to-right and right-to-left reactions IDs is downloaded to 
distinguish left-to-right and right-to-left RXN files;

    c. Only left-to-right reactions are kept;
3. The [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) is used for atom mapping of 
the RXN files;
4. Atom mapping tables are created with the use of 
[ChemmineR](https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html) 
from RXN files processed with the Reaction Decoder Tool;
5. The supplementary annotation files for the metabolite and lipid network and corresponding 
metabolites are created;
    - Metabolites annotation includes metabolites ChEBI ID & link extracted from ChEBI, 
   HMDB to ChEBI mapping obtained via metaboliteIDmapping R-package and ChEBI to KEGG mapping 
   extracted from KEGG;
    - Network annotation includes reactions Rhea ID & link extracted from Rhea, 
   reaction-enzyme mapping obtained from Rhea & processed atom mapping data;
6. The network & metabolites annotation objects for the metabolite and lipid network are created:

    a. Metabolites annotation object for metabolite network is created;

    b. Network-object for metabolite network is created;
7. The network and corresponding metabolites annotation objects for the lipid subnetwork are created.
