library(data.table)


met.rhea.db <- list()


## metabolites
metabolites <- fread("network_files/annotation_of_nodes_master_graph.csv")
metabolites <- metabolites[ , c("metabolite", "metabolite_name", "metabolite_url")]
met.rhea.db$metabolites <- metabolites


## mapFrom
mapFrom <- list()

HMDB2metabolite <- fread("network_files/mapFromHMDB.csv")
HMDB2metabolite <- HMDB2metabolite[, list(HMDB=HMDB, metabolite=ChEBI)]
setkey(HMDB2metabolite, HMDB)

KEGG2metabolite <- fread("network_files/mapFromKEGG.csv")
KEGG2metabolite <- as.data.frame(KEGG2metabolite)
KEGG2metabolite <- unique(KEGG2metabolite)
KEGG2metabolite <- as.data.table(KEGG2metabolite)
KEGG2metabolite <- KEGG2metabolite[, list(KEGG=KEGG, metabolite=ChEBI)]
setkey(KEGG2metabolite, KEGG)

SwissLipids2metabolite <- fread("network_files/mapFromSwissLipids.csv")
SwissLipids2metabolite <- as.data.table(SwissLipids2metabolite)
SwissLipids2metabolite <- SwissLipids2metabolite[, list(SwissLipids=SwissLipids, metabolite=metabolite)]
setkey(SwissLipids2metabolite, SwissLipids)

met.rhea.db$mapFrom <- list("HMDB"=HMDB2metabolite,
                            "KEGG"=KEGG2metabolite,
                            "SwissLipids"=SwissLipids2metabolite)
met.rhea.db$baseId <- "ChEBI"


## baseID
baseId <- "ChEBI"
met.rhea.db$baseId <- baseId


## anomers
anomers <- list()
anomers$metabolite2base_metabolite <- fread("network_files/anomers_suppl.csv")
anomers$base_metabolite2metabolite <- fread("network_files/anomers_suppl.csv")

anomers$metabolite2base_metabolite$metabolite <- as.character(anomers$metabolite2base_metabolite$metabolite)
anomers$metabolite2base_metabolite$base_metabolite <- as.character(anomers$metabolite2base_metabolite$base_metabolite)
anomers$base_metabolite2metabolite$metabolite <- as.character(anomers$base_metabolite2metabolite$metabolite)
anomers$base_metabolite2metabolite$base_metabolite <- as.character(anomers$base_metabolite2metabolite$base_metabolite)

met.rhea.db$anomers <- anomers
setkey(met.rhea.db$anomers$metabolite2base_metabolite, metabolite)
setkey(met.rhea.db$anomers$base_metabolite2metabolite, base_metabolite)

## saving
save(met.rhea.db, file="network/met.rhea.db.rda")
saveRDS(met.rhea.db, file="network/met.rhea.db.rds")
