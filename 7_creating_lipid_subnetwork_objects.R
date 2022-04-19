library(data.table)

#
## NEWORK.LIPIDS CREATION
#

enzyme2reaction <- fread("network_files/enzyme2reaction_lipids.tsv")
enzyme2reaction[] <- lapply(enzyme2reaction, function(x) {
  if(is.integer(x)) as.character(x) else x})

reactions <- fread("network_files/reaction_urls_lipids.tsv")
reactions[] <- lapply(reactions, function(x) {
  if(is.integer(x)) as.character(x) else x})

colnames(enzyme2reaction) <- c("reaction", "enzyme")
setkey(enzyme2reaction, enzyme)

network <- list()
network$reactions <- reactions
network$enzyme2reaction <- enzyme2reaction

network$reaction2align <- fread("network_files/reaction2align_lipids.tsv")
network$reaction2align[] <- lapply(network$reaction2align, function(x) {
  if(is.integer(x)) as.character(x) else x})

setkey(network$reaction2align, reaction)

network$atoms <- data.table(atom=union(network$reaction2align$atom.x,
                                            network$reaction2align$atom.y))

network$atoms <- fread("network_files/atoms_lipids.tsv")
setkey(network$atoms, atom)

network$metabolite2atom <- network$atoms[, list(metabolite, atom)]
setkey(network$metabolite2atom, metabolite)

save(network, file="network/network.rhea.lipids.rda")
saveRDS(network, file="network/network.rhea.lipids.rds")

#
## MET.LIPIDS.DB CREATION
#

met.lipids.db <- list()

# metabolites
met.lipids.db$metabolites <- fread("network_files/metabolites_lipids.tsv")

# mapFrom
mapFrom <- list()
LipidMaps2metabolite <- fread("network_files/mapFromLipidMaps.tsv")
LipidMaps2metabolite <- as.data.table(LipidMaps2metabolite)
LipidMaps2metabolite <- LipidMaps2metabolite[, list(LipidMaps=LipidMaps, metabolite=metabolite)]
setkey(LipidMaps2metabolite, LipidMaps)

SwissLipids2metabolite <- fread("network_files/mapFromSwissLipids_ver1.tsv")
SwissLipids2metabolite <- as.data.table(SwissLipids2metabolite)
SwissLipids2metabolite <- SwissLipids2metabolite[, list(SwissLipids=SwissLipids, metabolite=metabolite)]
setkey(SwissLipids2metabolite, SwissLipids)

species2metabolite <- fread("pre_data/mapFromSpecies.csv")
species2metabolite <- as.data.table(species2metabolite)
species2metabolite <- species2metabolite[, list(Species=Species, metabolite=metabolite)]
setkey(species2metabolite, Species)

met.lipids.db$mapFrom <- list("LipidMaps"=LipidMaps2metabolite,
                              "SwissLipids"=SwissLipids2metabolite,
                              "Species"=species2metabolite)

# baseID
met.lipids.db$baseId <- "ChEBI"

# anomers
anomers <- list()
anomers$metabolite2base_metabolite <- fread("network_files/anomers_suppl.csv")
anomers$base_metabolite2metabolite <- fread("network_files/anomers_suppl.csv")

anomers$metabolite2base_metabolite$metabolite <-
    as.character(anomers$metabolite2base_metabolite$metabolite)
anomers$metabolite2base_metabolite$base_metabolite <-
    as.character(anomers$metabolite2base_metabolite$base_metabolite)
anomers$base_metabolite2metabolite$metabolite <-
    as.character(anomers$base_metabolite2metabolite$metabolite)
anomers$base_metabolite2metabolite$base_metabolite <-
    as.character(anomers$base_metabolite2metabolite$base_metabolite)

met.lipids.db$anomers <- anomers
setkey(met.lipids.db$anomers$metabolite2base_metabolite, metabolite)
setkey(met.lipids.db$anomers$base_metabolite2metabolite, base_metabolite)

save(met.lipids.db, file="network/met.lipids.db.rda")
saveRDS(met.lipids.db, file="network/met.lipids.db.rds")

print("Pipeline finished")