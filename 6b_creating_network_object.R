library(data.table)


enzyme2reaction <- fread("network_files/enzyme2reaction.csv")
enzyme2reaction[] <- lapply(enzyme2reaction, function(x) {
  if(is.integer(x)) as.character(x) else x})

reactions <- fread("network_files/reaction_urls.tsv")
reactions[] <- lapply(reactions, function(x) {
  if(is.integer(x)) as.character(x) else x})

colnames(enzyme2reaction) <- c("reaction", "enzyme")
setkey(enzyme2reaction, enzyme)

network <- list()
network$reactions <- reactions
network$enzyme2reaction <- enzyme2reaction

network$reaction2align <- fread("network_files/reaction2align.csv")
network$reaction2align[] <- lapply(network$reaction2align, function(x) {
  if(is.integer(x)) as.character(x) else x})

setkey(network$reaction2align, reaction)

network$atoms <- data.table(atom=union(network$reaction2align$atom.x,
                                            network$reaction2align$atom.y))

network$atoms <- fread("network_files/atoms.csv")
setkey(network$atoms, atom)

network$metabolite2atom <- network$atoms[, list(metabolite, atom)]
setkey(network$metabolite2atom, metabolite)

save(network, file="network/network.rhea.rda")
saveRDS(network, file="network/network.rhea.rds")