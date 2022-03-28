library(data.table)
library(metaboliteIDmapping)
library(KEGGREST)

HMDB2metabolite <- metabolitesMapping[(!is.na(metabolitesMapping$ChEBI)),
                                      c("ChEBI", "HMDB")]
HMDB2metabolite <- na.omit(HMDB2metabolite)
HMDB2metabolite <- as.data.frame(HMDB2metabolite)
HMDB2metabolite <- unique(HMDB2metabolite)
HMDB2metabolite <- as.data.table(HMDB2metabolite)
HMDB2metabolite$ChEBI <- paste0("CHEBI:", HMDB2metabolite$ChEBI)

write.csv(HMDB2metabolite, "data/HMDB2metabolite.csv", row.names = F)

reaction2enzyme <- keggLink("enzyme", "reaction")
reaction2enzyme <- data.table(
    reaction=gsub("rn:", "", names(reaction2enzyme), fixed = T),
    enzyme=gsub("ec:", "", reaction2enzyme, fixed = T))

write.csv(reaction2enzyme, file="data/react2enz_kegg.csv", row.names = FALSE)