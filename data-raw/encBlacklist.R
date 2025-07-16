## Code to load the encBlacklist:
# TODO: make accessible for other genomes

# source of the data: encBlacklist - Ameniya et al 2019 Scientific Reports
# https://pubmed.ncbi.nlm.nih.gov/31249361/

# load encBlacklist
library(rtracklayer)
session <- browserSession("UCSC")
genome(session) <- "hg19"
query <- ucscTableQuery(session, track = "Problematic Regions", table = "encBlacklist")
encBlacklist <- getTable(query)
encBlacklist_gr = makeGRangesFromDataFrame(encBlacklist, keep.extra.columns=TRUE)

# get the trinucleotide counts:
usethis::use_data(encBlacklist_gr, overwrite = TRUE)


