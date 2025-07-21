## Code to prepare all the variables for dndscv internal values

# some code taken from: https://github.com/ToolsVanBox/MutationalPatterns/blob/master/R/MutationalPatterns.R
library(data.table)
library(tidyverse)
chroms = paste0("chr", c(1:22, "X"))

SUBSTITUTIONS = mut_types = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
SUBSTITUTIONS_96 <- rep(mut_types, each = 16)

C_TRIPLETS <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

T_TRIPLETS <- c(
  "ATA", "ATC", "ATG", "ATT",
  "CTA", "CTC", "CTG", "CTT",
  "GTA", "GTC", "GTG", "GTT",
  "TTA", "TTC", "TTG", "TTT"
)

CONTEXTS_96 <- c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

# combine substitutions and context in one
TRIPLETS_96 <- paste(substr(CONTEXTS_96, 1, 1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96, 3, 3), sep = "")
TRIPLETS_192 = paste0(rep(TRIPLETS_96, each = 2 ), c('_transcribed', '_untranscribed'))


# map the mutational signatures colors. Also include the "other category in the grey part, to include mutations assigned to other, less frequently occurring signatures
signature_names =  c("SBS1","SBS2","SBS4", "other", "SBS5", "SBS7a", "SBS7b",
                "SBS7d", "SBS10a", "SBS10b", "SBS10c", "SBS10d", "SBS11",
                "SBS13", "SBS15", "SBS16", "SBS18", "SBS19", "SBS25", "SBS31",
                "SBS88", "SBS89", "SBS92", "unassigned")
sig_colors = ggsci::pal_igv()(length(signature_names) + 1)[-4]

names(sig_colors) = signature_names

# Default colours for mutation spectrum plotting
COLORS6 <- c("#2EBAED", "#000000", "#DE1C14",
             "#D4D2D2", "#ADCC54", "#F0D0CE")

# load dndscv substmodel
data(list = sprintf("submod_%s", "192r_3w"), package = "dndscv") # load dndscv substmodel

# match dnds mutation types with the canonical mutation types
mut_types = rownames(substmodel)

library(Biostrings)
# model dndscv substmodel types:
triplet_match_substmodel = data.frame(mut_type = mut_types) |>
  mutate(strand = ifelse(substr(mut_type, 2,2) %in% c("A", "G"), "-", "+"),
         trinuc = substr(mut_type, 1,3),
         mut = substr(mut_type, 6,6)) |>
  mutate(trinuc = case_when(strand == "-" ~ as.character(reverseComplement(DNAStringSet(trinuc))), .default = trinuc),
         mut = case_when(strand == "-" ~ as.character(reverseComplement(DNAStringSet(mut))), .default = mut)) |>
  mutate(triplet = paste0(substr(trinuc, 1,1), "[", substr(trinuc, 2,2), ">", mut, "]", substr(trinuc, 3,3)), type = substr(triplet, 3,5),
         triplet = factor(triplet, levels = TRIPLETS_96))

usethis::use_data(triplet_match_substmodel, overwrite = TRUE)
usethis::use_data(mut_types, overwrite = TRUE)

# use the data internally:
usethis::use_data(triplet_match_substmodel, mut_types, COLORS6, TRIPLETS_96,
                  TRIPLETS_192, chroms, SUBSTITUTIONS, SUBSTITUTIONS_96,
                  cell_muts, metadata,
                  overwrite = TRUE, internal = TRUE)






