## code to prepare `signatures` dataset

# Signatures file downloaded from COSMIC:
# https://cancer.sanger.ac.uk/signatures/downloads/ - SBS , 96, v3.4
signatures = read.delim("data-raw/COSMIC_v3.4_SBS_GRCh37.txt") |> tibble::column_to_rownames("Type") |> as.matrix()
signatures = signatures[mutrisk:::TRIPLETS_96,]
usethis::use_data(signatures, overwrite = TRUE)



