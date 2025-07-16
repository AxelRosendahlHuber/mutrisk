## Code to prepare `hg38_trinuc_counts` dataset for the hg19 and the non-blacklisted genome
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

# function to get the trinucleotide counts from DNA sequences:
genome <- BSgenome.Hsapiens.UCSC.hg19
genome_ranges <- GRanges(seqnames=seqnames(genome),
                         ranges=IRanges(start=1, end=seqlengths(genome)))

chroms = paste0("chr", c(1:22, "X"))
seqlevels(genome_ranges, pruning.mode = "coarse") = chroms  # trim the genome to only contain the default chromosomes 1 to 22 and X

# get the non-blacklisted genome
nonbc_genome <- setdiff(genome_ranges, reduce(encBlacklist_gr))

genome_sequence = getSeq(Hsapiens, nonbc_genome_ranges)
hg19_trinuc_counts = get_trinuc_counts(genome_sequence)

usethis::use_data(hg19_trinuc_counts, overwrite = TRUE)
usethis::use_data(nonbc_genome, overwrite = TRUE)
