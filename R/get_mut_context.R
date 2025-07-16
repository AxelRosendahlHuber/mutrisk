# get signatures for individual mutations:
library(GenomicRanges)

#' Get mutation context
#'
#' @param input_muts list with mutations, 5 columns corresponding to: sampleID, chromosome, position,
#'  ref = reference base, alt = alternative allele -> where the mutation is mutated to
#' @param ref_genome Character vector indicating the reference genome of the mutations. Either hg19 or hg38
#'
#' @returns The same mutation list as input, now with an added column containing the trinucleotide (reference base including flanking bases), and triplet (96-trinucelotide format)
#' @export
#'
#' @examples
get_mut_context = function(input_muts, ref_genome = "hg19") {

  chroms = c(1:22, "X", "Y")
  nucs = c("A", "C", "G", "T")
  input_muts = input_muts[,1:5]
  colnames(input_muts) = c("sampleID", "chr", "pos", "ref", "alt")

  input_muts = input_muts |>
    mutate(pos = as.numeric(pos)) |>
    filter(chr %in% chroms) |>
    filter(ref %in% nucs & alt %in% nucs)
  gr = GRanges(seqnames = input_muts$chr, ranges = IRanges(input_muts$pos, width = 1),
               ref = input_muts$ref, alt = input_muts$alt)
  seqlevelsStyle(gr) = "UCSC"

  strand(gr) = ifelse(input_muts$ref %in% c("A", "G"), "-", "+")


  if (ref_genome == "hg19") {
    ref_genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
  } else if (ref_genome == "hg38") {
    ref_genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
  } else (stop("ref_genome variable must either be hg19 or hg38"))

  gr$trinuc = getSeq(ref_genome, gr + 1, as.character = TRUE)

  gr = gr |>
    plyranges::mutate(alt_strand = case_when(ref %in% c("A", "G") & alt == "C" ~ "G",
                                  ref %in% c("A", "G") & alt == "G" ~ "C",
                                  ref %in% c("A", "G") & alt == "T" ~ "A",
                                  ref %in% c("A", "G") & alt == "A" ~ "T",
                           .default = alt),
           ref_strand = case_when(ref == "G" ~ "C",
                                  ref == "A" ~ "T",
                                  .default = ref))

  gr[gr$ref_strand != substr(gr$trinuc,2,2)]

  input_muts$alt_strand = gr$alt_strand
  output_muts = input_muts |>
    mutate(trinuc = gr$trinuc,
           triplet = paste0(substr(trinuc, 1,1), "[", substr(trinuc, 2,2), ">", alt_strand, "]", substr(trinuc, 3,3)))

  return(output_muts)
}
