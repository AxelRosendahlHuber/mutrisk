# collection of scripts to calculate mutation rates for whole genome and whole exome

#' Estimate the mutation rates in the whole genome
#'
#' @param trinuc_muts include the trinucleotide
#' @param trinuc_counts Optional: Enter path to alternative trinucleotide counts file for the mappable genome. [default: hg19_trinuc_counts] Using default non-blacklisted hg19 genome.
#' For details, see data-raw/hg38_trinuc_counts.R and data-raw/encBlacklist.R scripts how to obtain the input data
#'
#' @returns
#' @export
#'
#' @examples
estimate_whole_genome_rates = function(trinuc_muts, trinuc_counts = "hg19") {

  if (trinuc_counts == "hg19") {
    # load the whole-genome trinuc counts
    trinuc_counts = hg19_trinuc_counts
  }

  # restrict analyses of whole-genome mutation rates to the mappable genome
  gr_muts = GRanges(trinuc_muts$chr, IRanges(trinuc_muts$pos, width = 1)) |>
    `seqlevelsStyle<-`("UCSC")

  mut_overlaps = findOverlaps(gr_muts, nonbc_genome) # find the overlaps between the mutations and the non-blacklisted genome
  trinuc_muts_ov = trinuc_muts[queryHits(mut_overlaps),]

  trinuc_mut_means = table(trinuc_muts_ov$sampleID, trinuc_muts_ov$triplet) |> colSums()
  trinuc_mut_rates = data.table(triplet = names(trinuc_mut_means), wg_muts = trinuc_mut_means) |>
    mutate(trinucleotide = paste0(substr(triplet, 1,1), substr(triplet, 3,3), substr(triplet, 7,7))) |>
    left_join(trinuc_counts, by = "trinucleotide") |>
    mutate(wg_rate = wg_muts / trinuc_counts)
  trinuc_mut_rates
}

estimate_exome_rates = function(dnds) {

  trinuc_counts = trinuc_counts_dnds(dnds$L)
  mut_counts = trinuc_counts_dnds(dnds$N) |>
    dplyr::rename(mut_counts = trinuc_counts)

  left_join(trinuc_counts, mut_counts) |>
    group_by(triplet, trinuc) |>
    summarize(exome_muts = sum(mut_counts),
              exome_trinucs = sum(trinuc_counts)) |>
    mutate(exome_rate = exome_muts / exome_trinucs) |>
    dplyr::rename(trinucleotide = trinuc)
}


get_wg_exome_rate_estimations = function(trinuc_muts, dnds) {

  if (!("trinuc" %in% colnames(trinuc_muts) & "triplet" %in% colnames(trinuc_muts))) {
    stop("trinuc_muts input needs 'trinuc' and 'triplet' columns.
         These can be obtained using the mutrisk::get_context() function")
  }

  wg_rates = estimate_whole_genome_rates(trinuc_muts = trinuc_muts)
  exome_rates = estimate_exome_rates(dnds)
  left_join(exome_rates, wg_rates)
}

