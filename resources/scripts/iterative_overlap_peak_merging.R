# Functions for iterative overlap peak merging algorithm from Corces, M.R. et al. (2018). The chromatin accessibility landscape of primary human cancers. Science, 362, eaav1898.
# Peak score normalisation method modified from original "score per million" to quantile-rank normalisation
# Original code writted by Jeff Granja and Ryan Corces (https://github.com/corceslab/ATAC_IterativeOverlapPeakMerging)
# Modified by Redwan Farooq


readSummits <- function(file) {
  df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr", "start", "end", "name", "score"))))
  df <- df[, c(1, 2, 3, 5)] # remove name column
  return(GenomicRanges::makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE))
}


getQuantiles <- function(x) {
  return(floor(rank(x)) / length(x))
}


clusterGRanges <- function(gr,
                           filter = TRUE,
                           by = "score",
                           decreasing = TRUE,
                           verbose = TRUE) {
  # find overlapping sets of peaks
  gr <- sort(GenomeInfoDb::sortSeqlevels(gr))
  r <- GenomicRanges::reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
  o <- IRanges::findOverlaps(gr, r)
  S4Vectors::mcols(gr)$cluster <- S4Vectors::subjectHits(o)
  if (verbose) message(sprintf("Found %s overlaps...", length(gr) - max(S4Vectors::subjectHits(o))))
  # filter overlapping peaks by metadata column (if available) or by order (first peak in overlapping set is retained)
  if (filter) {
    if (length(gr) - max(S4Vectors::subjectHits(o)) > 0) {
      if (by %in% colnames(S4Vectors::mcols(gr))) {
        if (verbose) message(sprintf("Filtering overlaps by %s...", by))
        gr <- gr[order(S4Vectors::mcols(gr)[, by], decreasing = decreasing), ]
        grn <- gr[!duplicated(S4Vectors::mcols(gr)$cluster), ]
        gr <- sort(GenomeInfoDb::sortSeqlevels(grn))
      } else {
        if (verbose) message(sprintf("Filtering overlaps by order..."))
        gr <- gr[!duplicated(S4Vectors::mcols(gr)$cluster), ]
      }
    }
    S4Vectors::mcols(gr)$cluster <- NULL
  }
  return(gr)
}


nonOverlappingGRanges <- function(gr,
                                  by = "score",
                                  decreasing = TRUE,
                                  verbose = TRUE) {
  if (!by %in% colnames(S4Vectors::mcols(gr))) warning("'", by, "' is not a metadata column in 'gr'; will merge overlapping peaks by order.")
  i <- 0
  gr_initial <- gr
  if (verbose) message("Merging overlapping peaks")
  while (length(gr_initial) > 0) {
    i <- i + 1
    if (verbose) message(sprintf("Iteration %d", i))
    gr_clustered <- clusterGRanges(gr = gr_initial, filter = TRUE, by = by, decreasing = decreasing, verbose = verbose)
    gr_initial <- IRanges::subsetByOverlaps(gr_initial, gr_clustered, invert = TRUE)
    if (i == 1) {
      gr_all <- gr_clustered
    } else {
      gr_all <- c(gr_all, gr_clustered)
    }
  }
  if (verbose) message("Done")
  return(gr_all)
}
