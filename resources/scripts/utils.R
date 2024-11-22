# ==============================
# Helper functions
# ==============================


require(Rcpp)


# C++ function for efficiently downsampling fragments files
Rcpp::sourceCpp("downsample_fragments.cpp")


#' Downsample fragments
#'
#' This function downsamples fragments from a fragments file for a given set of cells.
#'
#' @param input Path to the input fragments file.
#' @param output Path to save the downsampled fragments file.
#' @param cells Character vector of cell barcodes to downsample.
#' @param proportion Proportion of fragments to keep for each cell.
#' @param seed Random seed for reproducibility.
#' @param verbose Print progress messages.
#'
#' @return Path to the downsampled fragments file.
#'
#' @examples
#' downsample.fragments(
#'   input = "fragments.tsv.gz",
#'   output = "downsampled_fragments.tsv.gz",
#'   cells = c("cell1", "cell2", "cell3"),
#'   proportion = 0.5,
#'   seed = 42,
#'   verbose = TRUE
#' )
#'
#' @export
downsample.fragments <- function(input,
                                 output = tempfile(pattern = "fragments", fileext = ".tsv.gz"),
                                 cells,
                                 proportion,
                                 seed = 42L,
                                 verbose = TRUE) {
  if (!file.exists(input)) stop("File not found: ", input)
  if (!is.character(cells)) stop("Argument <cells>: must be a character vector")
  if (!is.numeric(proportion) || proportion < 0 || proportion > 1 || length(proportion) != 1) stop("Argument <proportion>: must be a numeric scalar between 0 and 1")
  if (!is.integer(seed) || length(seed) != 1) stop("Argument <seed>: must be an integer scalar")
  if (verbose) message("Downsampling fragments from ", input, " for ", length(cells), " cells to ", round(proportion * 100), "%")
  # Efficiently downsample fragments file by cell barcode with low memory usage
  downsample_fragments(
    input_file = input,
    output_file = output,
    valid_barcodes = cells,
    proportion = proportion,
    seed = seed,
    verbose = verbose
  )
  message("Output saved to ", output)
  return(output)
}
