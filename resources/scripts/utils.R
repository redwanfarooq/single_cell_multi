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


#' Get 10x count matrix
#'
#' Loads 10x count matrix (in HDF5 format).
#'
#' @param file Path to file.
#' @param version Character scalar. 10x HDF5 version ('auto', 'v2', or 'v3').
#' @param cells Character vector. If specified, will subset to matching cell
#' barcodes.
#' @param type Character vector. If specified, will subset to features with matching
#' values in feature type field (v3 only).
#' @param remove.suffix Logical scalar (default `TRUE`). Remove '-1' suffix automatically
#' appended to cell barcodes by Cell Ranger.
#' @param group Optional character scalar. Group name in HDF5 file containing count matrix.
#' If not specified, will use 'matrix' for v3 and the only group for v2.
#'
#' @returns A sparse matrix of counts with features as row names and cell barcodes
#' as column names.
#'
#' @export
get.10x.h5 <- function(file,
                       version = c("auto", "v2", "v3"),
                       cells = NULL,
                       type = NULL,
                       remove.suffix = TRUE,
                       group = NULL) {
  if (!file.exists(file)) stop(file, " does not exist")

  infile <- hdf5r::H5File$new(filename = file, mode = "r")

  version <- match.arg(version)
  if (version == "auto") {
    version <- if (hdf5r::existsGroup(infile, "matrix")) "v3" else "v2"
    message("Detected 10x HDF5 version: ", version)
  }
  features <- if (version == "v2") "genes" else "features/id"

  if (is.null(group)) {
    if (version == "v3") {
      group <- "matrix"
    } else {
      group <- names(infile)
      if (length(group) > 1) stop("Multiple groups detected: ", paste(group, collapse = ", "), ". Please specify 'group' parameter.")
    }
  }
  counts <- infile[[paste(group, "data", sep = "/")]]
  indices <- infile[[paste(group, "indices", sep = "/")]]
  indptr <- infile[[paste(group, "indptr", sep = "/")]]
  shape <- infile[[paste(group, "shape", sep = "/")]]
  features <- infile[[paste(group, features, sep = "/")]]
  barcodes <- infile[[paste(group, "barcodes", sep = "/")]]
  matrix <- Matrix::sparseMatrix(
    i = indices[],
    p = indptr[],
    x = as.numeric(counts[]),
    dims = shape[],
    dimnames = list(features[], barcodes[]),
    index1 = FALSE,
    repr = "C"
  )

  if (remove.suffix) colnames(matrix) <- gsub(pattern = "-1$", replacement = "", colnames(matrix))

  if (!is.null(cells)) matrix <- matrix[, match(cells, colnames(matrix))]
  if (any(is.na(colnames(matrix)))) {
    warning(sum(is.na(colnames(matrix))), " barcode(s) in 'cells' not present in count matrix")
    matrix <- matrix[, !is.na(colnames(matrix))]
  }

  if (version == "v3") {
    feature.type <- infile[["matrix/features/feature_type"]]
    if (!is.null(type)) {
      matrix <- matrix[feature.type[] %in% type, ]
    } else {
      if (length(unique(feature.type[])) > 1) message("Multiple feature types detected: ", paste(unique(feature.type[]), collapse = ", "), ". Returning list of matrices; please specify 'type' parameter to return a single matrix.")
      matrix <- lapply(
        unique(feature.type[]),
        function(type, matrix, feature.type) matrix[grep(pattern = type, x = feature.type), ],
        matrix = matrix,
        feature.type = feature.type[]
      ) |>
        setNames(unique(feature.type[]))
    }

    return(matrix)
  }
}
