# ==============================
# Helper functions
# ==============================


# C++ function for efficiently downsampling fragments files
Rcpp::sourceCpp(
  code = r"(
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <zlib.h>
#include <cstring>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
void downsample_fragments(const string& input_file, const string& output_file, const vector<string>& valid_barcodes, double proportion, int seed = 42, bool verbose = true) {
    // create a set of valid barcodes for quick lookup
    unordered_set<string> barcode_set(valid_barcodes.begin(), valid_barcodes.end());

    // open the input file
    gzFile gz_input = gzopen(input_file.c_str(), "rb");
    if (!gz_input) {
        Rcerr << "Failed to open input file: " << input_file << endl;
        return;
    }

    // create a map to store the indices of rows for each valid barcode
    unordered_map<string, vector<long>> barcode_indices;

    // read the input file line by line
    char buffer[4096];
    long line_number = 0;
    while (gzgets(gz_input, buffer, sizeof(buffer)) != Z_NULL) {
        char* token;
        char* rest = buffer;
        vector<char*> fields;

        // split the line into fields using tab as the delimiter
        while ((token = strtok_r(rest, "\t", &rest))) {
            fields.push_back(token);
        }

        // check if the line contains exactly 5 columns
        if (fields.size() != 5) {
            Rcerr << "Invalid input file on line " << line_number + 1 << endl;
            gzclose(gz_input);
            return;
        }

        // extract the 4th column (barcode)
        string barcode(fields[3]);

        // store row index if the barcode is valid
        if (barcode_set.find(barcode) != barcode_set.end()) {
            barcode_indices[barcode].push_back(line_number);
        }

        // print progress message for every million lines if verbose is true
        if (verbose && ++line_number % 1000000 == 0) {
            Rcout << "Processed " << line_number << " lines" << endl;
        }
    }

    // check if any valid barcodes were found
    if (barcode_indices.empty()) {
        Rcerr << "None of the barcodes were found in the input file" << endl;
        gzclose(gz_input);
        return;
    }

    // create a set to store the sampled indices
    unordered_set<long> sampled_indices;

    // create a random number generator with the given seed
    mt19937 gen(seed);

    // sample the indices for each barcode
    for (const auto& [barcode, indices] : barcode_indices) {
        size_t sample_size = static_cast<size_t>(indices.size() * proportion);
        vector<long> sampled;
        sample(indices.begin(), indices.end(), back_inserter(sampled), sample_size, gen);
        sampled_indices.insert(sampled.begin(), sampled.end());
    }

    // reset the input file pointer to the beginning
    gzrewind(gz_input);

    // open the output file
    gzFile gz_output = gzopen(output_file.c_str(), "wb");
    if (!gz_output) {
        Rcerr << "Failed to open output file: " << output_file << endl;
        gzclose(gz_input);
        return;
    }

    // read the input file again and write the selected lines to the output file
    line_number = 0;
    while (gzgets(gz_input, buffer, sizeof(buffer)) != Z_NULL) {
        // check if the row index is in the set of sampled indices
        if (sampled_indices.find(line_number) != sampled_indices.end()) {
            gzwrite(gz_output, buffer, strlen(buffer));
        }
        line_number++;
    }

    // close the files
    gzclose(gz_input);
    gzclose(gz_output);
}
)"
)


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
