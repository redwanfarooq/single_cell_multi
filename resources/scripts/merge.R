#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Merge counts and metadata from multisample multimodal single cell experiments

Usage:
  merge.R --samples=<sample[;sample...]> --output=<output> --hdf5=<path[;path...]> --metadata=<path[;path...]> [--fragments=<path[;path...]>] [options]

Arguments:
  REQUIRED
  --samples=<sample[;sample...]>          Semicolon-separated list of sample IDs
  --output=<output>                       Path to output file (must include file extension 'qs' or 'rds')

  --hdf5=<path[;path...]>                 Path(s) to 10x-formatted HDF5 files containing combined multimodal count matrices
  --metadata=<path[;path...]>             Path(s) to TSV files containing cell metadata
  --fragments=<path[;path...]>            Path(s) to ATAC fragments files (if applicable)

                                          <path[;path...]> arguments must be
                                          EITHER semicolon-separated list of paths (1 per <sample>)
                                          OR template string for path using '{sample}' as placeholder for sample ID

  OPTIONAL
  -t --threads=<int>                      Number of threads [default: 1]
  -l --log=<file>                         Path to log file [default: merge.log]
  --rna-assay=<assay>                     RNA assay name [default: RNA]
  --atac-assay=<assay>                    ATAC assay name [default: ATAC]
  --adt-assay=<assay>                     ADT assay name [default: ADT]
  --gene-types=<type[;type]...>           Semicolon-separated list of Ensembl gene biotypes to keep [default: protein_coding;lincRNA;IG_C_gene;TR_C_gene]
  --control=<name>                        Name of control sample (cross-batch technical replicate)


Options:
  -h --help                               Show this screen
  -b --binarize                           Binarize ATAC counts
  -B --bpcells                            Use BPCells to create compressed on-disk count matrices
  -d --downsample                         Downsample counts to match sequencing depth across libraries
  -e --exclude-control                    Exclude control sample (cross-batch technical replicate)
  -q --quiet                              Do not print logging messages to console
"

# Parse options
opt <- docopt::docopt(DOC)

# Logging options
logger::log_appender(logger::appender_file(opt[["log"]]), index = 1)
logger::log_layout(logger::layout_glue, index = 1)
if (!opt[["quiet"]]) {
  logger::log_appender(logger::appender_console, index = 2)
  logger::log_layout(logger::layout_glue_colors, index = 2)
}
logger::log_warnings()
logger::log_errors()

# Extract parameters from command line options
logger::log_info("Parsing command line arguments")
params <- lapply(opt, function(x) if (is.character(x)) stringr::str_split_1(string = x, pattern = ";") else x)

# Check parameters
for (x in c("hdf5", "metadata", "fragments")) {
  if (!is.null(params[[x]]) && length(params[[x]]) != length(params$samples)) {
    if (length(params[[x]]) != 1) stop("Argument <", x, ">: number of paths must be equal to number of samples") # mismatch between number of paths and samples, if exactly 1 path provided assume it is a template string
    if (!grepl(pattern = "\\{sample\\}", x = params[[x]])) stop("Argument <", x, ">: invalid path template string '", params[[x]], "'") # template string does not contain placeholder
  }
  params[[x]] <- glue::glue(params[[x]], sample = params$samples) # replace placeholder(s) if path provided is template string
}
if (!is.null(params$fragments)) {
  missing.files <- vapply(params$fragments, function(x) !file.exists(x), logical(1))
  for (file in params$fragments[missing.files]) {
    warning("File not found: ", file)
  }
  params$fragments <- params$fragments[!missing.files]
}
if (params$exclude_control && is.null(params$control)) stop("Argument <control>: must be provided if running with '--exclude-control' option")
if (!grepl(pattern = "\\.qs$|\\.rds$", x = params$output, ignore.case = TRUE)) stop("Argument <output>: invalid file extension; must be either 'qs' or 'rds'")

# Log options
for (x in grep(pattern = "^--|^help", names(params), value = TRUE, invert = TRUE)) {
  y <- if (!is.null(names(params[[x]])) && length(names(params[[x]]))) {
    mapply(
      function(k, v) paste(k, v, sep = "="),
      names(params[[x]]),
      params[[x]],
      SIMPLIFY = TRUE
    )
  } else {
    params[[x]]
  }
  logger::log_info("{x}: {paste(y, collapse = ';')}")
}


# ==============================
# SETUP
# ==============================
logger::log_info("Initialising")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(furrr)
  library(DropletUtils)
  library(Seurat)
  library(Signac)
  library(BPCells)
})

options(future.globals.maxSize = 1000000 * 1024^2, UCSC.goldenPath.url = "http://hgdownload.soe.ucsc.edu/goldenPath")
furrr.options <- furrr_options(seed = 42, scheduling = FALSE)

if (!dir.exists(dirname(params$output))) dir.create(dirname(params$output), recursive = TRUE)


# ==============================
# SCRIPT
# ==============================
plan(multicore, workers = as.integer(params$threads))

logger::log_info("Loading and processing data")
# Load cell metadata and subset to filtered barcodes
cell.metadata <- future_map(
  .x = params$metadata,
  .f = function(file) {
    if (!params$quiet) message("Loading metadata: ", file)
    df <- vroom::vroom(
      file,
      delim = "\t",
      col_names = TRUE,
      show_col_types = FALSE
    ) %>%
      tibble::column_to_rownames("barcode") %>%
      as.data.frame()
    return(df)
  },
  .options = furrr.options
) %>%
  setNames(params$samples)


# Load count matrices and peform preprocessing steps
mat <- future_map(
  .x = params$hdf5,
  .f = function(file) {
    if (!params$quiet) message("Loading count matrices: ", file)
    x <- Read10X_h5(file, use.names = FALSE, unique.features = FALSE) %>% suppressMessages()
    return(x)
  },
  .options = furrr.options
) %>%
  setNames(params$samples) %>%
  transpose() %>%
  setNames(case_match(names(.), "Gene Expression" ~ "RNA", "Peaks" ~ "ATAC", "Antibody Capture" ~ "ADT"))
# Clean up feature names in RNA assays
# convert RNA feature names from Ensembl ID to gene symbol and extract feature metadata
# only keep genes annotated as protein-coding or lincRNA
# adapted from source code for Azimuth:::ConvertEnsemblToSymbol but modified to use EnsDb.Hsapiens.v86
# NOTE: this method leads to loss of features if Ensembl IDs do not map to a gene symbol in database
#       an alternative approach would be to use the feature metadata in features.tsv.gz (from CellRanger/STARsolo output)
logger::log_info("Cleaning up feature names")
res <- future_map(
  .x = mat$RNA,
  .f = function(x, gene.types = params$gene_types) {
    df <- data.frame(rownames = rownames(x))
    df$ensembl_id <- sub(pattern = "[.][0-9]*", replacement = "", x = df$rownames)
    mapping <- ensembldb::select(
      EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
      keys = df$ensembl_id,
      keytype = "GENEID",
      columns = c("GENEID", "SYMBOL", "GENEBIOTYPE")
    ) %>%
      rename(
        GENEID = "ensembl_id",
        SYMBOL = "gene_symbol",
        GENEBIOTYPE = "gene_type"
      )
    df <- left_join(df, mapping, by = "ensembl_id") %>% tibble::column_to_rownames("rownames")
    df <- df[rownames(x), ] %>%
      filter(!is.na(gene_symbol), grepl(pattern = paste(gene.types, collapse = "|"), x = gene_type)) %>%
      mutate(gene_symbol = make.unique(gsub(pattern = "_", replacement = "", x = gene_symbol))) %>%
      arrange(gene_symbol)
    x <- x[rownames(df), ]
    rownames(x) <- df$gene_symbol
    rownames(df) <- df$gene_symbol
    return(list(counts = x, metadata = df))
  },
  .options = furrr.options
)
mat$RNA <- future_map(
  .x = res,
  .f = function(x) x$counts,
  .options = furrr.options
)
gene.metadata <- future_map(
  .x = res,
  .f = function(x) x$metadata,
  .options = furrr.options
)
# Remove isotype control antibody tags
if (!is.null(mat$ADT)) {
  logger::log_info("Removing isotype control antibody tags")
  mat$ADT <- future_map(
    .x = mat$ADT,
    .f = function(x) {
      x <- x[!grepl(pattern = "isotype|control", x = rownames(x), ignore.case = TRUE), ]
      return(x)
    },
    .options = furrr.options
  )
}
# Additional preprocessing steps for ATAC assays
if (!is.null(mat$ATAC)) {
  # Combine peaks across samples
  logger::log_info("Finding and quantifying merged ATAC peak set")
  # extract genomic ranges of peaks from rownames and merge overlapping peaks across samples
  peaks <- future_map(
    .x = mat$ATAC,
    .f = function(x) {
      stringr::str_split(string = rownames(x), pattern = "[:-]", simplify = TRUE) %>%
        as.data.frame() %>%
        setNames(c("chr", "start", "stop")) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") # remove peaks on non-standard chromosomes
    },
    .options = furrr.options
  ) %>%
    Reduce(f = c, x = .) %>%
    GenomicRanges::reduce()
  # filter out bad peaks based on width
  peaks <- peaks[width(peaks) < 10000 & width(peaks) > 20]
  # create fragments objects
  fragments <- future_pmap(
    .l = list(
      x = mat$ATAC,
      sample = names(mat$ATAC),
      path = params$fragments
    ),
    .f = function(x, sample, path) {
      if (!params$quiet) message("Getting fragments for ", sample)
      CreateFragmentObject(
        path = path,
        cells = colnames(x),
        verbose = !params$quiet
      )
    },
    .options = furrr.options
  ) %>%
    setNames(names(mat$ATAC))
  # count fragments in merged peaks
  mat$ATAC <- pmap(
    .l = list(
      x = mat$ATAC,
      sample = names(mat$ATAC),
      fragment.object = fragments
    ),
    .f = function(x, sample, fragment.object, merged.peaks = peaks) {
      if (!params$quiet) message("Counting fragments in merged peaks for ", sample)
      FeatureMatrix(
        fragments = fragment.object,
        features = merged.peaks,
        cells = colnames(x),
        sep = c(":", "-"),
        verbose = !params$quiet
      )
    }
  )
  if (params$binarize) {
    # Binarize counts
    logger::log_info("Binarizing ATAC counts")
    mat$ATAC <- future_map(
      .x = mat$ATAC,
      .f = BinarizeCounts,
      .options = furrr.options
    )
  }
}
if (params$downsample) {
  # Downsample counts to match sequencing depth across libraries
  mat <- map2(
    .x = mat,
    .y = names(mat),
    .f = function(x, type) {
      if (type %in% c("RNA", "ADT") || !params$binarize) { # do not downsample ATAC counts if binarized
        set.seed(42)
        logger::log_info("Downsampling {type} counts to match sequencing depth across {length(x)} libraries")
        x <- as.list(downsampleBatches(x))
      }
      return(x)
    }
  )
}
mat <- transpose(mat)

if (params$bpcells) {
  # Convert to BPCells matrices
  logger::log_info("Creating compressed on-disk count matrices")
  mat <- future_map2(
    .x = mat,
    .y = names(mat),
    .f = function(x, sample) {
      map2(
        .x = x,
        .y = names(x),
        .f = function(m, type) {
          if (type == "RNA") { # only convert RNA matrices (IterableMatrix not currently supported for ChromatinAssay or CLR normalization)
            if (!params$quiet) message("Writing matrix: ", glue::glue(file.path(dirname(params$output), "{sample}_{type}")))
            bpcells <- write_matrix_dir(
              mat = convert_matrix_type(m, type = "uint32_t"), # convert from dgCMatrix to IterableMatrix
              dir = glue::glue(file.path(dirname(params$output), "{sample}_{type}")),
              overwrite = TRUE
            )
            return(bpcells)
          }
          return(m)
        }
      )
    },
    .options = furrr.options
  )
}


# Create merged Seurat object
logger::log_info("Creating merged Seurat object")
seu <- future_pmap(
  .l = list(
    x = mat,
    sample = names(mat),
    cell.metadata = cell.metadata,
    gene.metadata = gene.metadata
  ),
  function(x, sample, cell.metadata, gene.metadata, fragment.objects = eval(if (exists("fragments")) fragments else NULL)) {
    obj <- CreateSeuratObject(counts = x$RNA, assay = "RNA", project = sample, meta.data = cell.metadata)
    obj[["RNA"]] <- AddMetaData(object = obj[["RNA"]], metadata = gene.metadata)
    if (!is.null(x$ATAC)) {
      if (is.null(fragment.objects[[sample]])) stop("Missing fragments file for ", sample)
      obj[["ATAC"]] <- CreateChromatinAssay(
        counts = x$ATAC,
        sep = c(":", "-"),
        fragments = fragment.objects[[sample]]
      )
    }
    if (!is.null(x$ADT)) obj[["ADT"]] <- CreateAssay5Object(counts = x$ADT)
    return(obj)
  },
  .options = furrr.options
)
seu <- if (length(seu) > 1) merge(x = seu[[1]], y = seu[seq.int(2, length(seu), by = 1)], add.cell.ids = params$samples) else seu[[1]]

# Join layers
for (x in Assays(seu)) {
  if (is(seu[[x]], "Assay5")) { # joining and splitting layers only works for Assay v5
    seu[[x]] <- JoinLayers(object = seu[[x]], layers = "counts")
  }
}

# Remove cells from control sample (if specified)
if (params$exclude_control) {
  logger::log_info("Removing cells from control sample: {params$control}")
  if (!any(seu$hash_id == params$control)) warning("No matching samples found")
  seu <- subset(seu, hash_id != params$control)
}


# Save output
logger::log_info("Saving output")
if (grepl(pattern = "\\.qs$", x = params$output, ignore.case = TRUE)) {
  qs::qsave(
    seu,
    file = params$output,
    nthreads = as.integer(params$threads)
  ) %>%
    suppressWarnings()
} else {
  saveRDS(
    seu,
    file = params$output
  ) %>%
    suppressWarnings()
}
logger::log_success("Output file: {params$output}")
