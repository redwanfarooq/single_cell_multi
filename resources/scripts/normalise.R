#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Normalise counts and select variable features from multisample multimodal single cell experiments

Usage:
  normalise.R --input=<input> --output=<output> [options]

Arguments:
  REQUIRED
  --input=<input>                         Path to RDS or QS file containing merged Seurat object with multimodal count matrices and metadata
  --output=<output>                       Path to output file (must include file extension 'qs' or 'rds')

  OPTIONAL
  -t --threads=<int>                      Number of threads [default: 1]
  -l --log=<file>                         Path to log file [default: normalise.log]
  --exp=<name>                            Name of metadata field specifying the original experiment [default: orig.ident]
  --rna-assay=<assay>                     RNA assay name [default: RNA]
  --atac-assay=<assay>                    ATAC assay name [default: ATAC]
  --adt-assay=<assay>                     ADT assay name [default: ADT]
  --n-features=<int>                      Maximum number of variable features to select [default: 3000]
  --rna-filter-features=<pattern>         Regular expression pattern to filter RNA variable features (matching features will be removed) [default: ^MT-|^RP[SL]]
  --atac-filter-features=<pattern>        Regular expression pattern to filter ATAC variable features (matching features will be removed) [default: ^chrM|^chrY]
  --clr=<method>                          Centred-log ratio (CLR) transformation method for ADT data ('seurat', 'original') [default: seurat]


Options:
  -h --help                               Show this screen
  -M --regress-percent-mitochondrial      Regress out percentage mitochondrial gene expression
  -C --regress-cell-cycle                 Regress out cell cycle scores
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
  library(Seurat)
  library(Signac)
})

options(future.globals.maxSize = 1000000 * 1024^2)
furrr.options <- furrr_options(seed = 42, scheduling = FALSE)

if (!dir.exists(dirname(params$output))) dir.create(dirname(params$output), recursive = TRUE)


# ==============================
# SCRIPT
# ==============================
plan(multicore, workers = as.integer(params$threads))

logger::log_info("Loading data")
# Load data
if (grepl(pattern = "\\.qs$", x = params$input, ignore.case = TRUE)) {
  seu <- qs::qread(file = params$input, nthreads = as.integer(params$threads))
} else {
  seu <- readRDS(file = params$input)
}


# Get experiment variable
params$exp <- params$exp %>% match.arg(choices = colnames(seu[[]]))


# Perform normalisation and variable feature selection
# Normalise RNA data
if (params$rna_assay %in% Assays(seu)) {
  logger::log_info("Normalising RNA data")
  assay <- params$rna_assay
  vars.to.regress <- NULL
  if (params$regress_percent_mitochondrial) vars.to.regress <- c(vars.to.regress, "pctMito_RNA")
  if (params$regress_cell_cycle) vars.to.regress <- c(vars.to.regress, "S.Score", "G2M.Score")
  suppressWarnings({
    seu <- seu %>%
      NormalizeData(assay = assay, normalization.method = "LogNormalize", verbose = !params$quiet) %>%
      PercentageFeatureSet(pattern = "^MT-", col.name = "pctMito_RNA") %>%
      CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, search = TRUE, verbose = !params$quiet) %>%
      split(f = seu[[params$exp, drop = TRUE]], assay = assay, layers = c("counts", "data")) %>%
      FindVariableFeatures(assay = assay, nfeatures = as.integer(params$n_features), selection.method = "vst", verbose = !params$quiet) %>%
      ScaleData(assay = assay, vars.to.regress = vars.to.regress, verbose = !params$quiet)
  }) # suppress unhelpful warnings
  min.depth <- map_dbl(
    .x = Layers(object = seu, assay = assay, search = "counts"),
    .f = function(x, obj = seu) median(colSums(LayerData(object = obj, assay = assay, layer = x)))
  ) %>%
    min() %>%
    floor()
  seu <- seu %>%
    SCTransform(
      assay = assay,
      variable.features.n = as.integer(params$n_features),
      vars.to.regress = vars.to.regress,
      vst.flavor = "v2",
      min_cells = 0, # use all genes for normalization
      scale_factor = min.depth, # use minimum sequencing depth to ensure corrected counts are comparable across batches
      verbose = !params$quiet
    ) %>%
    SetAssayData(
      assay = "SCT",
      layer = "counts",
      new.data = MatrixExtra::filterSparse(LayerData(object = ., assay = "SCT", layer = "counts"), fn = function(x) !is.nan(x))
    ) %>%
    SetAssayData(
      assay = "SCT",
      layer = "data",
      new.data = log1p(LayerData(object = ., assay = "SCT", layer = "counts"))
    ) # set NaN values in corrected counts matrix to 0 (from genes with zero counts in all cells in a batch)
  seu[["SCT"]] <- AddMetaData(seu[["SCT"]], seu[[assay]][[c("ensembl_id", "gene_symbol", "gene_type")]])
  for (x in c(assay, "SCT")) {
    VariableFeatures(seu, assay = x) <- grep(pattern = params$rna_filter_features, x = VariableFeatures(seu, assay = x), value = TRUE, invert = TRUE) # remove unwanted features
  }
}

# Normalise ATAC data
if (params$atac_assay %in% Assays(seu)) {
  logger::log_info("Normalising ATAC data")
  assay <- params$atac_assay
  # Temporarily switch default assay to 'ATAC' (DietSeurat does not allow removal of default assay)
  default <- DefaultAssay(seu)
  DefaultAssay(seu) <- assay
  # Split object by experiment prior to normalisation/feature selection (ChromatinAssay does not support layers)
  object.list <- seu %>%
    DietSeurat(assays = assay, layers = "counts") %>%
    SplitObject(split.by = params$exp) %>%
    future_map2(
      .x = .,
      .y = names(.),
      .f = function(obj, name) {
        if (!params$quiet) message("Normalizing layer: counts.", name)
        obj <- RunTFIDF(object = obj, assay = assay, verbose = !params$quiet) %>% suppressWarnings() # suppress unhelpful warnings
        return(obj)
      },
      .options = furrr.options
    )
  LayerData(object = seu, assay = assay, layer = "data") <- LayerData(
    object = if (length(object.list) > 1) merge(x = object.list[[1]], y = object.list[seq.int(2, length(object.list), by = 1)]) else object.list[[1]],
    assay = assay,
    layer = "data"
  )
  VariableFeatures(seu, assay = assay) <- grep(pattern = params$atac_filter_features, x = Features(seu[[assay]]), value = TRUE, invert = TRUE) # use all features (except unwanted features) as variable features
  DefaultAssay(seu) <- default
}

# Normalise ADT data
if (params$adt_assay %in% Assays(seu)) {
  logger::log_info("Normalising ADT data")
  params$clr <- params$clr %>% match.arg(choices = c("seurat", "original"))
  assay <- params$adt_assay
  seu <- split(seu, f = seu[[params$exp, drop = TRUE]], assay = assay, layers = c("counts", "data"))
  VariableFeatures(object = seu, assay = assay) <- map(
    .x = Layers(object = seu, assay = assay, search = "counts"),
    .f = function(x, obj = seu) {
      counts <- LayerData(object = obj, assay = assay, layer = x)
      return(rownames(counts[rowSums(counts) > 0, ]))
    }
  ) %>%
    Reduce(f = intersect, x = ., init = Features(seu[[assay]])) # use all common features as variable features (to account for different ADT panels across experiments)
  if (params$clr == "seurat") {
    seu <- seu %>%
      NormalizeData(assay = assay, normalization.method = "CLR", margin = 2, verbose = !params$quiet) %>%
      ScaleData(assay = assay, verbose = !params$quiet)
  } else if (params$clr == "original") {
    for (x in Layers(object = seu, assay = assay, search = "counts")) {
      if (!params$quiet) message("Normalizing layer: ", x)
      LayerData(
        object = seu,
        assay = assay,
        layer = sub(pattern = "counts", replacement = "data", x = x)
      ) <- apply(
        X = LayerData(object = seu, assay = assay, layer = x),
        MARGIN = 2,
        FUN = function(x) log(x + 1 / exp(sum(log(x + 1)) / length(x))) # original centred-log ratio (CLR) transformation (without Seurat modifications)
      )
    }
    seu <- ScaleData(object = seu, assay = assay, verbose = !params$quiet)
  }
}


# Clean up Seurat object
logger::log_info("Cleaning up merged Seurat object")
for (assay in Assays(seu)) {
  if (is(seu[[assay]], "Assay5")) {
    seu[[assay]] <- JoinLayers(object = seu[[assay]], layers = c("counts", "data"))
  }
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
