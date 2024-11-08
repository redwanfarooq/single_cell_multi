#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Integrate data from multisample multimodal single cell experiments

Usage:
  integrate.R --input=<input> --output=<output> [options]

Arguments:
  REQUIRED
  --input=<input>                         Path to RDS or QS file containing merged Seurat object with multimodal count matrices and metadata
  --output=<output>                       Path to output file (must include file extension 'qs' or 'rds')

  OPTIONAL
  -t --threads=<int>                      Number of threads [default: 1]
  -l --log=<file>                         Path to log file [default: integration.log]
  --batch=<name[;name...]>                Name(s) of metadata field(s) to use for batch variable [default: orig.ident]
  --rna-assay=<assay>                     RNA assay name [default: RNA]
  --atac-assay=<assay>                    ATAC assay name [default: ATAC]
  --adt-assay=<assay>                     ADT assay name [default: ADT]
  --normalisation-method=<method>         Normalisation method used for RNA assay ('LogNormalize' or 'SCT') [default: LogNormalize]
  --integration-method=<method>           Integration method function name ('RPCAIntegration', 'JointPCAIntegration', 'CCAIntegration', 'HarmonyIntegration', 'FastMNNIntegration' or 'scVIIntegration')
  --integration-args=<str>                Additional named arguments to pass to integration method function


Options:
  -h --help                               Show this screen
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
  library(future)
  library(Seurat)
  library(SeuratWrappers)
  library(Signac)
})

source("integration_functions.R")

options(future.globals.maxSize = 1000000 * 1024^2)

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


# Get/create batch variable
params$batch <- params$batch %>% match.arg(choices = colnames(seu[[]]), several.ok = TRUE)
if (length(params$batch) > 1) {
  seu$batch <- map(.x = params$batch, .f = function(x) seu[[x]][[1]]) %>% Reduce(f = function(x, y) paste(x, y, sep = "-"), x = .)
} else {
  seu$batch <- seu[[params$batch]][[1]]
}
# Merge and re-split counts layer by batch variable
for (x in Assays(seu)) {
  if (is(seu[[x]], "Assay5")) { # splitting layers only works for Assay v5
    seu[[x]] <- JoinLayers(object = seu[[x]], layers = "counts")
    seu[[x]] <- split(
      seu[[x]],
      f = seu$batch,
      layers = c("counts", "data")
    )
  }
}


# Perform PCA and integration steps
# RNA data
logger::log_info("Performing PCA for RNA data")
params$normalisation_method <- params$normalisation_method %>% match.arg(choices = c("LogNormalize", "SCT"))
assay <- if (params$normalisation_method == "SCT") "SCT" else params$rna_assay
seu <- RunPCA(object = seu, assay = assay, reduction.name = "pca", reduction.key = "PC_", verbose = FALSE)
if (!is.null(params$integration_method)) {
  logger::log_info("Integrating RNA data across {length(unique(seu$batch))} batches")
  args <- list(
    object = "seu",
    assay = "assay",
    method = "match.fun(params$integration_method)",
    normalization.method = "params$normalisation_method",
    orig.reduction = "'pca'",
    new.reduction = "'integrated.rna'",
    features = "VariableFeatures(object = seu, assay = assay)",
    dims = "1:50",
    verbose = "!params$quiet"
  ) %>%
    paste(names(.), ., sep = " = ", collapse = ", ") %>%
    paste(., params$integration_args, sep = ", ")
  seu <- eval(expr = parse(text = glue::glue("IntegrateLayers({args})")))
}

# ATAC data
if (params$atac_assay %in% Assays(seu)) {
  logger::log_info("Performing LSI for ATAC data")
  assay <- params$atac_assay
  seu <- RunSVD(object = seu, assay = assay, reduction.name = "lsi", reduction.key = "LSI_", verbose = FALSE)
  if (!is.null(params$integration_method)) {
    logger::log_info("Integrating ATAC data across {length(unique(seu$batch))} batches")
    # Temporarily switch ATAC assay from ChromatinAssay to Assay v5 to enable use of IntegrateLayers
    suppressWarnings({
      tmp <- seu[[assay]]
      seu[[assay]] <- CreateAssay5Object(data = GetAssayData(object = seu[[assay]], slot = "data")) %>%
        split(f = seu$batch)
      LayerData(object = seu, assay = assay, layer = "scale.data") <- SparseEmptyMatrix(nrow = nrow(tmp), ncol = ncol(tmp), rownames = rownames(tmp), colnames = colnames(tmp))
    }) # suppress unhelpful warnings
    args <- list(
      object = "seu",
      assay = "assay",
      method = "match.fun(gsub(pattern = 'PCA', replacement = 'LSI', x = params$integration_method))", # use custom LSI variant of integration function for ATAC data (if applicable)
      orig.reduction = "'lsi'",
      new.reduction = "'integrated.atac'",
      features = "VariableFeatures(object = tmp)",
      dims = "2:50",
      verbose = "!params$quiet"
    ) %>%
      paste(names(.), ., sep = " = ", collapse = ", ") %>%
      paste(., params$integration_args, sep = ", ")
    seu <- eval(expr = parse(text = glue::glue("IntegrateLayers({args})"))) %>% suppressWarnings() # suppress unhelpful warnings
    suppressWarnings({
      seu[[assay]] <- tmp
    }) # suppress unhelpful warnings
  }
}

# ADT data
if (params$adt_assay %in% Assays(seu)) {
  logger::log_info("Performing PCA for ADT data")
  assay <- params$adt_assay
  seu <- RunPCA(object = seu, assay = assay, reduction.name = "apca", reduction.key = "APC_", approx = FALSE, verbose = FALSE)
  if (!is.null(params$integration_method)) {
    logger::log_info("Integrating ADT data across {length(unique(seu$batch))} batches")
    args <- list(
      object = "seu",
      assay = "assay",
      method = "match.fun(params$integration_method)",
      orig.reduction = "'apca'",
      new.reduction = "'integrated.adt'",
      features = "VariableFeatures(object = seu, assay = assay)",
      dims = "1:50",
      verbose = "!params$quiet"
    ) %>%
      paste(names(.), ., sep = " = ", collapse = ", ") %>%
      paste(., params$integration_args, sep = ", ")
    seu <- eval(expr = parse(text = glue::glue("IntegrateLayers({args})"))) %>% suppressWarnings() # suppress unhelpful warnings
    seu[[paste0(assay, "C")]] <- CreateAssayObject(data = t(Embeddings(seu, reduction = "integrated.adt") %*% t(Loadings(seu, reduction = "integrated.adt"))))
    seu <- ScaleData(object = seu, assay = paste0(assay, "C"))
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
