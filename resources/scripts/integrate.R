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
  --integration-method=<method>           Integration method function name [default: RPCAIntegration]
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
  library(BPCells)
})

options(future.globals.maxSize = 50000 * 1024^2)

if (!dir.exists(dirname(params$output))) dir.create(dirname(params$output), recursive = TRUE)

# Custom integration function for RLSI (RPCA for ATAC data)
RLSIIntegration <- function(
  object = NULL,
  assay = NULL,
  layers = NULL,
  orig = NULL,
  new.reduction = "integrated.dr",
  reference = NULL,
  features = NULL,
  normalization.method = "LogNormalize",
  dims = 1:30,
  k.filter = NA,
  scale.layer = NULL,
  dims.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE,
  verbose = TRUE,
  ...
) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  features <- features %||% SelectIntegrationFeatures5(object = object)
  assay <- assay %||% "ATAC"
  layers <- layers %||% Layers(object = object, search = "data")
  # check that there enough cells present
  ncells <- sapply(
    X = layers,
    FUN = function(x) {
      ncell <- dim(object[x])[2]
      return(ncell)
    }
  )
  if (min(ncells) < max(dims))  {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- suppressMessages(suppressWarnings(
      CreateSeuratObject(counts = NULL, data = object[layers[i]][features, ])
    ))
    VariableFeatures(object = object.list[[i]]) <- features
    object.list[[i]] <- RunSVD(object = object.list[[i]], verbose = FALSE, n = max(dims))
  }
  anchor <- FindIntegrationAnchors(
    object.list = object.list,
    anchor.features = features,
    scale = FALSE,
    reduction = "rlsi",
    normalization.method = normalization.method,
    dims = dims,
    k.filter = k.filter,
    reference = reference,
    verbose = verbose,
    ...
  )
  slot(object = anchor, name = "object.list") <- lapply(
    X = slot(object = anchor, name = "object.list"),
    FUN = function(x) {
      suppressWarnings(expr = x <- DietSeurat(x, features = features[1:2]))
      return(x)
    }
  )
  object_merged <- IntegrateEmbeddings(
    anchorset = anchor,
    reductions = orig,
    new.reduction.name = new.reduction,
    dims.to.integrate = dims.to.integrate,
    k.weight = k.weight,
    weight.reduction = weight.reduction,
    sd.weight = sd.weight,
    sample.tree = sample.tree,
    preserve.order = preserve.order,
    verbose = verbose
  )
  output.list <- list(object_merged[[new.reduction]])
  names(output.list) <- c(new.reduction)
  return(output.list)
}
attr(x = RLSIIntegration, which = "Seurat.method") <- "integration"

# Custom integration function for reducedMNN (fastMNN algorithm with precomputed dimensionality reduction)
ReducedMNNIntegration <- function(object,
                                  assay = NULL,
                                  layers = NULL,
                                  orig = NULL,
                                  new.reduction = "integrated.mnn",
                                  features = NULL,
                                  normalization.method = c("LogNormalize", "SCT"),
                                  dims = 1:30,
                                  scale.layer = "scale.data",
                                  verbose = TRUE,
                                  ...) {
  rlang::check_installed(
    pkg = "batchelor",
    reason = "for running integration with reducedMNN"
  )
  assay <- assay %||% DefaultAssay(object = object)
  layers <- layers %||% Layers(object = object, search = "data")
  groups <- Seurat:::CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)

  if (verbose) message("Running reducedMNN")
  out <- do.call(
    what = batchelor::reducedMNN,
    args = list(
      Embeddings(object = orig)[, dims],
      batch = groups[, 1],
      ...
    ),
  )
  merged <- CreateDimReducObject(
    embeddings = out$corrected,
    loadings = Loadings(object = orig)[, dims],
    assay = assay,
    key = paste0(gsub(pattern = "[^[:alpha:]]", replacement = "", x = new.reduction), "_")
  ) %>%
    suppressWarnings()
  output.list <- list(merged)
  names(output.list) <- c(new.reduction)
  return(output.list)
}
attr(x = ReducedMNNIntegration, which = "Seurat.method") <- "integration"

# Custom integration function for Harmony
HarmonyIntegration <- function(object,
                               assay = NULL,
                               layers = NULL,
                               orig = NULL,
                               new.reduction = "integrated.harmony",
                               features = NULL,
                               normalization.method = c("LogNormalize", "SCT"),
                               dims = 1:30,
                               scale.layer = "scale.data",
                               verbose = TRUE,
                               ...) {
  rlang::check_installed(
    pkg = "harmony",
    reason = "for running integration with Harmony"
  )
  assay <- assay %||% DefaultAssay(object = object)
  layers <- layers %||% Layers(object = object, search = "data")
  var.args <- list(...)
  if (!is.null(var.args$meta_data) && !is.null(var.args$vars_use)) {
    vars_use <- var.args$vars_use
    groups <- var.args$meta_data[, vars_use]
  } else {
    vars_use <- "group"
    groups <- Seurat:::CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
  }
  var.args$meta_data <- NULL
  var.args$vars_use <- NULL

  if (verbose) message("Running Harmony")
  out <- do.call(
    what = harmony::RunHarmony,
    args = c(
      list(
        data_mat = Embeddings(object = orig)[, dims],
        meta_data = groups,
        vars_use = vars_use,
        verbose = verbose
      ),
      var.args
    ),
  )
  merged <- CreateDimReducObject(
    embeddings = out,
    loadings = Loadings(object = orig)[, dims],
    assay = assay,
    key = paste0(gsub(pattern = "[^[:alpha:]]", replacement = "", x = new.reduction), "_")
  ) %>%
    suppressWarnings()
  output.list <- list(merged)
  names(output.list) <- c(new.reduction)
  return(output.list)
}
attr(x = HarmonyIntegration, which = "Seurat.method") <- "integration"

# Function to build sample tree matrix for hierarchical integration
build.sample.tree <- function(groups) {
  ngroups <- length(groups)
  nsamples <- map_int(groups, length)

  steps <- vector(mode = "list", length = ngroups + sum(nsamples - 1) - 1)
  breaks <- vector(mode = "integer", length = ngroups)

  step <- 1
  for (i in seq_len(ngroups)) {
    for (j in seq_len(nsamples[i] - 1)) {
      if (j == 1) {
        steps[[step]] <- c(-groups[[i]][j], -groups[[i]][j + 1])
      } else {
        steps[[step]] <- c(step - 1, -groups[[i]][j + 1])
      }
      step <- step + 1
    }
    breaks[i] <- step - 1
  }
  for (i in seq_len(ngroups - 1)) {
    if (i == 1) {
      steps[[step]] <- c(breaks[i], breaks[i + 1])
    } else {
      steps[[step]] <- c(step - 1, breaks[i + 1])
    }
    step <- step + 1
  }

  tree <- do.call(rbind, steps)

  return(tree)
}


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


# Perform integration steps
# Integrate RNA data across batches
logger::log_info("Integrating RNA data across {length(unique(seu$batch))} batches")
params$normalisation_method <- params$normalisation_method %>% match.arg(choices = c("LogNormalize", "SCT"))
assay <- if (params$normalisation_method == "SCT") "SCT" else params$rna_assay
set.seed(42)
if (!params$quiet) message("Running multi-batch PCA")
pca <- batchelor::multiBatchPCA(
  eval({
    if (assay == "SCT") Seurat:::PrepDR(object = seu[[assay]], verbose = FALSE) else Seurat:::PrepDR5(object = seu[[assay]], verbose = FALSE)
  }),
  batch = seu$batch,
  get.variance = TRUE,
  preserve.single = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(workers = as.integer(params$threads))
)
seu[["pca"]] <- CreateDimReducObject(
  embeddings = pca[[1]],
  loadings = S4Vectors::metadata(pca)$rotation,
  stdev = sqrt(S4Vectors::metadata(pca)$var.explained),
  assay = assay,
  key = "PC_"
) %>%
  suppressWarnings() # suppress unhelpful warnings
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

# Integrate ATAC data across batches
if (params$atac_assay %in% Assays(seu)) {
  logger::log_info("Integrating ATAC data across {length(unique(seu$batch))} batches")
  assay <- params$atac_assay
  set.seed(42)
  if (!params$quiet) message("Running multi-batch PCA")
  lsi <- batchelor::multiBatchPCA(
    Seurat:::PrepDR(object = seu[[assay]], slot = "data", verbose = FALSE),
    batch = seu$batch,
    get.variance = TRUE,
    preserve.single = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = as.integer(params$threads))
  )
  seu[["lsi"]] <- CreateDimReducObject(
    embeddings = lsi[[1]],
    loadings = S4Vectors::metadata(lsi)$rotation,
    stdev = sqrt(S4Vectors::metadata(lsi)$var.explained),
    assay = assay,
    key = "LSI_"
  ) %>%
    suppressWarnings() # suppress unhelpful warnings
  # Temporarily switch ATAC assay from ChromatinAssay to Assay v5 to enable use of IntegrateLayers
  suppressWarnings({
    tmp <- seu[[assay]]
    seu[[assay]] <- CreateAssay5Object(data = GetAssayData(object = seu[[assay]], slot = "data")) %>%
      split(f = seu$batch)
    LayerData(object = seu, assay = assay, layer = "scale.data") <- SparseEmptyMatrix(nrow = nrow(seu[[assay]]), ncol = ncol(seu[[assay]]), rownames = rownames(seu[[assay]]), colnames = colnames(seu[[assay]]))
  }) # suppress unhelpful warnings
  args <- list(
    object = "seu",
    assay = "assay",
    method = "match.fun(gsub(pattern = 'PCA', replacement = 'LSI', x = params$integration_method))", # use custom LSI variant of integration function for ATAC data (if applicable)
    orig.reduction = "'lsi'",
    new.reduction = "'integrated.atac'",
    features = "VariableFeatures(object = tmp)",
    dims = "1:50",
    verbose = "!params$quiet"
  ) %>%
    paste(names(.), ., sep = " = ", collapse = ", ") %>%
    paste(., params$integration_args, sep = ", ")
  seu <- eval(expr = parse(text = glue::glue("IntegrateLayers({args})"))) %>% suppressWarnings() # suppress unhelpful warnings
  suppressWarnings({
    seu[[assay]] <- tmp
  }) # suppress unhelpful warnings
}

# Integrate ADT data across batches
if (params$adt_assay %in% Assays(seu)) {
  logger::log_info("Integrating ADT data across {length(unique(seu$batch))} batches")
  assay <- params$adt_assay
  set.seed(42)
  if (!params$quiet) message("Running multi-batch PCA")
  apca <- batchelor::multiBatchPCA(
    Seurat:::PrepDR5(object = seu[[assay]], verbose = FALSE),
    batch = seu$batch,
    get.variance = TRUE,
    preserve.single = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers = as.integer(params$threads)),
    BSPARAM = BiocSingular::ExactParam() # compute exact SVD for ADT assay due to low number of features
  )
  seu[["apca"]] <- CreateDimReducObject(
    embeddings = apca[[1]],
    loadings = S4Vectors::metadata(apca)$rotation,
    stdev = sqrt(S4Vectors::metadata(apca)$var.explained),
    assay = assay,
    key = "APC_"
  ) %>%
    suppressWarnings() # suppress unhelpful warnings
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
  seu <- ScaleData(object = seu, assay = paste0(assay, "C"), do.center = FALSE, do.scale = TRUE)
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
