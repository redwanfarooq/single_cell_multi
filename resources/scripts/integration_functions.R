# ==============================
# Integration helper functions
# ==============================


# Explicitly define set if null operator
`%||%` <- rlang::`%||%`


# Customised general wrapper functions for integration methods
HarmonyIntegration <- function(object,
                               assay = NULL,
                               layers = NULL,
                               orig = NULL,
                               new.reduction = "integrated.dr",
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
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  var.args <- list(...)
  if (!is.null(var.args$meta_data) && !is.null(var.args$vars_use)) {
    vars_use <- var.args$vars_use
    groups <- var.args$meta_data[colnames(object), vars_use, drop = FALSE]
  } else {
    if (!is.null(var.args$meta_data) || !is.null(var.args$vars_use)) warning("Arguments 'meta_data' and 'vars_use' are not both provided; using Seurat::CreateIntegrationGroups() instead.")
    vars_use <- "group"
    groups <- Seurat:::CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
  }
  var.args$meta_data <- var.args$vars_use <- NULL

  if (verbose) message("Running Harmony")
  out <- do.call(
    what = harmony::RunHarmony,
    args = c(
      list(
        data_mat = SeuratObject::Embeddings(object = orig)[, dims],
        meta_data = groups,
        vars_use = vars_use,
        verbose = verbose
      ),
      var.args
    ),
  )
  merged <- SeuratObject::CreateDimReducObject(
    embeddings = out,
    loadings = SeuratObject::Loadings(object = orig)[, dims],
    assay = assay,
    key = paste0(gsub(pattern = "[^[:alpha:]]", replacement = "", x = new.reduction), "_")
  ) %>%
    suppressWarnings()
  output.list <- list(merged)
  names(output.list) <- c(new.reduction)
  return(output.list)
}
attr(x = HarmonyIntegration, which = "Seurat.method") <- "integration"


FastMNNIntegration <- function(object,
                               assay = NULL,
                               layers = NULL,
                               orig = NULL,
                               new.reduction = "integrated.dr",
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
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  groups <- Seurat:::CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)

  if (verbose) message("Running reducedMNN")
  out <- do.call(
    what = batchelor::reducedMNN,
    args = list(
      SeuratObject::Embeddings(object = orig)[, dims],
      batch = groups[, 1],
      ...
    ),
  )
  merged <- SeuratObject::CreateDimReducObject(
    embeddings = out$corrected,
    loadings = SeuratObject::Loadings(object = orig)[, dims],
    assay = assay,
    key = paste0(gsub(pattern = "[^[:alpha:]]", replacement = "", x = new.reduction), "_")
  ) %>%
    suppressWarnings()
  output.list <- list(merged)
  names(output.list) <- c(new.reduction)
  return(output.list)
}
attr(x = FastMNNIntegration, which = "Seurat.method") <- "integration"


scVIIntegration <- function(object,
                            assay = NULL,
                            layers = NULL,
                            orig = NULL,
                            new.reduction = "integrated.dr",
                            features = NULL,
                            normalization.method = c("LogNormalize", "SCT"),
                            dims = 1:30,
                            scale.layer = "scale.data",
                            empirical_protein_background_prior = NULL,
                            max_epochs = NULL,
                            verbose = TRUE,
                            ...) {
  features <- features %||% Seurat::SelectIntegrationFeatures5(object = object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  var.args <- list(...)
  if (!is.null(var.args$meta_data) && !is.null(var.args$batch_key)) {
    args.setup <- list(batch_key = var.args$batch_key)
    if (!is.null(var.args$categorical_covariate_keys) || !is.null(var.args$continuous_covariate_keys)) {
      args.setup$categorical_covariate_keys <- as.list(var.args$categorical_covariate_keys)
      args.setup$continuous_covariate_keys <- as.list(var.args$continuous_covariate_keys)
    }
    vars_use <- c(var.args$batch_key, var.args$categorical_covariate_keys, var.args$continuous_covariate_keys) |> unlist()
    groups <- var.args$meta_data[colnames(object), vars_use, drop = FALSE]
  } else {
    if (!is.null(var.args$meta_data) || !is.null(var.args$batch_key)) warning("Arguments 'meta_data' and 'batch_key' are not both provided; using Seurat::CreateIntegrationGroups() instead.")
    args.setup <- list(batch_key = "group")
    groups <- Seurat:::CreateIntegrationGroups(object, layers = layers, scale.layer = scale.layer)
  }
  var.args$meta_data <- var.args$batch_key <- var.args$categorical_covariate_keys <- var.args$continuous_covariate_keys <- NULL

  if (verbose) message("Loading Python libraries")
  ad <- reticulate::import("anndata", convert = FALSE)
  scipy <- reticulate::import("scipy", convert = FALSE)
  scvi <- reticulate::import("scvi", convert = FALSE)

  # infer modality from assay name to select correct scVI model
  model <- dplyr::case_when(
    grepl(pattern = "rna|sct|gene", x = assay, ignore.case = TRUE) ~ "SCVI",
    grepl(pattern = "atac|peak|chrom", x = assay, ignore.case = TRUE) ~ "PEAKVI",
    grepl(pattern = "adt|cite|prot", x = assay, ignore.case = TRUE) ~ "TOTALVI",
    TRUE ~ NA
  )
  if (is.na(model)) stop("Unable to infer modality from assay name: ", assay)
  # extract user-defined model parameters accounting for parametrisation differences
  args.model <- var.args
  for (param in c("n_hidden", "n_layers", "n_latent")) {
    if (!is.null(args.model[[param]])) {
      args.model[[param]] <- as.integer(args.model[[param]])
    }
  }
  if (model %in% c("PEAKVI", "TOTALVI")) {
    if (!is.null(args.model$n_layers)) {
      args.model$n_layers_encoder <- args.model$n_layers_decoder <- args.model$n_layers
      args.model$n_layers <- NULL
    }
    if (model == "TOTALVI" && !is.null(args.model$dropout_rate)) {
      args.model$dropout_rate_encoder <- args.model$dropout_rate_decoder <- args.model$dropout_rate
      args.model$dropout_rate <- NULL
    }
  }

  # extract counts, integration groups and feature metadata and create AnnData object
  if (is(object, "Assay5")) object <- SeuratObject::JoinLayers(object = object, layers = "counts")
  args.adata <- list(obs = groups, var = NULL)
  if (model == "TOTALVI") {
    # create a dummy gene count matrix
    args.adata$X <- scipy$sparse$csr_matrix(SeuratObject::SparseEmptyMatrix(nrow = ncol(object), ncol = 1, rownames = colnames(object), colnames = "XXX"))
    # get protein count matrix
    args.adata$obsm <- list(protein = as.data.frame(Matrix::t(LayerData(object = object, layer = "counts")[features, ])))
    args.setup$protein_expression_obsm_key <- "protein"
    args.model$empirical_protein_background_prior <- if (is.null(empirical_protein_background_prior)) NULL else as.logical(empirical_protein_background_prior)
  } else {
    # get gene or peak count matrix
    args.adata$X <- scipy$sparse$csr_matrix(Matrix::t(LayerData(object = object, layer = "counts")[features, ]))
  }
  adata <- do.call(ad$AnnData, args.adata)
  args.setup$adata <- args.model$adata <- adata

  # set up and train the model
  if (verbose) message("Running ", model)
  do.call(scvi$model[[model]]$setup_anndata, args.setup)
  vae <- do.call(scvi$model[[model]], args.model)
  vae$train(max_epochs = if (is.null(max_epochs)) NULL else as.integer(max_epochs))

  # extract the latent embedding
  latent <- vae$get_latent_representation() |> as.matrix()
  rownames(latent) <- reticulate::py_to_r(adata$obs$index$values)
  merged <- SeuratObject::CreateDimReducObject(
    embeddings = latent,
    assay = assay,
    key = paste0(gsub(pattern = "[^[:alpha:]]", replacement = "", x = new.reduction), "_")
  ) %>%
    suppressWarnings()
  output.list <- list(merged)
  names(output.list) <- c(new.reduction)

  # extract batch-corrected log-normalised counts for ADT data (useful for visualisation and gating populations using biplots)
  # must be obtained directly from the trained model as the latent embedding does not have loading scores for the original features in contrast to PCA-based batch correction methods
  if (model == "TOTALVI") {
    data <- reticulate::py_to_r(vae$get_normalized_expression(transform_batch = list(unique(groups[, 1, drop = TRUE])), n_samples = 25L, return_mean = TRUE))[[2]]
    data <- data |>
      t() |>
      as("CsparseMatrix") |>
      log1p()
    colnames(data) <- colnames(object)
    rownames(data) <- features
    output.list[[paste0(assay, "C")]] <- SeuratObject::CreateAssayObject(data = data)
  }

  return(output.list)
}
attr(x = scVIIntegration, which = "Seurat.method") <- "integration"


# Customised ChromatinAssay-specific wrapper functions for integration methods
RLSIIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = "integrated.dr",
    reference = NULL,
    features = NULL,
    normalization.method = "LogNormalize",
    dims = 2:30,
    k.filter = NA,
    scale.layer = NULL,
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...) {
  rlang::check_installed(
    pkg = "Signac",
    reason = "for running integration with RLSIIntegration"
  )
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  features <- features %||% Seurat::SelectIntegrationFeatures5(object = object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  # check that there enough cells present
  ncells <- sapply(
    X = layers,
    FUN = function(x) {
      ncell <- dim(object[x])[2]
      return(ncell)
    }
  )
  if (min(ncells) < max(dims)) {
    abort(message = "At least one layer has fewer cells than dimensions specified, please lower 'dims' accordingly.")
  }
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- suppressMessages(suppressWarnings(
      SeuratObject::CreateSeuratObject(counts = NULL, data = object[layers[i]][features, ])
    ))
    SeuratObject::VariableFeatures(object = object.list[[i]]) <- features
    object.list[[i]] <- Signac::RunSVD(object = object.list[[i]], verbose = FALSE, n = max(dims))
  }
  anchor <- Seurat::FindIntegrationAnchors(
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
      suppressWarnings(expr = x <- Seurat::DietSeurat(x, features = features[1:2]))
      return(x)
    }
  )
  object_merged <- Seurat::IntegrateEmbeddings(
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


JointLSIIntegration <- function(
    object = NULL,
    assay = NULL,
    layers = NULL,
    orig = NULL,
    new.reduction = "integrated.dr",
    reference = NULL,
    features = NULL,
    normalization.method = "LogNormalize",
    dims = 2:30,
    k.anchor = 20,
    scale.layer = "scale.data",
    dims.to.integrate = NULL,
    k.weight = 100,
    weight.reduction = NULL,
    sd.weight = 1,
    sample.tree = NULL,
    preserve.order = FALSE,
    verbose = TRUE,
    ...) {
  op <- options(Seurat.object.assay.version = "v3", Seurat.object.assay.calcn = FALSE)
  on.exit(expr = options(op), add = TRUE)
  features <- features %||% Seurat::SelectIntegrationFeatures5(object = object)
  features.diet <- features[1:2]
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layers <- layers %||% SeuratObject::Layers(object = object, search = "data")
  object.list <- list()
  for (i in seq_along(along.with = layers)) {
    object.list[[i]] <- CreateSeuratObject(counts = object[layers[i]][features.diet, ])
    object.list[[i]][["RNA"]]$counts <- NULL
    object.list[[i]][["joint.pca"]] <- CreateDimReducObject(
      embeddings = Embeddings(object = orig)[Cells(object.list[[i]]), ],
      assay = "RNA",
      loadings = Loadings(orig),
      key = "J_"
    )
  }
  anchor <- FindIntegrationAnchors(
    object.list = object.list,
    anchor.features = features.diet,
    scale = FALSE,
    reduction = "jpca",
    normalization.method = normalization.method,
    dims = dims,
    k.anchor = k.anchor,
    k.filter = NA,
    reference = reference,
    verbose = verbose,
    ...
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
attr(x = JointLSIIntegration, which = "Seurat.method") <- "integration"


# Function to build sample tree matrix for hierarchical integration using Seurat anchor-based methods
build.sample.tree <- function(groups, init = 0, subtree = FALSE) {
  if (!is.list(groups)) stop("Argument <groups>: must be a list")
  if (length(unique(sapply(groups, typeof))) > 1) stop("Argument <groups>: all elements must be of the same type (either lists or integer vectors)")

  sample.trees <- vector(mode = "list", length = length(groups))

  # Build sample tree matrix per group
  for (i in seq_along(groups)) {
    if (is.list(groups[[i]])) {
      init <- suppressWarnings(max(sapply(sample.trees, max, na.rm = TRUE) + 1, init))
      sample.trees[[i]] <- build.sample.tree(groups[[i]], init, subtree = TRUE)
    } else {
      steps <- vector(mode = "list", length = length(groups[[i]]) - 1)
      for (j in seq_len(length(groups[[i]]) - 1)) {
        if (j == 1) {
          steps[[j]] <- c(-groups[[i]][j], -groups[[i]][j + 1])
        } else {
          steps[[j]] <- c(init + j - 1, -groups[[i]][j + 1])
        }
      }
      sample.trees[[i]] <- do.call(rbind, steps)
      attr(sample.trees[[i]], "init") <- init
      attr(sample.trees[[i]], "subtree") <- FALSE
    }
  }

  # Build merge tree matrix across groups
  breaks <- vector(mode = "integer", length = length(groups))
  for (i in seq_along(sample.trees)) {
    if (i > 1 && !attr(sample.trees[[i]], "subtree")) {
      tmp <- sample.trees[[i]]
      sample.trees[[i]] <- t(apply(X = sample.trees[[i]], MARGIN = 1, FUN = function(x) ifelse(x > 0, x + sum(sapply(seq_len(i - 1), function(x) nrow(sample.trees[[x]]))), x)))
      attr(sample.trees[[i]], "init") <- attr(tmp, "init")
      attr(sample.trees[[i]], "subtree") <- attr(tmp, "subtree")
    }
    breaks[i] <- nrow(sample.trees[[i]])
  }
  breaks <- cumsum(breaks) + sapply(sample.trees, function(x) if (!attr(x, "subtree")) attr(x, "init") else 0)
  steps <- vector(mode = "list", length = length(breaks) - 1)
  for (i in seq_len(length(breaks) - 1)) {
    if (i == 1) {
      steps[[i]] <- c(breaks[i], breaks[i + 1])
      step <- max(breaks) + 1
    } else {
      steps[[i]] <- c(step, breaks[i + 1])
      step <- step + 1
    }
  }
  merge.trees <- do.call(rbind, steps)

  # Combine sample and merge tree matrices
  final.tree <- do.call(rbind, sample.trees) %>% rbind(merge.trees)
  attr(final.tree, "init") <- init
  attr(final.tree, "subtree") <- subtree
  return(final.tree)
}
