# ==============================
# Integration helper functions
# ==============================


# Explicitly define set if null operator
`%||%` <- SeuratObject::`%||%`


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
  if (min(ncells) < max(dims))  {
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
