#!/bin/env -S Rscript --vanilla


# ==============================
# COMMAND LINE OPTIONS
# ==============================
# Define options
DOC <- "
Merge counts and metadata from multisample multimodal single cell experiments

Usage:
  merge.R --samples=<sample[;sample...]> --output=<output> --hdf5=<path[;path...]> --metadata=<path[;path...]> [--fragments=<path[;path...]>] [--summits=<path[;path...]>] [options]

Arguments:
  REQUIRED
  --samples=<sample[;sample...]>          Semicolon-separated list of sample IDs
  --output=<output>                       Path to output file (must include file extension 'qs' or 'rds')

  --hdf5=<path[;path...]>                 Path(s) to 10x-formatted HDF5 files containing combined multimodal count matrices
  --metadata=<path[;path...]>             Path(s) to TSV files containing cell metadata
  --fragments=<path[;path...]>            Path(s) to ATAC fragments files (if applicable)
  --summits=<path[;path...]>              Path(s) to ATAC peak summit BED files from MACS2 (if applicable)

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
  --peak-method=<method>                  Method to select joint ATAC peaks across samples ('bulk', 'fixed', 'disjoin' or 'reduce')
  --filter-cells=<str>                    Logical expression to filter cells based on metadata columns (matching cells will be discarded)


Options:
  -h --help                               Show this screen
  -b --binarize                           Binarize ATAC counts
  -d --downsample                         Downsample counts to match sequencing depth across libraries
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
for (x in c("hdf5", "metadata", "fragments", "summits")) {
  if (!is.null(params[[x]])) {
    if (length(params[[x]]) == 1 && grepl(pattern = "\\{sample\\}", x = params[[x]])) params[[x]] <- glue::glue(params[[x]], sample = params$samples) # replace placeholder(s) if path provided is template string
    if (length(params[[x]]) != length(params$samples)) stop("Argument <", x, ">: number of paths must be equal to number of samples or be a valid template string using '{sample}' as a placeholder for sample ID")
  }
}
if (!is.null(params$fragments)) {
  names(params$fragments) <- params$samples
  missing.files <- vapply(params$fragments, function(x) !file.exists(x), logical(1))
  for (file in params$fragments[missing.files]) {
    warning("File not found: ", file)
  }
  params$fragments <- params$fragments[!missing.files]
}
if (!is.null(params$summits)) {
  names(params$summits) <- params$samples
  missing.files <- vapply(params$summits, function(x) !file.exists(x), logical(1))
  for (file in params$summits[missing.files]) {
    warning("File not found: ", file)
  }
  params$summits <- params$summits[!missing.files]
}
if (!is.null(params$peak_method) && params$peak_method == "fixed" && is.null(params$summits)) stop("Argument <summits>: must be provided if running with '--peak-method=fixed'")
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

source("utils.R")
source("iterative_overlap_peak_merging.R")

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
    x <- get.10x.h5(file)
    if (!is.list(x)) x <- list(x) # convert to list if only one modality is present
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
if (!is.null(mat$RNA)) {
  logger::log_info("Cleaning up feature names")
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  res <- future_map(
    .x = mat$RNA,
    .f = function(x, gene.types = params$gene_types, db = ensdb) {
      df <- data.frame(rownames = rownames(x))
      df$ensembl_id <- sub(pattern = "[.][0-9]*", replacement = "", x = df$rownames)
      mapping <- ensembldb::select(
        db,
        keys = df$ensembl_id,
        keytype = "GENEID",
        columns = c("GENEID", "SYMBOL", "GENEBIOTYPE")
      ) %>%
        dplyr::rename(
          ensembl_id = GENEID,
          gene_symbol = SYMBOL,
          gene_type = GENEBIOTYPE
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
}
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
  # Create fragments objects
  fragments <- future_pmap(
    .l = list(
      x = mat$ATAC,
      sample = names(mat$ATAC),
      path = params$fragments[names(mat$ATAC)]
    ),
    .f = function(x, sample, path) {
      if (!params$quiet) message("Getting fragments for ", sample)
      CreateFragmentObject(
        path = path,
        cells = colnames(x),
        verbose = FALSE
      )
    },
    .options = furrr.options
  ) %>%
    setNames(names(mat$ATAC))
  if (!is.null(params$peak_method) && length(mat$ATAC) > 1) {
    # Get joint peak set across samples
    logger::log_info("Finding and quantifying joint ATAC peak set")
    params$peak_method <- params$peak_method %>% match.arg(choices = c("bulk", "fixed", "disjoin", "reduce"))
    if (params$peak_method == "bulk") {
      # downsample fragments to match depth (median unique fragments per cell) across samples
      cells <- map(.x = mat$ATAC, .f = colnames)
      depth <- future_map2_dbl(
        .x = params$fragments,
        .y = cells,
        .f = function(file, cells) {
          metrics <- CountFragments(fragments = file, cells = cells, verbose = FALSE)
          return(median(metrics$frequency_count))
        },
        .options = furrr.options
      )
      downsampled.fragments <- future_pmap_chr(
        .l = list(
          file = params$fragments,
          cells = cells,
          proportion = min(depth) / depth
        ),
        .f = function(file, cells, proportion) {
          if (proportion < 1) {
            out <- downsample.fragments(
              input = file,
              cells = cells,
              proportion = proportion
            )
          } else {
            out <- file
          }
          return(out)
        },
        .options = furrr.options
      )
      # call bulk peaks from downsampled fragments using MACS2
      if (!params$quiet) message("Calling joint peaks from ", length(downsampled.fragments), " fragments files using MACS2")
      peaks <- CallPeaks(
        object = downsampled.fragments,
        format = "BED",
        extsize = 200,
        shift = -100,
        additional.args = "--nolambda --max-gap 0 --keep-dup all",
        cleanup = TRUE,
        verbose = FALSE
      )
    } else if (params$peak_method == "fixed") {
      # get fixed width peaks for each sample and combine using iterative overlap merging algorithm
      # see Corces, M.R. et al. (2018). The chromatin accessibility landscape of primary human cancers. Science, 362, eaav1898.
      genome <- GenomeInfoDb::Seqinfo(genome = "hg38") %>%
        GenomicRanges::GRanges(seqnames = names(.), ranges = IRanges::IRanges(start = 1, end = seqlengths(.)), seqinfo = .) %>%
        GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")
      peaks <- future_map2(
        .x = params$summits,
        .y = names(params$summits),
        .f = function(file, sample, valid.granges = genome) {
          if (!params$quiet) message("Getting fixed width peaks for ", sample)
          peaks <- readSummits(file) %>%
            IRanges::resize(width = 501, fix = "center") %>%
            IRanges::subsetByOverlaps(valid.granges, type = "within") %>% # remove peaks that extend past chromosome boundaries
            nonOverlappingGRanges(by = "score", decreasing = TRUE, verbose = FALSE) %>%
            GenomeInfoDb::sortSeqlevels() %>%
            sort()
          S4Vectors::mcols(peaks)$score <- round(getQuantiles(S4Vectors::mcols(peaks)$score), 3) # normalise peak scores
          return(peaks)
        },
        .options = furrr.options
      ) %>%
        unname()
      peaks <- peaks %>%
        GenomicRanges::GRangesList() %>%
        unlist() %>%
        nonOverlappingGRanges(by = "score", decreasing = TRUE, verbose = TRUE) %>%
        GenomeInfoDb::sortSeqlevels() %>%
        sort()
    } else {
      # merge overlapping peaks across samples
      if (!params$quiet) message("Merging peaks across samples with method '", params$peak_method, "'")
      merge.function <- if (params$peak_method == "disjoin") GenomicRanges::disjoin else GenomicRanges::reduce
      peaks <- future_map(
        .x = mat$ATAC,
        .f = function(x) {
          stringr::str_split(string = rownames(x), pattern = "[:-]", simplify = TRUE) %>%
            as.data.frame() %>%
            setNames(c("chr", "start", "stop")) %>%
            GenomicRanges::makeGRangesFromDataFrame()
        },
        .options = furrr.options
      ) %>%
        Reduce(f = c, x = .) %>%
        merge.function() %>%
        suppressWarnings() # suppress unhelpful warnings
    }
    peaks <- peaks %>%
      IRanges::subsetByOverlaps(blacklist_hg38_unified, invert = TRUE) %>% # exclude peaks overlapping genomic blacklist regions
      GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") # remove peaks on non-standard chromosomes
    # count fragments in merged peaks
    mat$ATAC <- pmap(
      .l = list(
        x = mat$ATAC,
        sample = names(mat$ATAC),
        fragments.object = fragments[names(mat$ATAC)]
      ),
      .f = function(x, sample, fragments.object, joint.peaks = peaks) {
        if (!params$quiet) message("Counting fragments in joint peaks for ", sample)
        FeatureMatrix(
          fragments = fragments.object,
          features = joint.peaks,
          cells = colnames(x),
          sep = c(":", "-"),
          verbose = !params$quiet
        )
      }
    )
  }
}
if (params$downsample && length(params$samples) > 1) {
  # Downsample counts to match sequencing depth across libraries
  mat <- map2(
    .x = mat,
    .y = names(mat),
    .f = function(x, type) {
      set.seed(42)
      logger::log_info("Downsampling {type} counts to match sequencing depth across {length(x)} libraries")
      x <- as.list(DropletUtils::downsampleBatches(x))
      return(x)
    }
  )
}
if (params$binarize && !is.null(mat$ATAC)) {
  # Binarize ATAC counts
  logger::log_info("Binarizing ATAC counts")
  mat$ATAC <- future_map(
    .x = mat$ATAC,
    .f = BinarizeCounts,
    .options = furrr.options
  )
}
mat <- transpose(mat)


# Create merged Seurat object
logger::log_info("Creating merged Seurat object")
seu <- future_pmap(
  .l = list(
    x = mat,
    sample = names(mat),
    cell.metadata = cell.metadata,
    gene.metadata = if (exists("gene.metadata")) gene.metadata else map(.x = seq_along(mat), .f = function(x) NULL)
  ),
  function(x, sample, cell.metadata, gene.metadata, fragments.objects = eval(if (exists("fragments")) fragments else NULL)) {
    assays <- map2(
      .x = x,
      .y = names(x),
      .f = function(counts, modality, fragment.object = fragments.objects[[sample]]) {
        switch(modality,
          RNA = CreateAssay5Object(counts = counts),
          ATAC = CreateChromatinAssay(
            counts = counts,
            sep = c(":", "-"),
            fragments = fragment.object
          ),
          ADT = CreateAssay5Object(counts = counts)
        )
      }
    )
    if ("RNA" %in% names(assays)) assays$RNA <- AddMetaData(object = assays$RNA, metadata = gene.metadata)
    obj <- CreateSeuratObject(counts = assays[[1]], assay = names(assays)[1], project = sample, meta.data = cell.metadata)
    if (length(assays) > 1) {
      for (i in seq.int(from = 2, to = length(assays))) {
        obj[[names(assays)[i]]] <- assays[[i]]
      }
    }
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

# Filter cells (if specified)
if (!is.null(params$filter_cells)) {
  logger::log_info("Filtering cells")
  discard <- with(seu[[]], eval(parse(text = params$filter)))
  if (!is.logical(discard)) stop("Filter must return a logical vector; check '--filter-cells' argument")
  if (sum(discard) == 0) {
    warning("No cells matching filter")
  } else if (sum(!discard) == 0) {
    stop("All cells removed by filter; check '--filter-cells' argument")
  } else {
    seu <- seu[, !discard]
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
