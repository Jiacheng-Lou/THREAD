# scrna.R
# -----------------------------------------------------------------------------
# Utilities for converting a quality-controlled single-cell expression object
# into the gene-by-label expression matrix used by DPM.
#
# Design principle:
#   DPM itself requires only a gene x label matrix X. This file provides a
#   lightweight helper for users who start from a Seurat object or from a
#   gene x cell matrix plus cell metadata. It does not attempt to perform a full
#   single-cell QC workflow, batch correction, integration, or cell annotation.
# -----------------------------------------------------------------------------


#' Preprocess single-cell expression into a DPM expression matrix
#'
#' @description
#' Converts a Seurat object or a gene-by-cell expression matrix into a
#' gene-by-label pseudobulk matrix. Labels can be cell types, subtypes, states,
#' clusters, domains, or any user-defined grouping variable. The returned matrix
#' has genes in rows and labels in columns.
#'
#' The default workflow is:
#' \enumerate{
#'   \item obtain normalized expression. If \code{already_normalized = FALSE},
#'     Seurat's \code{NormalizeData()} is used with \code{LogNormalize}; if
#'     \code{already_normalized = TRUE}, the input expression is used directly;
#'   \item average expression over all cells sharing the same label;
#'   \item optionally apply a second log transform to the pseudobulk means.
#' }
#'
#' This function is intentionally not a full scRNA-seq preprocessing pipeline.
#' Users should perform cell filtering, doublet removal, batch correction and
#' cell annotation upstream when needed. Users who already have a gene-by-label
#' matrix can pass it directly to \code{prepare_input()} and skip this function.
#'
#' @param object A Seurat object or a numeric matrix-like object with genes in
#'   rows and cells in columns. Sparse matrices from the \pkg{Matrix} package are
#'   supported for matrix input.
#' @param meta A data frame containing cell metadata. Required for matrix input;
#'   ignored for Seurat input unless supplied, in which case it overrides the
#'   object's metadata. Row names should match \code{colnames(object)}, unless
#'   \code{cell_col} is supplied.
#' @param subtype_col A single column name in the metadata containing the label
#'   to aggregate by, or a vector of labels with length equal to the number of
#'   cells.
#' @param cell_col Optional column name in \code{meta} containing cell IDs.
#' @param assay Seurat assay to use. If \code{NULL}, the object's default assay
#'   is used. Ignored for matrix input.
#' @param normalized_layer Assay layer/slot containing normalized expression when
#'   \code{already_normalized = TRUE}, and the layer/slot read after Seurat
#'   normalization when \code{already_normalized = FALSE}. Default \code{"data"}.
#' @param counts_layer Assay layer/slot containing raw counts for Seurat input
#'   when \code{already_normalized = FALSE}. Default \code{"counts"}.
#' @param already_normalized Logical. If \code{FALSE}, expression is treated as
#'   raw counts and normalized with Seurat's \code{NormalizeData()}. If
#'   \code{TRUE}, expression is treated as already normalized non-negative
#'   expression and is aggregated directly.
#' @param normalization_method Method passed to \code{Seurat::NormalizeData()}.
#'   Default \code{"LogNormalize"}.
#' @param scale_factor Scale factor passed to \code{Seurat::NormalizeData()}.
#'   Default \code{1e4}.
#' @param second_log Logical. Whether to apply a second log transform after
#'   label-level averaging. Default \code{TRUE}.
#' @param second_log_base Logarithm base for the second log transform. Supported
#'   values are \code{2} and \code{exp(1)}. Default \code{2}.
#' @param pseudocount Non-negative pseudocount used in the second log transform.
#'   Default \code{1}.
#' @param min_cells_per_subtype Minimum number of cells required for a label to
#'   be retained. Default \code{1}.
#' @param min_cells_per_gene Optional minimum number of cells in which a gene
#'   must have non-zero expression before aggregation. Default \code{0} disables
#'   this filter.
#' @param drop_na_subtype Logical. Whether to remove cells with missing labels.
#'   Default \code{TRUE}.
#' @param return_sparse Logical. Whether to return a sparse matrix when possible.
#'   The default \code{FALSE} returns a dense matrix because the number of label
#'   columns is usually small.
#' @param seurat_verbose Logical. Verbosity passed to Seurat normalization.
#'   Default \code{FALSE}.
#' @param verbose Logical. Whether to print concise progress messages. Default
#'   \code{TRUE}.
#'
#' @return A gene-by-label numeric matrix suitable for \code{prepare_input()}.
#'
#' @examples
#' \dontrun{
#' X <- preprocess_scrna(
#'   object = seurat_obj,
#'   subtype_col = "cell_subtype",
#'   already_normalized = FALSE
#' )
#' }
#'
#' @export
preprocess_scrna <- function(object,
                             meta = NULL,
                             subtype_col,
                             cell_col = NULL,
                             assay = NULL,
                             normalized_layer = "data",
                             counts_layer = "counts",
                             already_normalized = FALSE,
                             normalization_method = "LogNormalize",
                             scale_factor = 1e4,
                             second_log = TRUE,
                             second_log_base = 2,
                             pseudocount = 1,
                             min_cells_per_subtype = 1,
                             min_cells_per_gene = 0,
                             drop_na_subtype = TRUE,
                             return_sparse = FALSE,
                             seurat_verbose = FALSE,
                             verbose = TRUE) {
  validate_scrna_arguments(
    object = object,
    meta = meta,
    subtype_col = subtype_col,
    cell_col = cell_col,
    normalized_layer = normalized_layer,
    counts_layer = counts_layer,
    already_normalized = already_normalized,
    normalization_method = normalization_method,
    scale_factor = scale_factor,
    second_log = second_log,
    second_log_base = second_log_base,
    pseudocount = pseudocount,
    min_cells_per_subtype = min_cells_per_subtype,
    min_cells_per_gene = min_cells_per_gene,
    drop_na_subtype = drop_na_subtype,
    return_sparse = return_sparse,
    seurat_verbose = seurat_verbose,
    verbose = verbose
  )

  if (is_seurat_object(object)) {
    seurat_input <- extract_expression_from_seurat(
      object = object,
      assay = assay,
      normalized_layer = normalized_layer,
      counts_layer = counts_layer,
      already_normalized = already_normalized,
      normalization_method = normalization_method,
      scale_factor = scale_factor,
      seurat_verbose = seurat_verbose,
      verbose = verbose
    )
    expr <- seurat_input$expr
    meta_use <- if (is.null(meta)) seurat_input$meta else meta
  } else {
    if (is.null(meta)) {
      stop("'meta' is required when 'object' is a matrix-like expression object.")
    }
    validate_gene_cell_matrix(object, object_name = "object")
    if (!requireNamespace("Seurat", quietly = TRUE) && !isTRUE(already_normalized)) {
      stop(
        "Package 'Seurat' is required when already_normalized = FALSE. ",
        "Install Seurat or provide already normalized expression with already_normalized = TRUE."
      )
    }
    expr <- object
    meta_use <- meta

    if (!isTRUE(already_normalized)) {
      expr <- normalize_matrix_with_seurat(
        counts = expr,
        meta = meta_use,
        scale_factor = scale_factor,
        normalization_method = normalization_method,
        seurat_verbose = seurat_verbose
      )
      if (verbose) {
        THREAD_log("Input", "Applied Seurat::NormalizeData() to matrix input.")
      }
    } else {
      validate_nonnegative_matrix(expr, object_name = "object")
      if (verbose) {
        THREAD_log("Input", "Matrix input was treated as already normalized non-negative expression.")
      }
    }
  }

  validate_gene_cell_matrix(expr, object_name = "expression")

  labels <- extract_subtype_labels(
    cell_ids = colnames(expr),
    meta = meta_use,
    subtype_col = subtype_col,
    cell_col = cell_col,
    drop_na_subtype = drop_na_subtype,
    verbose = verbose
  )

  keep_cells <- !is.na(labels)
  if (!all(keep_cells)) {
    expr <- expr[, keep_cells, drop = FALSE]
    labels <- labels[keep_cells]
  }

  if (ncol(expr) == 0L) {
    stop("No cells remain after removing missing subtype labels.")
  }

  label_counts <- table(labels)
  keep_levels <- names(label_counts)[label_counts >= min_cells_per_subtype]
  if (length(keep_levels) == 0L) {
    stop("No label has at least min_cells_per_subtype cells.")
  }

  if (length(keep_levels) < length(label_counts)) {
    dropped <- setdiff(names(label_counts), keep_levels)
    if (verbose) {
      THREAD_log(
        "Input",
        "Dropping ", length(dropped),
        " label(s) with fewer than min_cells_per_subtype cells: ",
        paste(dropped, collapse = ", ")
      )
    }
    keep_cells <- labels %in% keep_levels
    expr <- expr[, keep_cells, drop = FALSE]
    labels <- labels[keep_cells]
  }

  labels <- factor(labels, levels = keep_levels)

  if (min_cells_per_gene > 0) {
    detected <- row_detected_counts(expr)
    keep_genes <- detected >= min_cells_per_gene
    if (!any(keep_genes)) {
      stop("No genes remain after applying min_cells_per_gene.")
    }
    if (verbose) {
      THREAD_log(
        "Input",
        "Keeping ", sum(keep_genes), " of ", length(keep_genes),
        " genes after min_cells_per_gene filtering."
      )
    }
    expr <- expr[keep_genes, , drop = FALSE]
  }

  validate_nonnegative_matrix(expr, object_name = "expression")
  X <- aggregate_expression_by_label(expr, labels = labels)

  if (isTRUE(second_log)) {
    X <- apply_second_log_transform(
      X = X,
      base = second_log_base,
      pseudocount = pseudocount
    )
    if (verbose) {
      THREAD_log("Input", "Applied second log transform to pseudobulk means.")
    }
  }

  if (!isTRUE(return_sparse) && is_sparse_matrix(X)) {
    X <- as.matrix(X)
  }

  validate_pseudobulk_matrix(X)

  if (verbose) {
    THREAD_log(
      "Input",
      "Generated THREAD expression matrix with ", nrow(X),
      " genes and ", ncol(X), " label column(s)."
    )
  }

  X
}


validate_scrna_arguments <- function(object,
                                     meta,
                                     subtype_col,
                                     cell_col,
                                     normalized_layer,
                                     counts_layer,
                                     already_normalized,
                                     normalization_method,
                                     scale_factor,
                                     second_log,
                                     second_log_base,
                                     pseudocount,
                                     min_cells_per_subtype,
                                     min_cells_per_gene,
                                     drop_na_subtype,
                                     return_sparse,
                                     seurat_verbose,
                                     verbose) {
  if (missing(object)) stop("'object' is required.")
  if (missing(subtype_col)) stop("'subtype_col' is required.")

  if (!is.null(meta) && !is.data.frame(meta)) {
    stop("'meta' must be a data frame when supplied.")
  }
  if (!is.null(cell_col) && (!is.character(cell_col) || length(cell_col) != 1L)) {
    stop("'cell_col' must be NULL or a single column name.")
  }
  if (!is.character(normalized_layer) || length(normalized_layer) != 1L) {
    stop("'normalized_layer' must be a single string.")
  }
  if (!is.character(counts_layer) || length(counts_layer) != 1L) {
    stop("'counts_layer' must be a single string.")
  }
  if (!is.logical(already_normalized) || length(already_normalized) != 1L) {
    stop("'already_normalized' must be TRUE or FALSE.")
  }
  if (!is.character(normalization_method) || length(normalization_method) != 1L) {
    stop("'normalization_method' must be a single string.")
  }
  if (!is.numeric(scale_factor) || length(scale_factor) != 1L || scale_factor <= 0) {
    stop("'scale_factor' must be a positive scalar.")
  }
  if (!is.logical(second_log) || length(second_log) != 1L) {
    stop("'second_log' must be TRUE or FALSE.")
  }
  if (!is.numeric(second_log_base) || length(second_log_base) != 1L ||
      !(isTRUE(all.equal(second_log_base, 2)) || isTRUE(all.equal(second_log_base, exp(1))))) {
    stop("'second_log_base' must be either 2 or exp(1).")
  }
  if (!is.numeric(pseudocount) || length(pseudocount) != 1L || pseudocount < 0) {
    stop("'pseudocount' must be a non-negative scalar.")
  }
  if (!is.numeric(min_cells_per_subtype) || length(min_cells_per_subtype) != 1L ||
      min_cells_per_subtype < 1) {
    stop("'min_cells_per_subtype' must be a positive scalar.")
  }
  if (!is.numeric(min_cells_per_gene) || length(min_cells_per_gene) != 1L ||
      min_cells_per_gene < 0) {
    stop("'min_cells_per_gene' must be a non-negative scalar.")
  }
  if (!is.logical(drop_na_subtype) || length(drop_na_subtype) != 1L) {
    stop("'drop_na_subtype' must be TRUE or FALSE.")
  }
  if (!is.logical(return_sparse) || length(return_sparse) != 1L) {
    stop("'return_sparse' must be TRUE or FALSE.")
  }
  if (!is.logical(seurat_verbose) || length(seurat_verbose) != 1L) {
    stop("'seurat_verbose' must be TRUE or FALSE.")
  }
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be TRUE or FALSE.")
  }
}


is_seurat_object <- function(x) {
  inherits(x, "Seurat")
}


get_seurat_default_assay <- function(object, assay = NULL) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for Seurat object input.")
  }
  if (!is.null(assay)) return(assay)
  Seurat::DefaultAssay(object)
}


get_seurat_assay_matrix <- function(object, assay, layer) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required to access Seurat assay data.")
  }

  out <- tryCatch(
    Seurat::GetAssayData(object = object, assay = assay, layer = layer),
    error = function(e_layer) {
      tryCatch(
        Seurat::GetAssayData(object = object, assay = assay, slot = layer),
        error = function(e_slot) {
          stop(
            "Could not extract layer/slot '", layer, "' from assay '", assay, "'. ",
            "Original errors: ", conditionMessage(e_layer), " | ", conditionMessage(e_slot)
          )
        }
      )
    }
  )
  out
}


extract_expression_from_seurat <- function(object,
                                           assay = NULL,
                                           normalized_layer = "data",
                                           counts_layer = "counts",
                                           already_normalized = FALSE,
                                           normalization_method = "LogNormalize",
                                           scale_factor = 1e4,
                                           seurat_verbose = FALSE,
                                           verbose = TRUE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for Seurat object input.")
  }

  assay <- get_seurat_default_assay(object, assay)

  if (!isTRUE(already_normalized)) {
    counts <- get_seurat_assay_matrix(object, assay = assay, layer = counts_layer)
    validate_gene_cell_matrix(counts, object_name = "Seurat counts")
    validate_nonnegative_matrix(counts, object_name = "Seurat counts")

    object <- Seurat::NormalizeData(
      object = object,
      assay = assay,
      normalization.method = normalization_method,
      scale.factor = scale_factor,
      verbose = seurat_verbose
    )
    if (verbose) {
      THREAD_log("Input", "Applied Seurat::NormalizeData() to Seurat object input.")
    }
  }

  expr <- get_seurat_assay_matrix(object, assay = assay, layer = normalized_layer)
  validate_gene_cell_matrix(expr, object_name = "Seurat normalized expression")
  validate_nonnegative_matrix(expr, object_name = "Seurat normalized expression")

  list(expr = expr, meta = object[[]], assay = assay)
}


normalize_matrix_with_seurat <- function(counts,
                                         meta,
                                         scale_factor = 1e4,
                                         normalization_method = "LogNormalize",
                                         seurat_verbose = FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for normalization.")
  }
  validate_gene_cell_matrix(counts, object_name = "counts")
  validate_nonnegative_matrix(counts, object_name = "counts")

  meta_aligned <- align_meta_to_cells(meta, colnames(counts), cell_col = NULL, verbose = FALSE)

  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = meta_aligned,
    min.cells = 0,
    min.features = 0
  )
  seu <- Seurat::NormalizeData(
    object = seu,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    verbose = seurat_verbose
  )

  assay <- Seurat::DefaultAssay(seu)
  get_seurat_assay_matrix(seu, assay = assay, layer = "data")
}

is_numeric_matrix_like <- function(x) {
  if (is.matrix(x)) {
    return(is.numeric(x))
  }

  if (is_sparse_matrix(x)) {
    if (!"x" %in% methods::slotNames(x)) {
      return(FALSE)
    }
    return(is.numeric(x@x) || is.integer(x@x))
  }

  FALSE
}

validate_gene_cell_matrix <- function(x, object_name = "matrix") {
  if (!(is.matrix(x) || is_sparse_matrix(x))) {
    stop("'", object_name, "' must be a matrix or a sparse Matrix object.")
  }

  if (!is_numeric_matrix_like(x)) {
    stop("'", object_name, "' must contain numeric values.")
  }

  if (nrow(x) == 0L || ncol(x) == 0L) {
    stop("'", object_name, "' must have at least one row and one column.")
  }

  if (is.null(rownames(x)) || anyNA(rownames(x)) || any(rownames(x) == "")) {
    stop("'", object_name, "' must have non-missing gene names as row names.")
  }

  if (anyDuplicated(rownames(x))) {
    stop("'", object_name, "' contains duplicated gene names. Please collapse or rename genes upstream.")
  }

  if (is.null(colnames(x)) || anyNA(colnames(x)) || any(colnames(x) == "")) {
    stop("'", object_name, "' must have non-missing cell IDs as column names.")
  }

  if (anyDuplicated(colnames(x))) {
    stop("'", object_name, "' contains duplicated cell IDs.")
  }

  invisible(TRUE)
}


extract_subtype_labels <- function(cell_ids,
                                   meta,
                                   subtype_col,
                                   cell_col = NULL,
                                   drop_na_subtype = TRUE,
                                   verbose = TRUE) {
  n_cells <- length(cell_ids)

  if (!is.character(subtype_col) || length(subtype_col) != 1L) {
    if (length(subtype_col) == n_cells) {
      labels <- as.character(subtype_col)
      names(labels) <- cell_ids
      if (!isTRUE(drop_na_subtype) && anyNA(labels)) {
        stop("Missing subtype labels were found and drop_na_subtype = FALSE.")
      }
      return(labels)
    }
    stop("'subtype_col' must be a single metadata column name or a vector of length ncol(expression).")
  }

  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame when 'subtype_col' is a column name.")
  }
  if (!subtype_col %in% colnames(meta)) {
    stop("Column '", subtype_col, "' was not found in metadata.")
  }

  meta_work <- align_meta_to_cells(meta, cell_ids, cell_col = cell_col, verbose = verbose)
  labels <- as.character(meta_work[[subtype_col]])
  names(labels) <- cell_ids

  if (!isTRUE(drop_na_subtype) && anyNA(labels)) {
    stop("Missing subtype labels were found and drop_na_subtype = FALSE.")
  }

  labels
}


align_meta_to_cells <- function(meta, cell_ids, cell_col = NULL, verbose = TRUE) {
  if (!is.data.frame(meta)) {
    stop("'meta' must be a data frame.")
  }

  meta_work <- meta
  if (!is.null(cell_col)) {
    if (!cell_col %in% colnames(meta_work)) {
      stop("Column '", cell_col, "' was not found in metadata.")
    }
    rownames(meta_work) <- as.character(meta_work[[cell_col]])
  }

  if (!is.null(rownames(meta_work)) && all(cell_ids %in% rownames(meta_work))) {
    meta_work <- meta_work[cell_ids, , drop = FALSE]
  } else if (nrow(meta_work) == length(cell_ids)) {
    if (verbose) {
      warning(
        "Metadata row names do not match all cell IDs. Assuming metadata rows are in the same order as expression columns.",
        call. = FALSE
      )
    }
    rownames(meta_work) <- cell_ids
  } else {
    stop(
      "Metadata cannot be aligned to cells. Provide row names matching expression column names, ",
      "a valid 'cell_col', or metadata with nrow(meta) equal to the number of cells."
    )
  }

  meta_work
}


is_sparse_matrix <- function(x) {
  inherits(x, "sparseMatrix")
}


matrix_row_sums <- function(x) {
  if (is_sparse_matrix(x)) Matrix::rowSums(x) else rowSums(x)
}


row_detected_counts <- function(x) {
  matrix_row_sums(x > 0)
}


validate_nonnegative_matrix <- function(x, object_name = "matrix") {
  if (is_sparse_matrix(x)) {
    vals <- x@x
    if (length(vals) > 0L && any(!is.finite(vals))) {
      stop("'", object_name, "' contains non-finite values.")
    }
    if (length(vals) > 0L && any(vals < 0)) {
      stop("'", object_name, "' contains negative values. Use non-negative counts or normalized expression.")
    }
  } else {
    if (any(!is.finite(x))) {
      stop("'", object_name, "' contains non-finite values.")
    }
    if (any(x < 0)) {
      stop("'", object_name, "' contains negative values. Use non-negative counts or normalized expression.")
    }
  }
  invisible(TRUE)
}


aggregate_expression_by_label <- function(x, labels) {
  labels <- factor(labels)
  if (length(labels) != ncol(x)) {
    stop("Length of 'labels' must equal ncol(x).")
  }

  group_sizes <- as.numeric(table(labels))
  group_names <- levels(labels)

  if (is_sparse_matrix(x)) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package 'Matrix' is required for sparse matrix aggregation.")
    }
    design <- Matrix::sparseMatrix(
      i = seq_along(labels),
      j = as.integer(labels),
      x = 1,
      dims = c(length(labels), length(group_names)),
      dimnames = list(colnames(x), group_names)
    )

    summed <- x %*% design
    averaged <- summed %*% Matrix::Diagonal(x = 1 / group_sizes)
    colnames(averaged) <- group_names
    rownames(averaged) <- rownames(x)
    averaged
  } else {
    summed <- rowsum(t(x), group = labels, reorder = FALSE)
    averaged <- sweep(summed, 1, group_sizes, "/")
    averaged <- t(averaged)
    colnames(averaged) <- group_names
    rownames(averaged) <- rownames(x)
    averaged
  }
}


apply_second_log_transform <- function(X, base = 2, pseudocount = 1) {
  if (pseudocount < 0) {
    stop("'pseudocount' must be non-negative.")
  }

  if (is_sparse_matrix(X)) {
    if (pseudocount > 0) {
      X <- as.matrix(X)
      return(apply_second_log_transform(X, base = base, pseudocount = pseudocount))
    }
    if (length(X@x) > 0L && any(X@x <= 0)) {
      stop("The second log transform with pseudocount = 0 requires positive non-zero entries.")
    }
    out <- X
    if (isTRUE(all.equal(base, 2))) {
      out@x <- log2(out@x)
    } else if (isTRUE(all.equal(base, exp(1)))) {
      out@x <- log(out@x)
    } else {
      stop("Unsupported log base.")
    }
    out
  } else {
    if (any(X + pseudocount <= 0)) {
      stop("The second log transform requires X + pseudocount > 0.")
    }
    if (isTRUE(all.equal(base, 2))) {
      log2(X + pseudocount)
    } else if (isTRUE(all.equal(base, exp(1)))) {
      log(X + pseudocount)
    } else {
      stop("Unsupported log base.")
    }
  }
}


validate_pseudobulk_matrix <- function(X) {
  if (!(is.matrix(X) || is_sparse_matrix(X))) {
    stop("The output expression matrix must be a matrix or sparse Matrix object.")
  }

  if (!is_numeric_matrix_like(X)) {
    stop("The output expression matrix must contain numeric values.")
  }

  if (nrow(X) == 0L || ncol(X) == 0L) {
    stop("The output expression matrix must have at least one gene and one label column.")
  }

  if (is_sparse_matrix(X)) {
    if (length(X@x) > 0L && any(!is.finite(X@x))) {
      stop("The output expression matrix contains non-finite values.")
    }
  } else {
    if (any(!is.finite(X))) {
      stop("The output expression matrix contains non-finite values.")
    }
  }

  if (is.null(rownames(X)) || anyNA(rownames(X)) || any(rownames(X) == "")) {
    stop("The output expression matrix must have gene names as row names.")
  }

  if (is.null(colnames(X)) || anyNA(colnames(X)) || any(colnames(X) == "")) {
    stop("The output expression matrix must have label names as column names.")
  }

  invisible(TRUE)
}