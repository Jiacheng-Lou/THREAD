#' Read a gene2vec embedding file
#'
#' Read an external gene embedding file for gene2vec-based initialization in DPM.
#' The file should contain one gene identifier column followed by numeric
#' embedding dimensions. A Word2Vec-style metadata header, such as "25442 200",
#' is automatically detected and skipped.
#'
#' @param path Character string. Path to the gene2vec embedding file.
#' @param n_max Integer or Inf. Maximum number of rows to read after skipping a
#'   possible metadata header. Use Inf to read the full file.
#' @param gene_col Integer. Column index containing gene identifiers.
#' @param skip Integer or NULL. Number of rows to skip. If NULL, the function
#'   automatically detects a Word2Vec-style metadata header.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A numeric matrix with genes as row names and embedding dimensions as columns.
#'
#' @export
read_gene2vec <- function(path,
                          n_max = Inf,
                          gene_col = 1L,
                          skip = NULL,
                          verbose = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("'path' must be a single non-missing character string.")
  }

  if (!file.exists(path)) {
    stop("gene2vec file does not exist: ", path)
  }

  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) || n_max <= 0) {
    stop("'n_max' must be a positive numeric value or Inf.")
  }

  if (!is.numeric(gene_col) || length(gene_col) != 1L || is.na(gene_col)) {
    stop("'gene_col' must be a single column index.")
  }

  gene_col <- as.integer(gene_col)

  if (verbose) {
    THREAD_log("Input", "Reading gene2vec file: ", path)
  }

  first_line <- readLines(path, n = 1L, warn = FALSE)
  if (length(first_line) == 0L) {
    stop("gene2vec file is empty: ", path)
  }

  if (is.null(skip)) {
    tokens <- strsplit(trimws(first_line), "\\s+")[[1L]]
    is_word2vec_header <- length(tokens) == 2L &&
      !is.na(suppressWarnings(as.numeric(tokens[1L]))) &&
      !is.na(suppressWarnings(as.numeric(tokens[2L])))

    skip_rows <- if (is_word2vec_header) 1L else 0L

    if (verbose && is_word2vec_header) {
      THREAD_log("Input", "Detected Word2Vec metadata header. Skipping first row.")
    }
  } else {
    if (!is.numeric(skip) || length(skip) != 1L || is.na(skip) || skip < 0) {
      stop("'skip' must be NULL or a non-negative integer.")
    }
    skip_rows <- as.integer(skip)
  }

  read_args <- list(
    file = path,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    skip = skip_rows,
    strip.white = TRUE,
    quote = "",
    comment.char = ""
  )

  if (!is.infinite(n_max)) {
    read_args$nrows <- as.integer(n_max)
  }

  dt <- do.call(utils::read.table, read_args)

  if (nrow(dt) == 0L) {
    stop("No embedding rows were read from gene2vec file.")
  }

  if (ncol(dt) < 2L) {
    stop("gene2vec file must contain one gene column and at least one embedding dimension.")
  }

  if (gene_col < 1L || gene_col > ncol(dt)) {
    stop("'gene_col' is outside the number of columns in the gene2vec file.")
  }

  embedding_cols <- setdiff(seq_len(ncol(dt)), gene_col)

  first_embedding_col <- suppressWarnings(as.numeric(dt[[embedding_cols[1L]]]))
  if (length(first_embedding_col) > 0L && is.na(first_embedding_col[1L])) {
    if (verbose) {
      THREAD_log("Input", "Detected a non-numeric first data row. Treating it as a header row and removing it.")
    }
    dt <- dt[-1L, , drop = FALSE]
  }

  genes <- as.character(dt[[gene_col]])
  emb <- as.matrix(dt[, embedding_cols, drop = FALSE])
  suppressWarnings(storage.mode(emb) <- "numeric")

  finite_rows <- apply(emb, 1L, function(x) all(is.finite(x)))
  keep <- !is.na(genes) & genes != "" & !duplicated(genes) & finite_rows
  dropped <- sum(!keep)

  genes <- genes[keep]
  emb <- emb[keep, , drop = FALSE]

  if (length(genes) == 0L) {
    stop("No valid gene embeddings remained after filtering.")
  }

  rownames(emb) <- genes
  colnames(emb) <- paste0("dim", seq_len(ncol(emb)))

  if (verbose) {
    THREAD_log("Input", "Loaded gene2vec matrix: ", nrow(emb), " genes x ", ncol(emb), " dimensions.")
    if (dropped > 0L) {
      THREAD_log("Input", "Dropped ", dropped, " invalid or duplicated row(s).")
    }
  }


  emb
}

validate_gene2vec_matrix <- function(gene2vec) {
  if (!is.matrix(gene2vec)) {
    stop("'gene2vec' must be a numeric matrix with genes as row names.")
  }

  if (!is.numeric(gene2vec)) {
    stop("'gene2vec' must be numeric.")
  }

  if (nrow(gene2vec) == 0L || ncol(gene2vec) == 0L) {
    stop("'gene2vec' must have at least one row and one column.")
  }

  if (is.null(rownames(gene2vec)) || anyNA(rownames(gene2vec)) || any(rownames(gene2vec) == "")) {
    stop("'gene2vec' must have non-missing gene names as row names.")
  }

  if (anyDuplicated(rownames(gene2vec))) {
    stop("'gene2vec' contains duplicated gene names.")
  }

  if (any(!is.finite(gene2vec))) {
    stop("'gene2vec' contains non-finite values.")
  }

  invisible(TRUE)
}


#' Read a gene annotation table
#'
#' Read a gene annotation table from a local text file. This is a lightweight
#' reader only; column standardization is performed by standardize_gene_annotation().
#'
#' @param path Character string. Path to a gene annotation file.
#' @param ... Additional arguments passed to data.table::fread().
#'
#' @return A data.frame containing the raw gene annotation table.
#'
#' @export
read_gene_annotation <- function(path, ...) {
  if (!is.character(path) || length(path) != 1L || is.na(path)) {
    stop("'path' must be a single non-missing character string.")
  }

  if (!file.exists(path)) {
    stop("Gene annotation file does not exist: ", path)
  }

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  ann <- data.table::fread(
    file = path,
    data.table = FALSE,
    showProgress = FALSE,
    ...
  )

  if (nrow(ann) == 0L) {
    stop("Gene annotation file contains no rows.")
  }

  ann
}