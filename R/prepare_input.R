# prepare_input.R
# -----------------------------------------------------------------------------
# Data assembly for the DPM model.
#
#   prepare_input()      Align X and GWAS response, create z_init and priors,
#                        and return a thread_data object consumed by run_thread().
#   init_clusters()      Build starting gene labels from gene2vec when supplied;
#                        otherwise fall back to k-means on cbind(y, X).
#   recommend_priors()   Empirical-Bayes recommendation for alpha0/beta0/r/delta
#                        using a weighted ridge approximation.
#
# The response y is kept on the raw LD-corrected scale. This is intentional:
# y_train, d0, d1 and d2 must remain on the same scale as the LD-derived
# heteroscedastic likelihood. Expression X can optionally be transformed or
# scaled through explicit arguments.
# -----------------------------------------------------------------------------

# raw or normalized gene × cell matrix + metadata labels
#         ↓ preprocess_scrna()
# log-normalized pseudobulk gene × label matrix X
#         ↓ prepare_input()
# thread_data


#' Detect a gene identifier column
#'
#' @param x A data frame-like object.
#' @param gene_col Optional user-specified gene column.
#' @return Column name.
#' @noRd
get_gene_column <- function(x, gene_col = NULL) {
  if (!is.null(gene_col)) {
    if (!gene_col %in% names(x)) stop("Column '", gene_col, "' was not found.")
    return(gene_col)
  }

  candidates <- c("gene_name", "gene", "Gene", "GENE", "symbol", "SYMBOL")
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0L) {
    stop("Could not find a gene column. Provide 'gene_col'.")
  }
  hit[1]
}


#' Detect a response column
#'
#' @param response A response data frame.
#' @param response_col Optional user-specified response column.
#' @return Column name.
#' @noRd
get_response_column <- function(response, response_col = "final_y") {
  if (!is.null(response_col) && response_col %in% names(response)) return(response_col)

  fallback <- c("final_y", "y", "Y", "response")
  hit <- fallback[fallback %in% names(response)]
  if (length(hit) == 0L) {
    stop("Could not find the response column. Expected 'final_y' or provide 'response_col'.")
  }
  if (!identical(hit[1], response_col)) {
    warning("Using response column '", hit[1], "'. For package-standard inputs, prefer 'final_y'.")
  }
  hit[1]
}


#' Convert an object to a numeric matrix and validate row names
#'
#' @param X Matrix-like object.
#' @param arg_name Object name used in error messages.
#' @return Numeric matrix.
#' @noRd
as_numeric_matrix_checked <- function(X, arg_name = "X") {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.matrix(X)) stop("'", arg_name, "' must be a matrix or data frame.")
  storage.mode(X) <- "double"

  if (is.null(rownames(X))) {
    stop("'", arg_name, "' must have row names containing gene identifiers.")
  }
  if (anyNA(rownames(X)) || any(rownames(X) == "")) {
    stop("'", arg_name, "' contains missing or empty gene row names.")
  }
  if (anyDuplicated(rownames(X))) {
    dup <- unique(rownames(X)[duplicated(rownames(X))])
    stop("'", arg_name, "' contains duplicated gene row names, for example: ", dup[1], ".")
  }
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Subtype", seq_len(ncol(X)))
  }
  X
}


#' Apply an explicit expression transform
#'
#' @param X Numeric matrix.
#' @param x_transform One of none/log2p1/log1p.
#' @return Numeric matrix.
#' @noRd
transform_expression_matrix <- function(X, x_transform = c("none", "log2p1", "log1p")) {
  x_transform <- match.arg(x_transform)

  if (x_transform == "none") return(X)

  if (any(X < 0, na.rm = TRUE)) {
    stop("Expression transformation requires non-negative values in 'X'.")
  }

  if (x_transform == "log2p1") {
    return(log2(X + 1))
  }
  if (x_transform == "log1p") {
    return(log1p(X))
  }

  X
}


#' Center and scale a numeric matrix safely
#'
#' @param X Numeric matrix.
#' @param center Logical.
#' @param scale Logical.
#' @return Numeric matrix.
#' @noRd
center_scale_matrix <- function(X, center = TRUE, scale = TRUE) {
  X <- as.matrix(X)

  if (center) {
    means <- colMeans(X, na.rm = TRUE)
    X <- sweep(X, 2L, means, FUN = "-")
  }

  if (scale) {
    sds <- apply(X, 2L, stats::sd, na.rm = TRUE)
    sds[!is.finite(sds) | sds == 0] <- 1
    X <- sweep(X, 2L, sds, FUN = "/")
  }

  X[!is.finite(X)] <- 0
  X
}


#' Build a compact aligned response data frame
#'
#' @param genes Character vector.
#' @param y,d0,d1,d2 Numeric vectors.
#' @return Data frame.
#' @noRd
make_response_frame <- function(genes, y, d0, d1, d2) {
  data.frame(
    gene_name = genes,
    final_y = as.numeric(y),
    d0 = as.numeric(d0),
    d1 = as.numeric(d1),
    d2 = as.numeric(d2),
    stringsAsFactors = FALSE
  )
}


#' Align expression and response objects
#'
#' @param X Gene-by-subtype matrix.
#' @param response Data frame with gene, final_y, d0, d1 and d2 columns.
#' @param gene_col Optional gene column name in response.
#' @param response_col Response column name. Default final_y.
#' @param min_genes Minimum number of overlapping genes.
#' @param drop_nonfinite Drop rows with non-finite values.
#' @param verbose Logical.
#' @return List of aligned X, y, d0, d1, d2 and gene names.
#' @noRd
align_x_response <- function(X, response,
                             gene_col = NULL,
                             response_col = "final_y",
                             min_genes = 10,
                             drop_nonfinite = TRUE,
                             verbose = TRUE) {
  X <- as_numeric_matrix_checked(X, "X")
  response <- as.data.frame(
    response,
    stringsAsFactors = FALSE
  )

  if ("status" %in% names(response)) {
    n_response_before <- nrow(response)

    response <- response[
      !is.na(response$status) &
        response$status == "ok", ,
      drop = FALSE
    ]

    if (verbose) {
      DPM_log(
        "Input",
        sprintf(
          "Using %d response row(s) with status == 'ok'; excluded %d non-ok row(s).",
          nrow(response), n_response_before - nrow(response)
        )
      )
    }

    if (nrow(response) == 0L) {
      stop(
        "No response rows with status == 'ok' remain."
      )
    }
  }

  gene_col <- get_gene_column(response, gene_col)


  response_col <- get_response_column(response, response_col)

  required <- c(response_col, "d0", "d1", "d2")
  missing <- setdiff(required, names(response))
  if (length(missing) > 0L) {
    stop("'response' is missing required column(s): ", paste(missing, collapse = ", "), ".")
  }

  response_genes <- as.character(response[[gene_col]])
  if (anyNA(response_genes) || any(response_genes == "")) {
    stop("The response gene column contains missing or empty gene names.")
  }
  if (anyDuplicated(response_genes)) {
    dup <- unique(response_genes[duplicated(response_genes)])
    stop("'response' contains duplicated genes, for example: ", dup[1], ".")
  }

  common_genes <- intersect(rownames(X), response_genes)
  if (length(common_genes) < min_genes) {
    stop(
      "Only ", length(common_genes), " genes overlap between X and response; ",
      "at least ", min_genes, " are required."
    )
  }

  X_aligned <- X[common_genes, , drop = FALSE]
  idx <- match(common_genes, response_genes)

  y <- as.numeric(response[[response_col]][idx])
  d0 <- as.numeric(response[["d0"]][idx])
  d1 <- as.numeric(response[["d1"]][idx])
  d2 <- as.numeric(response[["d2"]][idx])

  keep <- rep(TRUE, length(common_genes))
  if (drop_nonfinite) {
    keep <- is.finite(y) & is.finite(d0) & is.finite(d1) & is.finite(d2) &
      apply(X_aligned, 1L, function(v) all(is.finite(v)))

    dropped <- sum(!keep)
    if (dropped > 0L && verbose) {
      warning("Dropping ", dropped, " genes with non-finite X/y/d values.")
    }
  }

  common_genes <- common_genes[keep]
  X_aligned <- X_aligned[keep, , drop = FALSE]
  y <- y[keep]
  d0 <- d0[keep]
  d1 <- d1[keep]
  d2 <- d2[keep]

  if (length(common_genes) < min_genes) {
    stop(
      "After removing invalid rows, only ", length(common_genes),
      " genes remain; at least ", min_genes, " are required."
    )
  }

  names(y) <- common_genes
  names(d0) <- common_genes
  names(d1) <- common_genes
  names(d2) <- common_genes

  if (verbose) {
    DPM_log(
      "Input",
      sprintf(
        "Aligned %d genes across X (%d total) and response (%d total).",
        length(common_genes), nrow(X), nrow(response)
      )
    )
  }

  list(
    x_train = X_aligned,
    y_train = y,
    d0 = d0,
    d1 = d1,
    d2 = d2,
    gene_names = common_genes,
    response_frame = make_response_frame(common_genes, y, d0, d1, d2)
  )
}


#' Align user-supplied initial labels to model genes
#'
#' @param z_init Vector or data frame of starting labels.
#' @param genes Model gene order.
#' @param gene_col Optional gene column when z_init is a data frame.
#' @return Integer vector with contiguous labels.
#' @noRd
align_z_init <- function(z_init, genes, gene_col = NULL) {
  if (is.null(z_init)) return(NULL)

  if (is.data.frame(z_init)) {
    zdf <- as.data.frame(z_init, stringsAsFactors = FALSE)
    gcol <- get_gene_column(zdf, gene_col)
    label_candidates <- c("cluster", "clusters", "z", "z_init", "label", "module")
    lcol <- label_candidates[label_candidates %in% names(zdf)]
    if (length(lcol) == 0L) {
      stop("'z_init' data frame must contain a label column such as 'cluster' or 'z'.")
    }
    z_genes <- as.character(zdf[[gcol]])
    if (anyDuplicated(z_genes)) {
      dup <- unique(z_genes[duplicated(z_genes)])
      stop("'z_init' contains duplicated genes, for example: ", dup[1], ".")
    }
    idx <- match(genes, z_genes)
    if (anyNA(idx)) {
      missing <- genes[is.na(idx)]
      stop("'z_init' is missing labels for ", length(missing), " gene(s), for example: ", missing[1], ".")
    }
    z <- zdf[[lcol[1]]][idx]
  } else {
    z <- z_init
    if (!is.null(names(z))) {
      idx <- match(genes, names(z))
      if (anyNA(idx)) {
        missing <- genes[is.na(idx)]
        stop("Named 'z_init' is missing labels for ", length(missing), " gene(s), for example: ", missing[1], ".")
      }
      z <- z[idx]
    } else if (length(z) != length(genes)) {
      stop("Unnamed 'z_init' must have the same length as the aligned gene set.")
    }
  }

  if (anyNA(z)) stop("'z_init' contains missing labels after alignment.")
  z <- match(z, sort(unique(z)))
  z <- as.integer(z)
  names(z) <- genes
  z
}


#' Load or validate a gene2vec matrix
#'
#' @param gene2vec Matrix/data frame or a path to a whitespace-delimited file.
#' @return Numeric matrix with genes as row names.
#' @noRd
load_gene2vec_matrix <- function(gene2vec) {
  if (is.null(gene2vec)) {
    return(NULL)
  }

  if (is.character(gene2vec) && length(gene2vec) == 1L) {
    return(read_gene2vec(gene2vec, verbose = FALSE))
  }

  if (is.character(gene2vec) && length(gene2vec) == 1L) {
    if (!file.exists(gene2vec)) stop("gene2vec file does not exist: ", gene2vec)
    g <- utils::read.table(gene2vec, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(g) < 2L) stop("gene2vec file must contain one gene column and at least one embedding column.")
    rownames(g) <- as.character(g[[1]])
    g <- g[, -1L, drop = FALSE]
  } else {
    g <- gene2vec
  }

  if (is.data.frame(g)) {
    if (is.null(rownames(g)) || all(rownames(g) == as.character(seq_len(nrow(g))))) {
      first_col_is_gene <- ncol(g) >= 2L && !is.numeric(g[[1]])
      if (first_col_is_gene) {
        rownames(g) <- as.character(g[[1]])
        g <- g[, -1L, drop = FALSE]
      }
    }
    g <- as.matrix(g)
  }

  if (!is.matrix(g)) stop("'gene2vec' must be a matrix, data frame or file path.")
  if (is.null(rownames(g))) stop("'gene2vec' must have row names containing gene identifiers.")
  if (anyDuplicated(rownames(g))) {
    dup <- unique(rownames(g)[duplicated(rownames(g))])
    stop("'gene2vec' contains duplicated genes, for example: ", dup[1], ".")
  }

  storage.mode(g) <- "double"
  bad_rows <- apply(g, 1L, function(v) !all(is.finite(v)))
  if (any(bad_rows)) {
    warning("Dropping ", sum(bad_rows), " genes with non-finite gene2vec values.")
    g <- g[!bad_rows, , drop = FALSE]
  }

  g
}


#' Majority vote with deterministic or random tie handling
#'
#' @param labels Integer labels.
#' @param tie_method Tie handling.
#' @return One label.
#' @noRd
vote_label <- function(labels, tie_method = c("first", "random")) {
  tie_method <- match.arg(tie_method)
  labels <- labels[!is.na(labels)]
  if (length(labels) == 0L) return(NA_integer_)

  tab <- table(labels)
  winners <- as.integer(names(tab)[tab == max(tab)])
  winners <- sort(winners)

  if (length(winners) == 1L || tie_method == "first") {
    winners[1]
  } else {
    sample(winners, 1L)
  }
}


#' Generate initial gene-cluster labels
#'
#' @description
#' Uses gene2vec embeddings when supplied. Genes present in gene2vec are first
#' clustered by k-means in embedding space. Genes missing from gene2vec are
#' assigned by KNN voting in the expression-profile space. If gene2vec is not
#' supplied, the function falls back to k-means on cbind(y, X).
#'
#' @param X Gene-by-subtype matrix.
#' @param response Response data frame with gene, final_y, d0, d1 and d2.
#' @param gene2vec Optional matrix/data frame/path with genes as rows.
#' @param K_init Initial number of clusters. Default 10.
#' @param k_neighbors Number of neighbours for KNN voting. Default 100.
#' @param gene_col Optional gene column name in response.
#' @param response_col Response column name. Default final_y.
#' @param seed Optional random seed.
#' @param nstart Number of k-means starts. Default 25.
#' @param iter.max Maximum k-means iterations. Default 1000.
#' @param tie_method Tie handling for KNN voting, "first" or "random".
#' @param min_genes Minimum overlap required between X and response.
#' @param verbose Logical.
#'
#' @return Named integer vector of length N, aligned to the intersected genes.
#'
#' @importFrom stats kmeans
#' @export
init_clusters <- function(X, response,
                          gene2vec = NULL,
                          K_init = 10,
                          k_neighbors = 100,
                          gene_col = NULL,
                          response_col = "final_y",
                          seed = NULL,
                          nstart = 25,
                          iter.max = 1000,
                          tie_method = c("first", "random"),
                          min_genes = 10,
                          verbose = TRUE) {
  tie_method <- match.arg(tie_method)

  aligned <- align_x_response(
    X = X,
    response = response,
    gene_col = gene_col,
    response_col = response_col,
    min_genes = min_genes,
    drop_nonfinite = TRUE,
    verbose = FALSE
  )

  X_aligned <- aligned$x_train
  y <- aligned$y_train
  genes <- aligned$gene_names
  N <- length(genes)

  if (!is.numeric(K_init) || length(K_init) != 1L || K_init < 1) {
    stop("'K_init' must be a positive scalar.")
  }
  K_eff <- min(as.integer(K_init), N)
  if (K_eff < K_init && verbose) {
    warning("Reducing K_init from ", K_init, " to ", K_eff, " because there are only ", N, " genes.")
  }

  if (!is.null(seed)) set.seed(seed)

  g2v <- load_gene2vec_matrix(gene2vec)

  if (!is.null(g2v)) {
    genes_in_vec <- genes[genes %in% rownames(g2v)]
    if (length(genes_in_vec) >= K_eff) {
      if (verbose) {
        DPM_log(
          "Input",
          sprintf(
            "Initialisation: gene2vec k-means for %d/%d genes; KNN voting for %d missing genes.",
            length(genes_in_vec), N, N - length(genes_in_vec)
          )
        )
      }

      embedding <- g2v[genes_in_vec, , drop = FALSE]
      embedding <- center_scale_matrix(embedding, center = TRUE, scale = TRUE)

      km <- stats::kmeans(
        embedding,
        centers = K_eff,
        nstart = nstart,
        iter.max = iter.max
      )

      z <- rep(NA_integer_, N)
      names(z) <- genes
      z[genes_in_vec] <- as.integer(km$cluster)

      missing_genes <- genes[is.na(z)]
      if (length(missing_genes) > 0L) {
        if (!requireNamespace("FNN", quietly = TRUE)) {
          stop(
            "Package 'FNN' is required for KNN assignment of genes missing from gene2vec. ",
            "Install FNN or call init_clusters() without gene2vec to use k-means fallback."
          )
        }

        k_eff_nn <- min(as.integer(k_neighbors), length(genes_in_vec))
        if (k_eff_nn < 1L) stop("No labelled genes are available for KNN voting.")

        X_scaled <- center_scale_matrix(X_aligned, center = TRUE, scale = TRUE)
        labelled_matrix <- X_scaled[genes_in_vec, , drop = FALSE]
        missing_matrix <- X_scaled[missing_genes, , drop = FALSE]

        nn <- FNN::get.knnx(labelled_matrix, missing_matrix, k = k_eff_nn)
        labelled_z <- z[genes_in_vec]

        for (i in seq_along(missing_genes)) {
          neighbour_labels <- labelled_z[nn$nn.index[i, ]]
          z[missing_genes[i]] <- vote_label(neighbour_labels, tie_method = tie_method)
        }
      }

      if (anyNA(z)) stop("Internal error: some initial labels remain missing.")
      z <- match(z, sort(unique(z)))
      z <- as.integer(z)
      names(z) <- genes
      return(z)
    }

    warning(
      "Only ", length(genes_in_vec), " aligned genes are present in gene2vec, ",
      "which is fewer than K_init = ", K_eff, ". Falling back to k-means on cbind(y, X)."
    )
  }

  if (verbose) {
    DPM_log("Input", sprintf("Initialisation: k-means fallback on cbind(y, X), K = %d.", K_eff))
  }

  features <- cbind(y_train = as.numeric(y), X_aligned)
  features <- center_scale_matrix(features, center = TRUE, scale = TRUE)

  km <- stats::kmeans(
    features,
    centers = K_eff,
    nstart = nstart,
    iter.max = iter.max
  )

  z <- as.integer(km$cluster)
  z <- match(z, sort(unique(z)))
  z <- as.integer(z)
  names(z) <- genes
  z
}


#' Validate or complete prior hyperparameters
#'
#' @param priors Prior list.
#' @return Prior list.
#' @noRd
validate_priors <- function(priors) {
  if (is.null(priors)) return(NULL)
  needed <- c("alpha0", "beta0", "r", "delta")
  missing <- setdiff(needed, names(priors))
  if (length(missing) > 0L) {
    stop("'priors' is missing required field(s): ", paste(missing, collapse = ", "), ".")
  }

  priors$alpha0 <- as.numeric(priors$alpha0)
  priors$beta0 <- as.numeric(priors$beta0)
  priors$r <- as.numeric(priors$r)
  priors$delta <- as.numeric(priors$delta)

  vals <- unlist(priors[needed], use.names = TRUE)
  if (any(!is.finite(vals)) || any(vals <= 0)) {
    stop("Prior hyperparameters alpha0, beta0, r and delta must be finite and strictly positive.")
  }

  sigma2_prior_mode <- priors$beta0 / (priors$alpha0 + 1)
  tau2_prior_mean <- 2 / (priors$r / priors$delta)
  priors$prior_sd <- sqrt(sigma2_prior_mode * tau2_prior_mean)

  priors
}


#' Fit a weighted ridge regression and return coefficient variance
#'
#' @param X Numeric matrix.
#' @param y Numeric vector.
#' @param w Numeric weights.
#' @param lambda Ridge penalty.
#' @return Numeric scalar.
#' @noRd
ridge_beta_variance <- function(X, y, w, lambda) {
  p <- ncol(X)
  x_w <- X * as.numeric(w)
  y_w <- y * as.numeric(w)

  beta <- tryCatch({
    as.vector(solve(crossprod(x_w) + lambda * diag(p), crossprod(x_w, y_w)))
  }, error = function(e) {
    rep(0, p)
  })

  v <- stats::var(beta)
  if (!is.finite(v) || v <= 0) NA_real_ else v
}


#' Recommend empirical-Bayes prior hyperparameters
#'
#' @description
#' Estimates a rough coefficient scale through a weighted ridge approximation.
#' The response is kept on its raw LD-corrected scale; no log transform is
#' applied to y.
#'
#' @param X Gene-by-subtype matrix.
#' @param response Response data frame with gene, final_y, d0, d1 and d2.
#' @param z_init Optional initial labels, vector or data frame.
#' @param gene_col Optional gene column name.
#' @param response_col Response column name. Default final_y.
#' @param alpha0 Default inverse-gamma shape for sigma2. Default 2.
#' @param r Default gamma shape for lambda2. Default 1.
#' @param delta Default gamma rate for lambda2. Default 1.
#' @param min_beta0 Lower bound for beta0. Default 1e-16.
#' @param min_var Lower bound for approximate variance components. Default 1e-12.
#' @param ridge_scale Multiplicative scale for the ridge penalty. Default 1.
#' @param min_cluster_genes Minimum genes required to use a cluster-specific
#'   ridge fit. Defaults to max(p + 1, 10).
#' @param min_genes Minimum overlap required between X and response.
#' @param verbose Logical.
#'
#' @return List with alpha0, beta0, r, delta and diagnostic fields.
#'
#' @export
recommend_priors <- function(X, response,
                             z_init = NULL,
                             gene_col = NULL,
                             response_col = "final_y",
                             alpha0 = 2.0,
                             r = 1.0,
                             delta = 1.0,
                             min_beta0 = 1e-16,
                             min_var = 1e-12,
                             ridge_scale = 1.0,
                             min_cluster_genes = NULL,
                             min_genes = 10,
                             verbose = TRUE) {

  aligned <- align_x_response(
    X = X,
    response = response,
    gene_col = gene_col,
    response_col = response_col,
    min_genes = min_genes,
    drop_nonfinite = TRUE,
    verbose = FALSE
  )

  X_aligned <- aligned$x_train
  y <- aligned$y_train
  d0 <- aligned$d0
  genes <- aligned$gene_names

  p <- ncol(X_aligned)
  N <- nrow(X_aligned)
  if (is.null(min_cluster_genes)) min_cluster_genes <- max(p + 1L, 10L)

  weights <- 1 / sqrt(pmax(d0, min_var))
  weights[!is.finite(weights)] <- 1 / sqrt(min_var)

  base_lambda <- ridge_scale * sum(X_aligned^2) / max(N, 1L)
  if (!is.finite(base_lambda) || base_lambda <= 0) base_lambda <- 1

  z <- align_z_init(z_init, genes, gene_col = gene_col)
  beta_variances <- numeric(0)

  if (!is.null(z) && length(unique(z)) > 1L) {
    clusters <- sort(unique(z))
    if (verbose) {
      DPM_log("Input", sprintf("Prior recommendation: cluster-wise weighted ridge over K = %d initial clusters.", length(clusters)))
    }

    for (k in clusters) {
      idx <- which(z == k)
      if (length(idx) >= min_cluster_genes) {
        lambda_k <- base_lambda * length(idx) / N
        v_k <- ridge_beta_variance(
          X = X_aligned[idx, , drop = FALSE],
          y = y[idx],
          w = weights[idx],
          lambda = lambda_k
        )
        if (is.finite(v_k) && v_k > 0) beta_variances <- c(beta_variances, v_k)
      }
    }
  }

  if (length(beta_variances) > 0L) {
    emp_beta_var <- mean(beta_variances)
    mode <- "cluster-wise"
  } else {
    if (verbose) DPM_log("Input", "Prior recommendation: global weighted ridge fit.")
    emp_beta_var <- ridge_beta_variance(X_aligned, y, weights, base_lambda)
    mode <- "global"
  }

  if (!is.finite(emp_beta_var) || emp_beta_var <= 0) {
    emp_beta_var <- stats::var(y, na.rm = TRUE) / max(base_lambda, min_var)
  }
  if (!is.finite(emp_beta_var) || emp_beta_var <= 0) emp_beta_var <- 1e-12

  beta0 <- max(emp_beta_var / 2, min_beta0)

  priors <- list(
    alpha0 = as.numeric(alpha0),
    beta0 = as.numeric(beta0),
    r = as.numeric(r),
    delta = as.numeric(delta),
    beta_variance = as.numeric(emp_beta_var),
    base_lambda = as.numeric(base_lambda),
    mode = mode
  )
  priors <- validate_priors(priors)

  if (verbose) {
    DPM_log(
      "Input",
      sprintf(
        "Recommended priors: alpha0 = %.3g, beta0 = %.6e, r = %.3g, delta = %.3g, prior_sd = %.6e.",
        priors$alpha0, priors$beta0, priors$r, priors$delta, priors$prior_sd
      )
    )
  }

  priors
}


#' Prepare DPM model input
#'
#' @description
#' Aligns a gene-by-subtype expression matrix with a gene-level GWAS response,
#' creates or validates initial cluster labels, recommends prior hyperparameters
#' when they are not supplied, and returns a \code{thread_data} object consumed by
#' \code{run_thread()}.
#'
#' The response is kept on the raw LD-corrected scale. No log transform or
#' truncation is applied to y, because y, d0, d1 and d2 must remain on the same
#' scale as the heteroscedastic likelihood.
#'
#' @param X Gene-by-subtype expression matrix. Rows must be genes.
#' @param response Data frame with gene, final_y, d0, d1 and d2 columns.
#' @param z_init Optional initial labels. Can be a vector aligned to genes, a
#'   named vector, or a data frame with gene and cluster columns.
#' @param priors Optional list with alpha0, beta0, r and delta. When NULL,
#'   \code{recommend_priors()} is called.
#' @param gene2vec Optional gene2vec matrix/data frame/path used only when
#'   \code{z_init} is NULL.
#' @param K_init Initial number of clusters when z_init is not supplied.
#' @param k_neighbors Number of neighbours for gene2vec-missing genes.
#' @param gene_col Optional gene column name in response and z_init data frames.
#' @param response_col Response column name. Default final_y.
#' @param x_transform Expression transform: none, log2p1 or log1p. Default none.
#' @param center_x Center expression columns before model fitting. Default FALSE.
#' @param scale_x Scale expression columns before model fitting. Default FALSE.
#' @param min_genes Minimum number of overlapping valid genes. Default 50.
#' @param seed Optional seed used by init_clusters.
#' @param verbose Logical.
#' @param ... Reserved for future extensions.
#'
#' @return A list of class \code{thread_data} with x_train, y_train, d0, d1, d2,
#'   z_init, gene_names and priors.
#'
#' @export
prepare_input <- function(X, response,
                          z_init = NULL,
                          priors = NULL,
                          gene2vec = NULL,
                          K_init = 10,
                          k_neighbors = 100,
                          gene_col = NULL,
                          response_col = "final_y",
                          x_transform = c("none", "log2p1", "log1p"),
                          center_x = FALSE,
                          scale_x = FALSE,
                          min_genes = 50,
                          seed = NULL,
                          verbose = TRUE,
                          ...) {
  x_transform <- match.arg(x_transform)

  X <- as_numeric_matrix_checked(X, "X")
  X <- transform_expression_matrix(X, x_transform = x_transform)

  aligned <- align_x_response(
    X = X,
    response = response,
    gene_col = gene_col,
    response_col = response_col,
    min_genes = min_genes,
    drop_nonfinite = TRUE,
    verbose = verbose
  )

  x_train <- aligned$x_train
  if (center_x || scale_x) {
    x_train <- center_scale_matrix(x_train, center = center_x, scale = scale_x)
  }

  gene_names <- aligned$gene_names
  y_train <- aligned$y_train
  d0 <- aligned$d0
  d1 <- aligned$d1
  d2 <- aligned$d2

  response_aligned <- make_response_frame(gene_names, y_train, d0, d1, d2)

  if (is.null(z_init)) {
    z <- init_clusters(
      X = x_train,
      response = response_aligned,
      gene2vec = gene2vec,
      K_init = K_init,
      k_neighbors = k_neighbors,
      gene_col = "gene_name",
      response_col = "final_y",
      seed = seed,
      min_genes = min_genes,
      verbose = verbose
    )
  } else {
    z <- align_z_init(z_init, gene_names, gene_col = gene_col)
    if (verbose) {
      DPM_log("Input", sprintf("Using user-supplied z_init with K = %d.", length(unique(z))))
    }
  }

  names(z) <- gene_names

  if (is.null(priors)) {
    priors <- recommend_priors(
      X = x_train,
      response = response_aligned,
      z_init = z,
      gene_col = "gene_name",
      response_col = "final_y",
      min_genes = min_genes,
      verbose = verbose
    )
  } else {
    priors <- validate_priors(priors)
    if (verbose) {
      DPM_log("Input", "Using user-supplied priors.")
    }
  }

  names(y_train) <- gene_names
  names(d0) <- gene_names
  names(d1) <- gene_names
  names(d2) <- gene_names
  rownames(x_train) <- gene_names

  out <- list(
    x_train = x_train,
    y_train = y_train,
    d0 = d0,
    d1 = d1,
    d2 = d2,
    z_init = as.integer(z),
    gene_names = gene_names,
    priors = priors,
    settings = list(
      x_transform = x_transform,
      center_x = center_x,
      scale_x = scale_x,
      K_init = K_init,
      k_neighbors = k_neighbors,
      response_col = response_col
    )
  )

  names(out$z_init) <- gene_names

  class(out) <- "thread_data"

  if (verbose) {
    DPM_log(
      "Input",
      sprintf(
        "Prepared thread_data: %d genes x %d subtypes, K_init = %d.",
        nrow(out$x_train), ncol(out$x_train), length(unique(out$z_init))
      )
    )
  }

  out
}
