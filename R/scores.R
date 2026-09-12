# scores.R
# -----------------------------------------------------------------------------
# Prioritisation outputs derived from the dual-significance results.
#
#   gene_score()        rank genes within each module by an expression-weighted
#                       combination of the module's (positive) coefficients
#   cell_type_score()   integrated per-subtype importance summary
#
# Both consume the output of test_significance(): gene_score() uses the
# per-module FGLS coefficients ($coefficients); cell_type_score() uses the
# dual-significance table ($table).
#
# Notation: the coefficient is gamma (the manuscript's gamma_k).
# -----------------------------------------------------------------------------


#' Score and rank genes within each module
#'
#' @description
#' Scores every gene by an expression-weighted combination of its module's
#' coefficients, \eqn{h_{jk} = \sum_p \max(\hat\gamma_{kp}, 0)\, x_{jp}},
#' restricting the sum to positive coefficients so the score reflects only
#' positive contributions to genetic variance, and ranks genes within each
#' module in descending order.
#'
#' This score is an expression-weighted prioritisation, not a direct measure of
#' disease relevance: because it scales with expression magnitude, broadly and
#' highly expressed genes (including housekeeping genes) tend to rank highly.
#' The top-ranked genes should therefore be treated as a permissive candidate
#' set for set-level analysis rather than as individually validated targets, and
#' no conclusion should be drawn from the single highest-ranked gene. Restricting
#' the set to protein-coding genes is left to the user (see the package
#' documentation), as it requires an external gene-type annotation.
#'
#' @param partition A partition from \code{\link{get_partition}} (uses
#'   \code{$assignments}).
#' @param coefficients The per-module coefficient list returned by
#'   \code{\link{test_significance}} as \code{$coefficients} (named by module
#'   id, each a vector over cell subtypes).
#' @param x_train The gene-by-subtype matrix used to fit the model
#'   (\code{thread_data$x_train}); row names are gene identifiers and column order
#'   must match the coefficient vectors.
#' @param top_n Number of top genes to keep per module. Default \code{200}.
#' @param out_dir Optional directory; when supplied, one CSV of the top genes is
#'   written per module. When \code{NULL} (default) no file is written and the
#'   data are only returned.
#' @param prefix File-name prefix used when \code{out_dir} is supplied. Default
#'   \code{"DPM"}.
#' @param verbose Logical; print a short summary. Default \code{TRUE}.
#'
#' @return A list with
#'   \describe{
#'     \item{scores}{named numeric vector of the score of every gene (\code{NA}
#'       for genes in a module that could not be fit).}
#'     \item{top_genes}{data frame with \code{cluster}, \code{rank},
#'       \code{gene} and \code{score}, holding the top \code{top_n} genes per
#'       module.}
#'   }
#'
#' @examples
#' \dontrun{
#' sig <- test_significance(fit, partition, thread_data)
#' gs <- gene_score(partition, sig$coefficients, thread_data$x_train, top_n = 200)
#' head(gs$top_genes)
#' }
#'
#' @importFrom utils write.csv
#' @export
gene_score <- function(partition, coefficients, x_train,
                       top_n = 200, out_dir = NULL, prefix = "DPM",
                       verbose = TRUE) {
  assignments <- partition$assignments
  if (is.null(assignments)) stop("'partition' has no $assignments; use get_partition().")
  if (is.null(coefficients) || length(coefficients) == 0L) {
    stop("'coefficients' is empty; pass test_significance()$coefficients.")
  }

  x_train <- as.matrix(x_train)
  N <- nrow(x_train)
  if (length(assignments) != N) {
    stop(
      "partition$assignments length (", length(assignments),
      ") does not match x_train rows (", N, ")."
    )
  }

  gene_names <- rownames(x_train)
  if (is.null(gene_names)) {
    gene_names <- as.character(seq_len(N))
    rownames(x_train) <- gene_names
  }

  # ---- per-gene expression-weighted score (positive coefficients only) -------
  scores <- numeric(N)
  names(scores) <- gene_names
  for (i in seq_len(N)) {
    cl <- assignments[i]
    gamma_vec <- coefficients[[as.character(cl)]]
    if (is.null(gamma_vec)) {
      scores[i] <- NA_real_
    } else {
      gamma_pos <- pmax(gamma_vec, 0) # clamp negative coefficients to 0
      scores[i] <- sum(gamma_pos * x_train[i, ])
    }
  }

  # ---- top genes per module --------------------------------------------------
  top_list <- list()
  for (cl in sort(unique(assignments))) {
    idx <- which(assignments == cl)
    g <- gene_names[idx]
    s <- scores[idx]
    keep <- !is.na(s)
    g <- g[keep]
    s <- s[keep]
    if (length(s) == 0L) next

    o <- order(s, decreasing = TRUE)
    n_sel <- min(length(g), top_n)
    top_list[[length(top_list) + 1L]] <- data.frame(
      cluster = cl,
      rank = seq_len(n_sel),
      gene = g[o][seq_len(n_sel)],
      score = s[o][seq_len(n_sel)],
      stringsAsFactors = FALSE
    )
  }
  top_genes <- if (length(top_list)) {
    do.call(rbind, top_list)
  } else {
    data.frame(
      cluster = integer(0), rank = integer(0),
      gene = character(0), score = numeric(0), stringsAsFactors = FALSE
    )
  }

  # ---- optional CSV output (only when explicitly requested) ------------------
  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (cl in unique(top_genes$cluster)) {
      sub <- top_genes[top_genes$cluster == cl, c("rank", "gene", "score"), drop = FALSE]
      fn <- file.path(out_dir, sprintf("%s_Cluster%s_TopGenes.csv", prefix, cl))
      utils::write.csv(sub, file = fn, row.names = FALSE)
    }
    if (verbose) {
      DPM_log(
        "Scores",
        sprintf(
          "Wrote top-%d gene tables for %d module(s) to %s.",
          top_n, length(unique(top_genes$cluster)), out_dir
        )
      )
    }

  }

  if (verbose) {
    n_na <- sum(is.na(scores))
    if (n_na > 0L) {
      DPM_log(
        "Scores",
        sprintf("%d gene(s) belong to a module without coefficients (score = NA).", n_na)
      )
    }
  }

  list(scores = scores, top_genes = top_genes)
}


#' Integrated cell-subtype importance score (PLACEHOLDER)
#'
#' @description
#' Summarises the overall importance of each cell subtype across the whole model
#' as the total absolute coefficient across modules,
#' \eqn{S_p = \sum_k |\hat\gamma_{kp}|}. Absolute values are used because, unlike
#' the within-module gene score, a negative coefficient still encodes a genuine
#' association and therefore contributes to a subtype's overall effect. The
#' score is a magnitude summary only; statistical prioritisation is carried by
#' the dual-significance criterion, so a subtype is flagged as significant only
#' when it is dual-significant in at least one module, graded by its strongest
#' evidence (smallest FDR among its dual-significant modules).
#'
#' NOTE: this scoring scheme is a placeholder following the manuscript and is
#' not yet finalised; it may be revised. Only this function need change.
#'
#' @param sig_result The list returned by \code{\link{test_significance}} (uses
#'   \code{$table}).
#'
#' @return A data frame ordered by \code{S} (descending), with columns
#'   \code{subtype}, \code{S} (integrated score), \code{significant} (dual in at
#'   least one module), \code{min_fdr} (smallest FDR among dual-significant
#'   modules, \code{NA} if none) and \code{grade} (0/1/2/3 for
#'   \code{q >= 0.05}, \code{< 0.05}, \code{< 0.01}, \code{< 0.001}).
#'
#' @examples
#' \dontrun{
#' sig <- test_significance(fit, partition, thread_data)
#' cts <- cell_type_score(sig)
#' head(cts)
#' }
#'
#' @export
cell_type_score <- function(sig_result) {
  if (is.null(sig_result$table)) {
    stop("'sig_result' has no $table; pass the output of test_significance().")
  }
  tab <- sig_result$table

  empty <- data.frame(
    subtype = character(0), S = numeric(0),
    significant = logical(0), min_fdr = numeric(0),
    grade = integer(0), stringsAsFactors = FALSE
  )
  if (nrow(tab) == 0L) {
    return(empty)
  }

  subtypes <- unique(tab$subtype)
  rows <- lapply(subtypes, function(st) {
    r <- tab[tab$subtype == st, , drop = FALSE]

    S_p <- sum(abs(r$gamma_fgls), na.rm = TRUE) # S_p = sum_k |gamma_kp|
    is_dual <- isTRUE(any(r$dual_significant %in% TRUE))

    dual_rows <- r[r$dual_significant %in% TRUE, , drop = FALSE]
    min_fdr <- if (nrow(dual_rows) > 0L) min(dual_rows$fdr, na.rm = TRUE) else NA_real_

    grade <- if (!is.na(min_fdr)) {
      if (min_fdr < 0.001) 3L else if (min_fdr < 0.01) 2L else if (min_fdr < 0.05) 1L else 0L
    } else {
      0L
    }

    data.frame(
      subtype = st, S = S_p, significant = is_dual,
      min_fdr = min_fdr, grade = grade, stringsAsFactors = FALSE
    )
  })

  res <- do.call(rbind, rows)
  res <- res[order(res$S, decreasing = TRUE), , drop = FALSE]
  rownames(res) <- NULL
  res
}
