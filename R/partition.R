# partition.R
# -----------------------------------------------------------------------------
# Point estimate of the gene partition from the MCMC label draws.
#
# Strategy:
#   1. Build the posterior dissimilarity matrix 1 - PSM from retained label draws.
#   2. Estimate the trajectory mode of the number of major clusters, where
#      micro-clusters with <= min_cluster_size genes are treated as CRP noise.
#   3. Run PAM on 1 - PSM only in a local K window around the trajectory mode:
#        K in [max(2, mode_k - mode_window), min(N - 1, mode_k + mode_window)].
#   4. Choose the K with the maximum average silhouette width in this local
#      posterior-supported window.
#
# This makes the final partition depend on both the MCMC posterior trajectory
# and the topology of the posterior similarity matrix, without allowing the
# silhouette sweep to select an arbitrary K far from the posterior-supported
# cluster count.
# -----------------------------------------------------------------------------


#' Most frequent value (mode) of a vector
#'
#' @param v A vector.
#' @return The value that occurs most often (the first one on ties).
#' @noRd
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#' Point estimate of the gene partition from a fitted DPM model
#'
#' @description
#' Turns posterior label draws into a single point-estimate partition. Because
#' mixture labels are exchangeable, the function first constructs the posterior
#' dissimilarity matrix \code{1 - PSM}. It then computes the trajectory mode of
#' the number of major clusters, ignoring transient micro-clusters with at most
#' \code{min_cluster_size} genes. The final K is selected by PAM silhouette
#' within a local window around this trajectory mode.
#'
#' @param fit A \code{thread_fit} object from \code{\link{run_thread}}. For backward
#'   compatibility, objects with \code{$results} instead of \code{$samples} are
#'   also accepted.
#' @param min_cluster_size Clusters with at most this many genes are excluded
#'   when counting the number of major modules along the chain. Default
#'   \code{10}.
#' @param mode_window Non-negative integer. Candidate K values are restricted to
#'   \code{mode_k - mode_window} through \code{mode_k + mode_window}, clipped to
#'   the valid range. Default \code{3}.
#' @param min_silhouette Optional minimum average silhouette width. If supplied
#'   and the best silhouette is below this threshold, the result collapses to
#'   \code{K = 1}. Default \code{NULL}.
#' @param return_dissimilarity Logical; whether to include the full gene-by-gene
#'   dissimilarity matrix in the output. This can be large, so the default is
#'   \code{FALSE}.
#' @param verbose Logical; print trajectory and selection diagnostics. Default
#'   \code{TRUE}.
#'
#' @return A list with
#'   \describe{
#'     \item{assignments}{Integer vector of length N giving each gene's module.}
#'     \item{K}{Number of modules in the point estimate.}
#'     \item{silhouette}{List containing candidate K values, silhouette widths,
#'       selected K and maximum width.}
#'     \item{mode_k}{Trajectory mode of the major-cluster count.}
#'     \item{major_cluster_count}{Major-cluster count for every retained draw.}
#'     \item{dissimilarity}{The full \code{1 - PSM} matrix, only when
#'       \code{return_dissimilarity = TRUE}.}
#'   }
#'
#' @examples
#' \dontrun{
#' partition <- get_partition(fit, min_cluster_size = 10, mode_window = 3)
#' partition$K
#' table(partition$assignments)
#' }
#'
#' @importFrom cluster pam
#' @importFrom stats as.dist
#' @importFrom utils tail
#' @export
get_partition <- function(fit,
                          min_cluster_size = 10,
                          mode_window = 3,
                          min_silhouette = NULL,
                          return_dissimilarity = FALSE,
                          verbose = TRUE) {
  # ---- validate --------------------------------------------------------------
  samples <- fit$samples

  # backward compatibility with old scripts
  if (is.null(samples) && !is.null(fit$results)) {
    samples <- fit$results
    if (verbose) {
      DPM_log("Partition", "Using legacy fit$results as posterior samples. Consider converting to fit$samples.")
    }
  }

  if (is.null(samples)) {
    stop("'fit' has no $samples; pass a thread_fit object produced by run_thread().")
  }

  n_samples <- length(samples)
  if (n_samples == 0L) {
    stop("'fit$samples' is empty; nothing to summarise.")
  }

  N <- length(samples[[1]]$z)
  if (N == 0L) {
    stop("The label vectors in 'fit$samples' are empty.")
  }

  if (!is.numeric(min_cluster_size) || length(min_cluster_size) != 1L || min_cluster_size < 0) {
    stop("'min_cluster_size' must be a non-negative scalar.")
  }
  min_cluster_size <- as.integer(min_cluster_size)

  if (!is.numeric(mode_window) || length(mode_window) != 1L || mode_window < 0) {
    stop("'mode_window' must be a non-negative scalar.")
  }
  mode_window <- as.integer(mode_window)

  # ---- 1. assemble the label matrix: draws x genes ---------------------------
  z_matrix <- matrix(NA_integer_, nrow = n_samples, ncol = N)

  for (i in seq_len(n_samples)) {
    zi <- samples[[i]]$z
    if (length(zi) != N) {
      stop(sprintf(
        "Draw %d has %d labels but draw 1 has %d; inconsistent sample sizes.",
        i, length(zi), N
      ))
    }
    z_matrix[i, ] <- as.integer(zi)
  }

  # ---- 2. trajectory mode of major clusters ---------------------------------
  major_clusters_count <- integer(n_samples)

  for (i in seq_len(n_samples)) {
    cluster_sizes <- table(z_matrix[i, ])
    n_garbage <- sum(cluster_sizes <= min_cluster_size)
    major_clusters_count[i] <- length(cluster_sizes) - n_garbage
  }

  mode_major_k <- as.integer(get_mode(major_clusters_count))

  mode_freq_total <- sum(major_clusters_count == mode_major_k) / n_samples * 100
  last_n <- min(1000L, n_samples)
  last_samples <- utils::tail(major_clusters_count, last_n)
  mode_freq_last <- sum(last_samples == mode_major_k) / last_n * 100

  if (verbose) {
    DPM_log(
      "Partition",
      sprintf(
        "Trajectory mode of major clusters: %d (overall %.1f%%, last %d draws %.1f%%).",
        mode_major_k, mode_freq_total, last_n, mode_freq_last
      )
    )
  }

  # ---- 3. posterior dissimilarity matrix: 1 - PSM ----------------------------
  dissim_matrix <- calcDissimilarityMatrixCpp(z_matrix)

  collapse_to_one <- function(sil = NULL) {
    assignments <- rep(1L, N)

    gene_names <- NULL
    if (!is.null(fit$dims$gene_names)) gene_names <- fit$dims$gene_names
    if (!is.null(gene_names) && length(gene_names) == N) names(assignments) <- gene_names

    out <- list(
      assignments = assignments,
      K = 1L,
      silhouette = sil,
      mode_k = mode_major_k,
      major_cluster_count = major_clusters_count
    )

    if (return_dissimilarity) out$dissimilarity <- dissim_matrix
    out
  }

  if (mode_major_k <= 1L) {
    if (verbose) {
      DPM_log("Partition", "Trajectory collapsed onto a single major cluster: returning K = 1.")
    }
    return(collapse_to_one())
  }

  if (N < 3L) {
    warning("Too few genes to form a stable multi-cluster partition; returning K = 1.")
    return(collapse_to_one())
  }

  # ---- 4. local silhouette sweep around mode_k -------------------------------
  dist_mat <- stats::as.dist(dissim_matrix)

  lower_k <- max(2L, mode_major_k - mode_window)

  upper_target <- mode_major_k + mode_window
  upper_k <- min(N - 1L, upper_target)

  if (upper_k < lower_k) {
    upper_k <- lower_k
  }

  k_grid <- seq.int(lower_k, upper_k)

  if (verbose) {
    DPM_log(
      "Partition",
      sprintf(
        "Silhouette K search window: %s (mode_k = %d, mode_window = %d).",
        paste(k_grid, collapse = ", "), mode_major_k, mode_window
      )
    )
  }

  sil_scores <- vapply(k_grid, function(k) {
    pm <- tryCatch(
      cluster::pam(dist_mat, k, diss = TRUE),
      error = function(e) NULL
    )

    if (is.null(pm) || is.null(pm$silinfo)) {
      return(-Inf)
    }

    pm$silinfo$avg.width
  }, numeric(1))

  if (all(!is.finite(sil_scores))) {
    warning("No valid clustering could be evaluated in the mode-local K window; returning K = 1.")
    return(collapse_to_one())
  }

  best_idx <- which.max(sil_scores)
  best_k <- as.integer(k_grid[best_idx])
  max_sil <- sil_scores[best_idx]

  sil_out <- list(
    k = k_grid,
    width = sil_scores,
    best_k = best_k,
    max_width = max_sil,
    mode_k = mode_major_k,
    mode_window = mode_window
  )

  if (!is.null(min_silhouette) && max_sil < min_silhouette) {
    if (verbose) {
      DPM_log(
        "Partition",
        sprintf(
          "Max silhouette %.3f < min_silhouette %.3f; collapsing to K = 1.",
          max_sil, min_silhouette
        )
      )
    }
    return(collapse_to_one(sil = sil_out))
  }

  if (verbose) {
    DPM_log(
      "Partition",
      sprintf(
        "Selected K = %d by maximum average silhouette width within the mode-local window (%.3f).",
        best_k, max_sil
      )
    )
  }

  pm <- cluster::pam(dist_mat, best_k, diss = TRUE)
  assignments <- as.integer(pm$clustering)

  gene_names <- NULL
  if (!is.null(fit$dims$gene_names)) gene_names <- fit$dims$gene_names
  if (!is.null(gene_names) && length(gene_names) == N) names(assignments) <- gene_names

  out <- list(
    assignments = assignments,
    K = length(unique(assignments)),
    silhouette = sil_out,
    mode_k = mode_major_k,
    major_cluster_count = major_clusters_count
  )

  if (return_dissimilarity) out$dissimilarity <- dissim_matrix

  out
}
