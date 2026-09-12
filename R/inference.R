# inference.R
# -----------------------------------------------------------------------------
# Significance of module / cell-subtype associations, combining a frequentist
# and a Bayesian criterion on the point-estimate partition.
#
#   test_significance()   public entry point; returns the dual-significance
#                         table plus the per-module FGLS coefficients used for
#                         gene scoring
#
# Internal helpers (not exported):
#   .ols_cluster()        ordinary least squares with sandwich SEs
#   .fgls_cluster()       feasible GLS: iteratively reweight by 1 / v_jk, with
#                         heteroscedasticity-consistent (sandwich) SEs
#   .fit_clusters_fgls()  per-module refit + one-sided test + BH FDR
#   .align_gamma_draws()  align posterior gamma draws to the point-estimate
#                         modules by majority vote (handles label switching)
#   .compute_lfsr()       local false sign rate (positive direction) per module
#                         / subtype, with HPD intervals
#
# A module / cell-subtype association is reported as dual-significant only when
# it satisfies both criteria: FDR < fdr_thr (frequentist) and lFSR < lfsr_thr
# (Bayesian). Because the Bayesian posterior mean of gamma is shrunk toward zero
# by the Lasso prior, the FGLS estimate (gamma_fgls) is reported as the effect
# size; the shrunk posterior mean is also returned for reference.
#
# Notation: the coefficient is gamma (the manuscript's gamma_k); the GWAS SNP
# effect named `beta` lives only in gwas.R and is unrelated.
# -----------------------------------------------------------------------------


#' Ordinary least squares for one module, with sandwich standard errors
#'
#' @param X_k Module design matrix (genes by subtypes).
#' @param Y_k Module response vector.
#' @return List with \code{gamma_hat}, \code{se} (robust) and \code{df_resid}.
#' @noRd
ols_cluster <- function(X_k, Y_k) {
  p_k <- ncol(X_k)
  n_k <- nrow(X_k)
  XtX <- crossprod(X_k)
  XtY <- crossprod(X_k, Y_k)

  gamma_hat <- tryCatch(
    as.vector(solve(XtX + 1e-12 * diag(p_k), XtY)),
    error = function(e) rep(0, p_k)
  )

  resid <- Y_k - as.vector(X_k %*% gamma_hat)
  df_resid <- max(n_k - p_k, 1)
  XtX_inv <- tryCatch(solve(XtX + 1e-12 * diag(p_k)),
    error = function(e) matrix(0, p_k, p_k)
  )

  meat_hc <- crossprod(X_k * resid)
  hc_vcov <- XtX_inv %*% meat_hc %*% XtX_inv * (n_k / df_resid)
  se_robust <- sqrt(pmax(diag(hc_vcov), 0))

  list(gamma_hat = gamma_hat, se = se_robust, df_resid = df_resid)
}


#' Feasible GLS for one module (iteratively reweighted by the derived variance)
#'
#' @param X_k Module design matrix (genes by subtypes).
#' @param Y_k Module response vector.
#' @param d0_k,d1_k,d2_k Module heteroscedastic-variance coefficients.
#' @param max_iter Maximum reweighting iterations. Default \code{100}.
#' @param tol Convergence tolerance on the coefficient change. Default
#'   \code{1e-10}.
#' @return List with \code{gamma_hat}, \code{se} (sandwich) and \code{df_resid}.
#' @noRd
fgls_cluster <- function(X_k, Y_k, d0_k, d1_k, d2_k, max_iter = 100, tol = 1e-10) {
  p_k <- ncol(X_k)
  n_k <- nrow(X_k)

  XtX <- crossprod(X_k)
  XtY <- crossprod(X_k, Y_k)
  gamma_hat <- tryCatch(as.vector(solve(XtX + 1e-12 * diag(p_k), XtY)),
    error = function(e) rep(0, p_k)
  )

  for (iter in seq_len(max_iter)) {
    gamma_old <- gamma_hat
    tau_hat <- as.vector(X_k %*% gamma_hat)
    v_hat <- pmax(d0_k + d1_k * tau_hat + d2_k * tau_hat^2, 1e-14)
    sqw <- sqrt(1.0 / v_hat)
    Xw <- X_k * sqw
    Yw <- Y_k * sqw
    gamma_hat <- tryCatch(
      as.vector(solve(crossprod(Xw) + 1e-12 * diag(p_k), crossprod(Xw, Yw))),
      error = function(e) gamma_old
    )
    if (max(abs(gamma_hat - gamma_old)) < tol) break
  }

  tau_final <- as.vector(X_k %*% gamma_hat)
  v_final <- pmax(d0_k + d1_k * tau_final + d2_k * tau_final^2, 1e-14)
  w_final <- 1.0 / v_final
  resid <- Y_k - tau_final
  df_resid <- max(n_k - p_k, 1)

  sqw_f <- sqrt(w_final)
  Xw_f <- X_k * sqw_f
  XtWX_inv <- tryCatch(solve(crossprod(Xw_f) + 1e-12 * diag(p_k)),
    error = function(e) matrix(0, p_k, p_k)
  )
  meat_mat <- crossprod(X_k * (w_final * resid))
  sandwich_vcov <- XtWX_inv %*% meat_mat %*% XtWX_inv
  se_sandwich <- sqrt(pmax(diag(sandwich_vcov), 0))

  list(gamma_hat = gamma_hat, se = se_sandwich, df_resid = df_resid)
}


#' Refit every module and run the one-sided FDR test
#'
#' @param x_train Gene-by-subtype matrix.
#' @param Y Response vector (same scale used to fit the model).
#' @param d0,d1,d2 Heteroscedastic-variance coefficient vectors.
#' @param assignments Integer module label per gene.
#' @param method Either "FGLS" (default) or "OLS".
#' @param fdr_thr FDR threshold for the frequentist call. Default \code{0.05}.
#' @return List with \code{table} (one row per module / subtype) and
#'   \code{coefficients} (named list of per-module \code{gamma_hat}).
#' @importFrom stats pt p.adjust
#' @noRd
fit_clusters_fgls <- function(x_train,
                              Y,
                              d0,
                              d1,
                              d2,
                              assignments,
                              method = "FGLS",
                              fdr_thr = 0.05,
                              verbose = TRUE) {
  method <- toupper(method)

  if (!method %in% c("FGLS", "OLS")) {
    stop("Unknown regression method: choose 'FGLS' or 'OLS'.")
  }

  subtype_names <- colnames(x_train)

  if (is.null(subtype_names)) {
    subtype_names <- paste0("V", seq_len(ncol(x_train)))
  }

  p_k <- ncol(x_train)
  clusters <- sort(unique(assignments))

  rows <- list()
  coefficients <- list()
  module_diagnostics <- list()

  for (cl in clusters) {
    idx <- which(assignments == cl)
    n_k <- length(idx)

    if (n_k == 0L) {
      next
    }

    X_c <- x_train[idx, , drop = FALSE]
    Y_c <- Y[idx]

    design_rank <- qr(X_c, tol = 1e-10)$rank

    fgls_estimable <- (
      n_k > p_k &&
        design_rank == p_k
    )

    fit <- switch(method,
      "FGLS" = fgls_cluster(
        X_k = X_c,
        Y_k = Y_c,
        d0_k = d0[idx],
        d1_k = d1[idx],
        d2_k = d2[idx]
      ),
      "OLS" = ols_cluster(
        X_k = X_c,
        Y_k = Y_c
      )
    )

    gamma_hat <- as.numeric(fit$gamma_hat)
    names(gamma_hat) <- subtype_names

    coefficients[[as.character(cl)]] <- gamma_hat

    if (!fgls_estimable && verbose) {
      reason <- if (n_k <= p_k) {
        sprintf(
          "%d genes and %d predictors",
          n_k,
          p_k
        )
      } else {
        sprintf(
          "design rank %d is smaller than %d predictors",
          design_rank,
          p_k
        )
      }

      message(
        sprintf(
          paste0(
            "[DPM:Inference] Module %s has %s; ",
            "frequentist inference is underdetermined. ",
            "Coefficients are retained, but SE, P value and FDR are set to NA."
          ),
          cl,
          reason
        )
      )
    }

    if (fgls_estimable) {
      se_out <- as.numeric(fit$se)

      p_out <- vapply(
        seq_along(subtype_names),
        function(j) {
          if (!is.finite(gamma_hat[j]) ||
            !is.finite(se_out[j]) ||
            se_out[j] <= 0) {
            return(NA_real_)
          }

          t_stat <- gamma_hat[j] / se_out[j]

          stats::pt(
            t_stat,
            df = fit$df_resid,
            lower.tail = FALSE
          )
        },
        numeric(1L)
      )
    } else {
      se_out <- rep(NA_real_, p_k)
      p_out <- rep(NA_real_, p_k)
    }

    rows[[length(rows) + 1L]] <- data.frame(
      cluster = rep(cl, p_k),
      subtype = subtype_names,
      module_n_genes = rep(n_k, p_k),
      n_predictors = rep(p_k, p_k),
      design_rank = rep(design_rank, p_k),
      fgls_estimable = rep(fgls_estimable, p_k),
      gamma_fgls = gamma_hat,
      se_fgls = se_out,
      pval_fgls = p_out,
      stringsAsFactors = FALSE
    )

    module_diagnostics[[length(module_diagnostics) + 1L]] <- data.frame(
      cluster = cl,
      module_n_genes = n_k,
      n_predictors = p_k,
      design_rank = design_rank,
      fgls_estimable = fgls_estimable,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    empty <- data.frame(
      cluster = integer(0),
      subtype = character(0),
      module_n_genes = integer(0),
      n_predictors = integer(0),
      design_rank = integer(0),
      fgls_estimable = logical(0),
      gamma_fgls = numeric(0),
      se_fgls = numeric(0),
      pval_fgls = numeric(0),
      fdr = numeric(0),
      sig_fdr = logical(0),
      stringsAsFactors = FALSE
    )

    return(
      list(
        table = empty,
        coefficients = coefficients,
        module_diagnostics = data.frame()
      )
    )
  }

  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  valid_test <- (
    tab$fgls_estimable &
      is.finite(tab$pval_fgls)
  )

  tab$fdr <- NA_real_

  if (any(valid_test)) {
    # Count non-estimable comparisons conservatively as P = 1
    # in the overall multiple-testing family.
    p_for_adjustment <- rep(1, nrow(tab))
    p_for_adjustment[valid_test] <- tab$pval_fgls[valid_test]

    adjusted_all <- stats::p.adjust(
      p_for_adjustment,
      method = "BH"
    )

    tab$fdr[valid_test] <- adjusted_all[valid_test]
  }

  tab$sig_fdr <- (
    valid_test &
      !is.na(tab$fdr) &
      tab$fdr < fdr_thr &
      tab$gamma_fgls > 0
  )

  module_diagnostics <- do.call(
    rbind,
    module_diagnostics
  )

  rownames(module_diagnostics) <- NULL

  if (verbose) {
    n_estimable <- sum(
      module_diagnostics$fgls_estimable
    )

    message(
      sprintf(
        paste0(
          "[DPM:Inference] Frequentist inference was estimable for ",
          "%d of %d module(s)."
        ),
        n_estimable,
        nrow(module_diagnostics)
      )
    )
  }

  list(
    table = tab,
    coefficients = coefficients,
    module_diagnostics = module_diagnostics
  )
}


#' Align posterior gamma draws to the point-estimate modules
#'
#' @description
#' Because module labels switch across MCMC iterations, each draw's coefficient
#' rows are matched to the point-estimate modules by majority vote: for every
#' point-estimate module, the draw's most frequent label among that module's
#' genes selects which coefficient row to take.
#'
#' @param assignments Integer module label per gene (point estimate).
#' @param fit A \code{thread_fit}; each draw in \code{$samples} carries \code{$z}
#'   and \code{$gamma}.
#' @return List with \code{gamma_aligned} (a draws-by-module-by-subtype array)
#'   and \code{ref_clusters} (the module labels in row order).
#' @noRd
align_gamma_draws <- function(assignments, fit) {
  ref_clusters <- sort(unique(assignments))
  C <- length(ref_clusters)
  optAlloc <- lapply(ref_clusters, function(cc) which(assignments == cc))

  samples <- fit$samples
  Tsave <- length(samples)
  p <- ncol(samples[[1]]$gamma)
  gamma_aligned <- array(NA_real_, dim = c(Tsave, C, p))

  for (t in seq_len(Tsave)) {
    z_t <- samples[[t]]$z
    G_t <- samples[[t]]$gamma
    if (is.null(rownames(G_t))) rownames(G_t) <- as.character(seq_len(nrow(G_t)))

    for (c_i in seq_len(C)) {
      idx <- optAlloc[[c_i]]
      if (length(idx) == 0L) next
      label_counts <- table(z_t[idx])
      best_match_k <- names(which.max(label_counts))
      if (best_match_k %in% rownames(G_t)) {
        gamma_aligned[t, c_i, ] <- G_t[best_match_k, ]
      }
    }
  }

  list(gamma_aligned = gamma_aligned, ref_clusters = ref_clusters)
}


#' Local false sign rate (positive direction) per module / subtype
#'
#' @param gamma_aligned Draws-by-module-by-subtype array from
#'   \code{.align_gamma_draws}.
#' @param ref_clusters Module labels in array-row order.
#' @param subtype_names Cell-subtype names.
#' @param lfsr_thr lFSR threshold for the Bayesian call. Default \code{0.05}.
#' @param cred_mass Credible mass for the HPD interval. Default \code{0.95}.
#' @return Data frame with one row per module / subtype.
#' @importFrom HDInterval hdi
#' @noRd
compute_lfsr <- function(gamma_aligned, ref_clusters, subtype_names,
                          lfsr_thr = 0.05, cred_mass = 0.95) {
  nCluster <- dim(gamma_aligned)[2]
  nFeature <- dim(gamma_aligned)[3]
  if (is.null(subtype_names)) subtype_names <- as.character(seq_len(nFeature))

  rows <- list()
  for (k in seq_len(nCluster)) {
    for (j in seq_len(nFeature)) {
      s <- gamma_aligned[, k, j]
      if (all(is.na(s))) next

      # lFSR for a positive effect: posterior probability that gamma <= 0
      lfsr <- mean(s <= 0, na.rm = TRUE)
      pm <- mean(s, na.rm = TRUE)
      hpd <- tryCatch(HDInterval::hdi(s, credMass = cred_mass),
        error = function(e) c(NA_real_, NA_real_)
      )

      rows[[length(rows) + 1L]] <- data.frame(
        cluster = ref_clusters[k],
        subtype = subtype_names[j],
        gamma_post_mean = pm,
        hpd_lower = hpd[1],
        hpd_upper = hpd[2],
        lfsr = lfsr,
        sig_lfsr = (lfsr < lfsr_thr) & (pm > 0),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      cluster = integer(0), subtype = character(0),
      gamma_post_mean = numeric(0), hpd_lower = numeric(0),
      hpd_upper = numeric(0), lfsr = numeric(0),
      sig_lfsr = logical(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}


#' Dual-significance testing of module / cell-subtype associations
#'
#' @description
#' Evaluates each module / cell-subtype association under two independent
#' criteria on the point-estimate partition and reports those supported by
#' both. The frequentist criterion refits each module by feasible GLS
#' (reweighting genes by the precision of their own estimate, with sandwich
#' standard errors), runs a one-sided test for positive enrichment, and
#' controls the FDR by Benjamini-Hochberg. The Bayesian criterion aligns the
#' posterior coefficient draws to the modules and computes the local false sign
#' rate. An association is dual-significant when \code{fdr < fdr_thr} and
#' \code{lfsr < lfsr_thr}. The reported effect size is the FGLS coefficient,
#' since the Bayesian posterior mean is shrunk by the Lasso prior.
#'
#' @param fit A \code{thread_fit} from \code{\link{run_thread}}; each draw carries
#'   \code{$z} and \code{$gamma}.
#' @param partition A partition from \code{\link{get_partition}} (uses
#'   \code{$assignments}).
#' @param thread_data The \code{thread_data} used to fit the model; supplies
#'   \code{x_train}, \code{y_train} and \code{d0}, \code{d1}, \code{d2} on the
#'   same scale the sampler used.
#' @param method Frequentist engine, "FGLS" (default) or "OLS".
#' @param fdr_thr FDR threshold. Default \code{0.05}.
#' @param lfsr_thr lFSR threshold. Default \code{0.05}.
#' @param cred_mass Credible mass for HPD intervals. Default \code{0.95}.
#' @param verbose Logical; print progress. Default \code{TRUE}.
#'
#' @return A list with the following elements:
#'   \describe{
#'     \item{table}{A data frame with one row per module and cellular
#'       subtype. It contains module dimensions, FGLS estimates, standard
#'       errors, one-sided P values, FDR values, posterior summaries, local
#'       false sign rates, and the dual-significance indicator.}
#'     \item{coefficients}{A named list of per-module FGLS coefficient
#'       vectors. These coefficients are used by \code{\link{gene_score}}.}
#'     \item{module_diagnostics}{A data frame with one row per module,
#'       reporting the number of genes, number of predictors, design rank,
#'       and whether frequentist inference was estimable.}
#'     \item{method}{The frequentist inference method used.}
#'     \item{thresholds}{A list containing the FDR threshold, lFSR threshold,
#'       and credible mass.}
#'   }
#'
#' @examples
#' \dontrun{
#' sig <- test_significance(fit, partition, thread_data)
#' subset(sig$table, dual_significant)
#' }
#'
#' @export
test_significance <- function(fit,
                              partition,
                              thread_data,
                              method = "FGLS",
                              fdr_thr = 0.05,
                              lfsr_thr = 0.05,
                              cred_mass = 0.95,
                              verbose = TRUE) {
  # ---- validate posterior draws ---------------------------------------------
  if (is.null(fit$samples) ||
    length(fit$samples) == 0L) {
    stop("'fit' has no posterior draws in $samples.")
  }

  if (is.null(fit$samples[[1L]]$gamma)) {
    stop(
      "'fit$samples' draws do not contain $gamma; ",
      "was the model fit with run_thread()?"
    )
  }

  assignments <- partition$assignments

  if (is.null(assignments)) {
    stop(
      "'partition' has no $assignments; use get_partition()."
    )
  }

  required_fields <- c(
    "x_train",
    "y_train",
    "d0",
    "d1",
    "d2"
  )

  missing_fields <- setdiff(
    required_fields,
    names(thread_data)
  )

  if (length(missing_fields) > 0L) {
    stop(
      "'thread_data' is missing field(s): ",
      paste(missing_fields, collapse = ", ")
    )
  }

  if (!is.numeric(fdr_thr) ||
    length(fdr_thr) != 1L ||
    is.na(fdr_thr) ||
    fdr_thr <= 0 ||
    fdr_thr >= 1) {
    stop("'fdr_thr' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(lfsr_thr) ||
    length(lfsr_thr) != 1L ||
    is.na(lfsr_thr) ||
    lfsr_thr <= 0 ||
    lfsr_thr >= 1) {
    stop("'lfsr_thr' must lie strictly between 0 and 1.")
  }

  if (!is.numeric(cred_mass) ||
    length(cred_mass) != 1L ||
    is.na(cred_mass) ||
    cred_mass <= 0 ||
    cred_mass >= 1) {
    stop("'cred_mass' must lie strictly between 0 and 1.")
  }

  x_train <- as.matrix(thread_data$x_train)
  Y <- as.numeric(thread_data$y_train)
  d0 <- as.numeric(thread_data$d0)
  d1 <- as.numeric(thread_data$d1)
  d2 <- as.numeric(thread_data$d2)

  N <- nrow(x_train)

  if (length(assignments) != N) {
    stop(
      "partition$assignments length (",
      length(assignments),
      ") does not match the number of genes (",
      N,
      ")."
    )
  }

  if (length(Y) != N ||
    length(d0) != N ||
    length(d1) != N ||
    length(d2) != N) {
    stop(
      "y_train/d0/d1/d2 lengths must match nrow(x_train)."
    )
  }

  subtype_names <- colnames(x_train)

  if (is.null(subtype_names)) {
    subtype_names <- fit$dims$subtype_names
  }

  if (is.null(subtype_names)) {
    subtype_names <- paste0(
      "Subtype",
      seq_len(ncol(x_train))
    )
  }

  # ---- frequentist inference -----------------------------------------------
  if (verbose) {
    message(
      sprintf(
        "[DPM:Inference] Frequentist inference (%s) over %d module(s).",
        toupper(method),
        length(unique(assignments))
      )
    )
  }

  freq <- fit_clusters_fgls(
    x_train = x_train,
    Y = Y,
    d0 = d0,
    d1 = d1,
    d2 = d2,
    assignments = assignments,
    method = method,
    fdr_thr = fdr_thr,
    verbose = verbose
  )

  # ---- Bayesian inference ---------------------------------------------------
  if (verbose) {
    message(
      paste0(
        "[DPM:Inference] Bayesian inference (lFSR): ",
        "aligning posterior draws to the partition."
      )
    )
  }

  aligned <- align_gamma_draws(
    assignments = assignments,
    fit = fit
  )

  bayes <- compute_lfsr(
    gamma_aligned = aligned$gamma_aligned,
    ref_clusters = aligned$ref_clusters,
    subtype_names = subtype_names,
    lfsr_thr = lfsr_thr,
    cred_mass = cred_mass
  )

  # ---- combine criteria -----------------------------------------------------
  tab <- merge(
    freq$table,
    bayes,
    by = c("cluster", "subtype"),
    all = TRUE,
    sort = FALSE
  )

  tab$dual_significant <- (
    tab$sig_fdr %in% TRUE &
      tab$sig_lfsr %in% TRUE
  )

  column_order <- c(
    "cluster",
    "subtype",
    "module_n_genes",
    "n_predictors",
    "design_rank",
    "fgls_estimable",
    "gamma_fgls",
    "se_fgls",
    "pval_fgls",
    "fdr",
    "sig_fdr",
    "gamma_post_mean",
    "hpd_lower",
    "hpd_upper",
    "lfsr",
    "sig_lfsr",
    "dual_significant"
  )

  tab <- tab[
    ,
    intersect(column_order, names(tab)),
    drop = FALSE
  ]

  tab <- tab[
    order(tab$cluster, tab$subtype), ,
    drop = FALSE
  ]

  rownames(tab) <- NULL

  if (verbose) {
    n_dual <- sum(
      tab$dual_significant,
      na.rm = TRUE
    )

    message(
      sprintf(
        paste0(
          "[DPM:Inference] Dual-significant ",
          "(FDR < %.3g and lFSR < %.3g): ",
          "%d module-subtype association(s)."
        ),
        fdr_thr,
        lfsr_thr,
        n_dual
      )
    )
  }

  list(
    table = tab,
    coefficients = freq$coefficients,
    module_diagnostics = freq$module_diagnostics,
    method = toupper(method),
    thresholds = list(
      fdr = fdr_thr,
      lfsr = lfsr_thr,
      cred_mass = cred_mass
    )
  )
}
