# LD-based heteroscedastic gene-level response precomputation.
#
# This implementation provides one unified calculation path for:
#   final_y, d0, d1, d2 and n_snps.
#
# Variant-key mode is user-selectable:
#   1. match_mode = "snp_name" (default and recommended)
#   2. match_mode = "snp_pos"  (fallback for GWAS files with unreliable SNP IDs)
#
# IMPORTANT:
# - Chromosome-level LD files may store each off-diagonal SNP pair in only one
#   direction. This implementation therefore ALWAYS recovers the reciprocal
#   direction and constructs a symmetric gene-level LD matrix.
# - There is intentionally no switch for disabling LD symmetrization.
# - In position mode, GWAS matching is performed by chromosome + position.
#   Position-only matching across the whole genome is not allowed.
# - Missing LD pairs in both directions are treated as zero, consistent with a
#   sparse LD representation. Their fraction is retained as a diagnostic.
# - If the symmetric sparse LD matrix is not positive semidefinite, negative
#   eigenvalues are clipped to zero and the retained eigenvalues are rescaled
#   to sum to the SNP count. This preserves tr(R) = k, which is required by the
#   response derivation.

validate_compute_response_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, names(x))

  if (length(missing_cols) > 0L) {
    stop(
      "'", object_name, "' is missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}

resolve_response_column <- function(data,
                                    candidates,
                                    user_col = NULL,
                                    arg_name = "column") {
  if (!is.null(user_col)) {
    if (!is.character(user_col) ||
      length(user_col) != 1L ||
      is.na(user_col)) {
      stop("'", arg_name, "' must be a single column name.")
    }

    if (!user_col %in% names(data)) {
      stop(
        "Column '", user_col,
        "' specified by '", arg_name,
        "' was not found."
      )
    }

    return(user_col)
  }

  lower_names <- tolower(names(data))
  lower_candidates <- tolower(candidates)
  idx <- match(lower_candidates, lower_names)
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    stop(
      "Could not infer ", arg_name,
      ". Tried candidate columns: ",
      paste(candidates, collapse = ", "),
      "."
    )
  }

  names(data)[idx[1L]]
}

load_table_or_object <- function(x, object_name = "object") {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) {
      stop("File for '", object_name, "' does not exist: ", x)
    }

    if (grepl("\\.rds$", x, ignore.case = TRUE)) {
      return(readRDS(x))
    }

    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required. Please install it first.")
    }

    return(
      data.table::fread(
        x,
        data.table = FALSE,
        showProgress = FALSE
      )
    )
  }

  x
}

make_ld_file_path <- function(ld_dir,
                              chr,
                              ld_file_template = "{chr}.rds") {
  file_name <- gsub(
    "\\{chr\\}",
    as.character(chr),
    ld_file_template,
    fixed = FALSE
  )

  file.path(ld_dir, file_name)
}

normalize_chr_response <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x
}

normalize_position_key <- function(x) {
  x_num <- suppressWarnings(as.integer(x))
  as.character(x_num)
}

# Remove duplicate GWAS rows safely.
#
# Exact duplicates are collapsed. If one key still has multiple distinct
# Beta/SE combinations, that key is ambiguous and is removed rather than
# silently retaining an arbitrary record.
collapse_gwas_effect_keys <- function(gwas_eff,
                                      key_cols,
                                      verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  gwas_eff <- data.table::as.data.table(gwas_eff)

  exact_cols <- c(key_cols, "Beta", "SE")
  gwas_eff <- unique(gwas_eff, by = exact_cols)

  key_counts <- gwas_eff[, .N, by = key_cols]
  ambiguous_keys <- key_counts[N > 1L]

  if (nrow(ambiguous_keys) > 0L) {
    if (verbose) {
      warning(
        nrow(ambiguous_keys),
        " ambiguous GWAS key(s) had multiple distinct Beta/SE records and ",
        "were removed."
      )
    }

    ambiguous_keys[, ambiguous_key := TRUE]

    gwas_eff <- ambiguous_keys[
      gwas_eff,
      on = key_cols
    ][is.na(ambiguous_key)]

    gwas_eff[, ambiguous_key := NULL]

    if ("N" %in% names(gwas_eff)) {
      gwas_eff[, N := NULL]
    }
  }

  gwas_eff
}

standardize_matches_for_response <- function(matches,
                                             match_mode = c(
                                               "snp_name",
                                               "snp_pos"
                                             ),
                                             snp_col = "snp_name",
                                             snp_pos_col = "snp_pos",
                                             gene_col = "gene_name",
                                             chr_col = "chrom") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  match_mode <- match.arg(match_mode)
  matches <- data.table::as.data.table(matches)

  validate_compute_response_columns(
    matches,
    c(gene_col, chr_col),
    "matches"
  )

  out <- data.table::copy(matches)

  if (gene_col != "gene_name") {
    data.table::setnames(out, gene_col, "gene_name")
  }

  if (chr_col != "chrom") {
    data.table::setnames(out, chr_col, "chrom")
  }

  out[, gene_name := as.character(gene_name)]
  out[, chrom := normalize_chr_response(chrom)]

  if (match_mode == "snp_name") {
    validate_compute_response_columns(out, snp_col, "matches")

    if (snp_col != "snp_name") {
      data.table::setnames(out, snp_col, "snp_name")
    }

    out[, snp_name := as.character(snp_name)]
    out <- out[!is.na(snp_name) & snp_name != ""]
    out[, key_value := snp_name]
  } else {
    validate_compute_response_columns(out, snp_pos_col, "matches")

    if (snp_pos_col != "snp_pos") {
      data.table::setnames(out, snp_pos_col, "snp_pos")
    }

    out[, snp_pos := suppressWarnings(as.integer(snp_pos))]
    out <- out[!is.na(snp_pos) & snp_pos > 0L]
    out[, key_value := normalize_position_key(snp_pos)]
  }

  out <- out[
    !is.na(gene_name) & gene_name != "" &
      !is.na(chrom) & chrom != "" &
      !is.na(key_value) & key_value != ""
  ]

  # Multiple rows for the same gene and key do not add information for the
  # gene-level response and would otherwise duplicate one variant.
  out <- unique(out, by = c("gene_name", "chrom", "key_value"))

  out
}

merge_gwas_effects_for_response <- function(matches,
                                            gwas_clean = NULL,
                                            match_mode = c(
                                              "snp_name",
                                              "snp_pos"
                                            ),
                                            gwas_snp_col = "rsid",
                                            gwas_pos_col = NULL,
                                            gwas_chr_col = "chrom",
                                            beta_col = "beta",
                                            se_col = "se",
                                            verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  match_mode <- match.arg(match_mode)
  matches <- data.table::as.data.table(matches)

  if (!is.null(gwas_clean)) {
    gwas_clean <- load_table_or_object(
      gwas_clean,
      object_name = "gwas_clean"
    )

    gwas_clean <- data.table::as.data.table(gwas_clean)

    validate_compute_response_columns(
      gwas_clean,
      c(beta_col, se_col),
      "gwas_clean"
    )

    if (match_mode == "snp_name") {
      validate_compute_response_columns(
        gwas_clean,
        gwas_snp_col,
        "gwas_clean"
      )

      gwas_eff <- gwas_clean[
        ,
        .(
          key_value = as.character(get(gwas_snp_col)),
          Beta = suppressWarnings(as.numeric(get(beta_col))),
          SE = suppressWarnings(as.numeric(get(se_col)))
        )
      ]

      gwas_eff <- gwas_eff[
        !is.na(key_value) & key_value != "" &
          is.finite(Beta) & is.finite(SE) &
          SE > 0
      ]

      gwas_eff <- collapse_gwas_effect_keys(
        gwas_eff,
        key_cols = "key_value",
        verbose = verbose
      )

      if ("Beta" %in% names(matches)) {
        matches[, Beta := NULL]
      }

      if ("SE" %in% names(matches)) {
        matches[, SE := NULL]
      }

      matches <- gwas_eff[
        matches,
        on = "key_value"
      ]
    } else {
      # Position matching is only valid when chromosome is included.
      if (is.null(gwas_chr_col) ||
        !gwas_chr_col %in% names(gwas_clean)) {
        stop(
          "For match_mode = 'snp_pos', the GWAS file must contain a ",
          "chromosome column specified by 'gwas_chr_col'. Position-only ",
          "matching across the whole genome is not allowed."
        )
      }

      if (is.null(gwas_pos_col)) {
        gwas_pos_col <- resolve_response_column(
          gwas_clean,
          candidates = c(
            "pos_hg19",
            "snp_pos",
            "pos",
            "bp",
            "BP",
            "position"
          ),
          user_col = NULL,
          arg_name = "gwas_pos_col"
        )
      } else {
        validate_compute_response_columns(
          gwas_clean,
          gwas_pos_col,
          "gwas_clean"
        )
      }

      gwas_eff <- gwas_clean[
        ,
        .(
          chrom = normalize_chr_response(get(gwas_chr_col)),
          key_value = normalize_position_key(get(gwas_pos_col)),
          Beta = suppressWarnings(as.numeric(get(beta_col))),
          SE = suppressWarnings(as.numeric(get(se_col)))
        )
      ]

      gwas_eff <- gwas_eff[
        !is.na(chrom) & chrom != "" &
          !is.na(key_value) & key_value != "" &
          is.finite(Beta) & is.finite(SE) &
          SE > 0
      ]

      gwas_eff <- collapse_gwas_effect_keys(
        gwas_eff,
        key_cols = c("chrom", "key_value"),
        verbose = verbose
      )

      if ("Beta" %in% names(matches)) {
        matches[, Beta := NULL]
      }

      if ("SE" %in% names(matches)) {
        matches[, SE := NULL]
      }

      matches <- gwas_eff[
        matches,
        on = c("chrom", "key_value")
      ]
    }

    if (verbose) {
      THREAD_log(
        "Response",
        "Merged GWAS Beta/SE into matches: ",
        sum(!is.na(matches$Beta) & !is.na(matches$SE)),
        " matched SNP-gene row(s)."
      )
    }
  } else {
    if (!"Beta" %in% names(matches)) {
      if ("beta" %in% names(matches)) {
        matches[, Beta := suppressWarnings(as.numeric(beta))]
      } else if ("BETA" %in% names(matches)) {
        matches[, Beta := suppressWarnings(as.numeric(BETA))]
      } else {
        stop(
          "No 'Beta' column was found in matches and no 'gwas_clean' ",
          "object was provided."
        )
      }
    }

    if (!"SE" %in% names(matches)) {
      if ("se" %in% names(matches)) {
        matches[, SE := suppressWarnings(as.numeric(se))]
      } else {
        stop(
          "No 'SE' column was found in matches and no 'gwas_clean' ",
          "object was provided."
        )
      }
    }
  }

  matches[, Beta := suppressWarnings(as.numeric(Beta))]
  matches[, SE := suppressWarnings(as.numeric(SE))]

  matches <- matches[
    is.finite(Beta) & is.finite(SE) & SE > 0 &
      !is.na(gene_name) & gene_name != "" &
      !is.na(chrom) & chrom != "" &
      !is.na(key_value) & key_value != ""
  ]

  matches
}

standardize_ld_table_for_response <- function(ld_chr,
                                              match_mode = c(
                                                "snp_name",
                                                "snp_pos"
                                              ),
                                              ld_snp_a_col = "SNP_A",
                                              ld_snp_b_col = "SNP_B",
                                              ld_bp_a_col = "BP_A",
                                              ld_bp_b_col = "BP_B",
                                              ld_r_col = "R",
                                              verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  match_mode <- match.arg(match_mode)
  ld_chr <- data.table::as.data.table(ld_chr)

  if (match_mode == "snp_name") {
    validate_compute_response_columns(
      ld_chr,
      required = c(ld_snp_a_col, ld_snp_b_col, ld_r_col),
      object_name = "LD table"
    )

    out <- data.table::data.table(
      A = as.character(ld_chr[[ld_snp_a_col]]),
      B = as.character(ld_chr[[ld_snp_b_col]]),
      R = suppressWarnings(as.numeric(ld_chr[[ld_r_col]]))
    )
  } else {
    validate_compute_response_columns(
      ld_chr,
      required = c(ld_bp_a_col, ld_bp_b_col, ld_r_col),
      object_name = "LD table"
    )

    out <- data.table::data.table(
      A = normalize_position_key(ld_chr[[ld_bp_a_col]]),
      B = normalize_position_key(ld_chr[[ld_bp_b_col]]),
      R = suppressWarnings(as.numeric(ld_chr[[ld_r_col]]))
    )
  }

  out <- out[
    !is.na(A) & A != "" &
      !is.na(B) & B != "" &
      is.finite(R)
  ]

  bad_r <- out[abs(R) > 1 + 1e-8, .N]

  if (bad_r > 0L) {
    stop(
      "LD table contains ", bad_r,
      " record(s) with |R| > 1."
    )
  }

  # Duplicated directed pairs are uncommon but can occur after merging LD
  # sources. If present, aggregate them deterministically.
  if (anyDuplicated(out, by = c("A", "B")) > 0L) {
    if (verbose) {
      warning(
        "Duplicated directed LD pair(s) were detected and averaged."
      )
    }

    out <- out[
      ,
      .(R = mean(R)),
      by = .(A, B)
    ]
  }

  data.table::setkey(out, A, B)
  out
}

construct_complete_symmetric_ld <- function(g_keys,
                                            ld_chr) {
  k <- length(g_keys)

  if (k == 1L) {
    return(
      list(
        R_mat = matrix(1, nrow = 1L, ncol = 1L),
        forward_coverage = 1,
        reverse_coverage = 1,
        both_present_fraction = 1,
        missing_both_fraction = 0,
        max_reciprocal_difference = 0,
        max_asymmetry = 0
      )
    )
  }

  pairs <- data.table::CJ(
    A = g_keys,
    B = g_keys,
    sorted = TRUE
  )

  pairs[, pair_id := .I]

  forward <- ld_chr[
    pairs,
    on = .(A, B),
    .(
      pair_id = i.pair_id,
      R_forward = x.R
    )
  ]

  reverse_query <- pairs[
    ,
    .(
      pair_id,
      A = B,
      B = A
    )
  ]

  reverse <- ld_chr[
    reverse_query,
    on = .(A, B),
    .(
      pair_id = i.pair_id,
      R_reverse = x.R
    )
  ]

  data.table::setorder(forward, pair_id)
  data.table::setorder(reverse, pair_id)

  if (nrow(forward) != nrow(pairs) ||
    nrow(reverse) != nrow(pairs)) {
    stop("Internal LD reciprocal lookup length mismatch.")
  }

  r_forward <- forward$R_forward
  r_reverse <- reverse$R_reverse

  diagonal_idx <- pairs$A == pairs$B
  off_diagonal_idx <- !diagonal_idx

  both_present <- is.finite(r_forward) & is.finite(r_reverse)
  forward_only <- is.finite(r_forward) & !is.finite(r_reverse)
  reverse_only <- !is.finite(r_forward) & is.finite(r_reverse)
  missing_both <- !is.finite(r_forward) & !is.finite(r_reverse)

  # If both directed records exist, use their mean. Under a valid LD file the
  # two values should be equal up to numerical precision.
  r_complete <- numeric(length(r_forward))

  r_complete[both_present] <-
    (r_forward[both_present] + r_reverse[both_present]) / 2

  r_complete[forward_only] <- r_forward[forward_only]
  r_complete[reverse_only] <- r_reverse[reverse_only]
  r_complete[missing_both] <- 0
  r_complete[diagonal_idx] <- 1

  R_mat <- matrix(
    r_complete,
    nrow = k,
    ncol = k,
    byrow = TRUE,
    dimnames = list(g_keys, g_keys)
  )

  # This final averaging is deterministic protection against any residual
  # floating-point or duplicated-record disagreement.
  R_mat <- (R_mat + t(R_mat)) / 2
  diag(R_mat) <- 1

  reciprocal_difference <- abs(r_forward - r_reverse)
  reciprocal_difference <- reciprocal_difference[
    off_diagonal_idx & is.finite(reciprocal_difference)
  ]

  max_reciprocal_difference <- if (
    length(reciprocal_difference) > 0L
  ) {
    max(reciprocal_difference)
  } else {
    NA_real_
  }

  list(
    R_mat = R_mat,
    forward_coverage = mean(
      is.finite(r_forward[off_diagonal_idx])
    ),
    reverse_coverage = mean(
      is.finite(r_reverse[off_diagonal_idx])
    ),
    both_present_fraction = mean(
      both_present[off_diagonal_idx]
    ),
    missing_both_fraction = mean(
      missing_both[off_diagonal_idx]
    ),
    max_reciprocal_difference = max_reciprocal_difference,
    max_asymmetry = max(abs(R_mat - t(R_mat)))
  )
}

# Compute trace powers after a trace-preserving PSD spectral repair.
#
# For a valid LD correlation matrix, all eigenvalues are nonnegative and sum to
# k. Sparse pairwise storage can produce an indefinite symmetric matrix after
# missing pairs are filled with zero. We therefore:
#   1. clip negative eigenvalues to zero;
#   2. rescale the retained eigenvalues so that their sum is k.
#
# This ensures the effective matrix remains PSD and satisfies tr(R) = k, which
# is required by E[T_j] = theta_j tr(R^2) + k v_j.
compute_ld_trace_powers <- function(R_mat,
                                    psd_tolerance = 1e-8) {
  k <- nrow(R_mat)

  eigen_raw <- eigen(
    R_mat,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  min_eigenvalue_raw <- min(eigen_raw)
  n_negative_eigenvalues <- sum(eigen_raw < -psd_tolerance)

  negative_eigen_mass <- sum(abs(eigen_raw[eigen_raw < 0]))
  total_eigen_mass <- sum(abs(eigen_raw))

  negative_eigen_mass_ratio <- if (total_eigen_mass > 0) {
    negative_eigen_mass / total_eigen_mass
  } else {
    NA_real_
  }

  eigen_used <- pmax(eigen_raw, 0)
  eigen_sum <- sum(eigen_used)

  if (!is.finite(eigen_sum) || eigen_sum <= 0) {
    stop("LD eigenvalue repair failed: non-positive retained eigenvalue sum.")
  }

  # Preserve tr(R) = k after clipping.
  eigen_used <- eigen_used * (k / eigen_sum)

  list(
    tr_R2 = sum(eigen_used^2),
    tr_R3 = sum(eigen_used^3),
    tr_R4 = sum(eigen_used^4),
    min_eigenvalue_raw = min_eigenvalue_raw,
    n_negative_eigenvalues = n_negative_eigenvalues,
    negative_eigen_mass_ratio = negative_eigen_mass_ratio,
    psd_repaired = n_negative_eigenvalues > 0L
  )
}

compute_gene_hetero_metrics <- function(g_data,
                                        ld_chr,
                                        max_snps_per_gene = Inf,
                                        psd_tolerance = 1e-8,
                                        include_diagnostics = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  g_data <- data.table::as.data.table(g_data)

  g_unique <- g_data[
    !duplicated(key_value)
  ][order(key_value)]

  g_name <- g_unique$gene_name[1L]
  g_keys <- g_unique$key_value
  g_betas <- g_unique$Beta
  g_ses <- g_unique$SE
  k <- length(g_keys)

  empty_diagnostics <- list(
    ld_forward_coverage = NA_real_,
    ld_reverse_coverage = NA_real_,
    ld_both_present_fraction = NA_real_,
    ld_missing_both_fraction = NA_real_,
    ld_max_reciprocal_difference = NA_real_,
    ld_max_asymmetry = NA_real_,
    min_eigenvalue_raw = NA_real_,
    n_negative_eigenvalues = NA_integer_,
    negative_eigen_mass_ratio = NA_real_,
    psd_repaired = NA
  )

  make_result <- function(d0 = NA_real_,
                          d1 = NA_real_,
                          d2 = NA_real_,
                          final_y = NA_real_,
                          status,
                          diagnostics = empty_diagnostics) {
    out <- data.table::data.table(
      gene_name = g_name,
      n_snps = k,
      d0 = d0,
      d1 = d1,
      d2 = d2,
      final_y = final_y,
      status = status
    )

    if (include_diagnostics) {
      out[, `:=`(
        ld_forward_coverage = diagnostics$ld_forward_coverage,
        ld_reverse_coverage = diagnostics$ld_reverse_coverage,
        ld_both_present_fraction = diagnostics$ld_both_present_fraction,
        ld_missing_both_fraction = diagnostics$ld_missing_both_fraction,
        ld_max_reciprocal_difference =
          diagnostics$ld_max_reciprocal_difference,
        ld_max_asymmetry = diagnostics$ld_max_asymmetry,
        min_eigenvalue_raw = diagnostics$min_eigenvalue_raw,
        n_negative_eigenvalues = diagnostics$n_negative_eigenvalues,
        negative_eigen_mass_ratio =
          diagnostics$negative_eigen_mass_ratio,
        psd_repaired = diagnostics$psd_repaired
      )]
    }

    out
  }

  if (k == 0L) {
    return(NULL)
  }

  if (is.finite(max_snps_per_gene) &&
    k > max_snps_per_gene) {
    return(
      make_result(status = "skipped_too_many_snps")
    )
  }

  if (any(!is.finite(g_betas)) ||
    any(!is.finite(g_ses)) ||
    any(g_ses <= 0)) {
    return(
      make_result(status = "non_finite_beta_or_se")
    )
  }

  median_se_sq <- stats::median(g_ses^2)
  sum_beta2 <- sum(g_betas^2)

  if (!is.finite(median_se_sq) ||
    !is.finite(sum_beta2)) {
    return(
      make_result(status = "non_finite_beta_or_se")
    )
  }

  ld_info <- tryCatch(
    construct_complete_symmetric_ld(
      g_keys = g_keys,
      ld_chr = ld_chr
    ),
    error = function(e) e
  )

  if (inherits(ld_info, "error")) {
    return(
      make_result(status = "ld_matrix_construction_failed")
    )
  }

  trace_info <- tryCatch(
    compute_ld_trace_powers(
      R_mat = ld_info$R_mat,
      psd_tolerance = psd_tolerance
    ),
    error = function(e) e
  )

  if (inherits(trace_info, "error")) {
    return(
      make_result(status = "invalid_ld_eigendecomposition")
    )
  }

  diagnostics <- list(
    ld_forward_coverage = ld_info$forward_coverage,
    ld_reverse_coverage = ld_info$reverse_coverage,
    ld_both_present_fraction = ld_info$both_present_fraction,
    ld_missing_both_fraction = ld_info$missing_both_fraction,
    ld_max_reciprocal_difference =
      ld_info$max_reciprocal_difference,
    ld_max_asymmetry = ld_info$max_asymmetry,
    min_eigenvalue_raw = trace_info$min_eigenvalue_raw,
    n_negative_eigenvalues = trace_info$n_negative_eigenvalues,
    negative_eigen_mass_ratio =
      trace_info$negative_eigen_mass_ratio,
    psd_repaired = trace_info$psd_repaired
  )

  tr_R2 <- trace_info$tr_R2
  tr_R3 <- trace_info$tr_R3
  tr_R4 <- trace_info$tr_R4

  if (!is.finite(tr_R2) || tr_R2 <= 0 ||
    !is.finite(tr_R3) ||
    !is.finite(tr_R4)) {
    return(
      make_result(
        status = "invalid_trace",
        diagnostics = diagnostics
      )
    )
  }

  d0_val <- (2 * median_se_sq^2) / tr_R2
  d1_val <- (4 * median_se_sq * tr_R3) / (tr_R2^2)
  d2_val <- (2 * tr_R4) / (tr_R2^2)

  y_val <- (
    sum_beta2 - k * median_se_sq
  ) / tr_R2

  make_result(
    d0 = d0_val,
    d1 = d1_val,
    d2 = d2_val,
    final_y = y_val,
    status = "ok",
    diagnostics = diagnostics
  )
}

#' Compute LD-based heteroscedastic gene-level THREAD response
#'
#' @param gene_snp_matches A data.frame/data.table or file path to gene-SNP
#'   matches.
#' @param matches Deprecated alias for gene_snp_matches.
#' @param gwas_clean Optional data.frame/data.table or file path containing beta
#'   and SE.
#' @param ld_dir Directory containing chromosome-level LD RDS files.
#' @param match_mode Variant-key mode. Use "snp_name" (default and recommended)
#'   for rsID-based matching, or "snp_pos" for chromosome-position matching.
#'   The selected mode controls both GWAS matching and LD lookup.
#' @param genes Optional character vector restricting computation to selected
#'   genes.
#' @param chromosomes Optional chromosomes to process.
#' @param n_cores Number of forked workers for gene-level parallelism.
#' @param set_dt_threads Logical. If TRUE, data.table threads are set to 1 during
#'   computation.
#' @param max_genes_per_chr Optional maximum number of genes per chromosome for
#'   smoke tests.
#' @param max_snps_per_gene Maximum SNPs per gene. Larger genes are skipped.
#' @param psd_tolerance Eigenvalues below -psd_tolerance are counted as negative.
#' @param include_diagnostics Logical. Include LD completeness, symmetry and PSD
#'   diagnostic columns in the returned table.
#' @param output_file Optional final merged RDS output path.
#' @param output_dir Optional directory for chromosome-level checkpoint files.
#' @param save_by_chr Logical. Whether to save one response file per chromosome.
#' @param resume Logical. Whether to reuse existing chromosome-level files when
#'   save_by_chr is TRUE.
#' @param snp_col SNP-ID column in matches.
#' @param snp_pos_col SNP-position column in matches.
#' @param gene_col Gene column in matches.
#' @param chr_col Chromosome column in matches.
#' @param gwas_snp_col SNP-ID column in gwas_clean.
#' @param gwas_pos_col Position column in gwas_clean. If NULL, inferred for
#'   match_mode = "snp_pos".
#' @param gwas_chr_col Chromosome column in gwas_clean. Required for position
#'   matching.
#' @param beta_col Beta column in gwas_clean.
#' @param se_col Standard-error column in gwas_clean.
#' @param ld_file_template Character string defining chromosome-specific LD
#'   files under \code{ld_dir}. The token \code{chr}, enclosed in braces, is
#'   replaced by the chromosome value. The default corresponds to one RDS file
#'   named by chromosome.
#' @param ld_snp_a_col SNP_A column in SNP-name LD files.
#' @param ld_snp_b_col SNP_B column in SNP-name LD files.
#' @param ld_bp_a_col BP_A column in position LD files.
#' @param ld_bp_b_col BP_B column in position LD files.
#' @param ld_r_col LD correlation column.
#' @param drop_failed Logical. Whether to remove genes with non-ok status.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A data.frame containing gene_name, n_snps, d0, d1, d2, final_y,
#'   status and, optionally, LD diagnostic columns.
#'
#' @export
compute_response <- function(gene_snp_matches = NULL,
                             matches = NULL,
                             gwas_clean = NULL,
                             ld_dir,
                             match_mode = c("snp_name", "snp_pos"),
                             genes = NULL,
                             chromosomes = NULL,
                             n_cores = 1L,
                             set_dt_threads = TRUE,
                             max_genes_per_chr = NULL,
                             max_snps_per_gene = Inf,
                             psd_tolerance = 1e-8,
                             include_diagnostics = TRUE,
                             output_file = NULL,
                             output_dir = NULL,
                             save_by_chr = FALSE,
                             resume = FALSE,
                             snp_col = "snp_name",
                             snp_pos_col = "snp_pos",
                             gene_col = "gene_name",
                             chr_col = "chrom",
                             gwas_snp_col = "rsid",
                             gwas_pos_col = NULL,
                             gwas_chr_col = "chrom",
                             beta_col = "beta",
                             se_col = "se",
                             ld_file_template = "{chr}.rds",
                             ld_snp_a_col = "SNP_A",
                             ld_snp_b_col = "SNP_B",
                             ld_bp_a_col = "BP_A",
                             ld_bp_b_col = "BP_B",
                             ld_r_col = "R",
                             drop_failed = FALSE,
                             verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  match_mode <- match.arg(match_mode)

  if (is.null(gene_snp_matches)) {
    gene_snp_matches <- matches
  }

  if (is.null(gene_snp_matches)) {
    stop("'gene_snp_matches' must be provided.")
  }

  if (missing(ld_dir) ||
    !is.character(ld_dir) ||
    length(ld_dir) != 1L ||
    is.na(ld_dir)) {
    stop("'ld_dir' must be a single directory path.")
  }

  if (!dir.exists(ld_dir)) {
    stop("'ld_dir' does not exist: ", ld_dir)
  }

  if (!is.numeric(n_cores) ||
    length(n_cores) != 1L ||
    is.na(n_cores) ||
    n_cores < 1) {
    stop("'n_cores' must be a positive integer.")
  }

  n_cores <- as.integer(n_cores)

  if (n_cores > 1L && .Platform$OS.type == "windows") {
    warning(
      "Forked parallelism is not available on Windows; using one core."
    )
    n_cores <- 1L
  }

  old_dt_threads <- NULL

  if (set_dt_threads) {
    old_dt_threads <- data.table::getDTthreads()
    data.table::setDTthreads(1L)

    on.exit(
      data.table::setDTthreads(old_dt_threads),
      add = TRUE
    )
  }

  if (!is.null(max_genes_per_chr)) {
    if (!is.numeric(max_genes_per_chr) ||
      length(max_genes_per_chr) != 1L ||
      is.na(max_genes_per_chr) ||
      max_genes_per_chr < 1) {
      stop(
        "'max_genes_per_chr' must be NULL or a positive integer."
      )
    }

    max_genes_per_chr <- as.integer(max_genes_per_chr)
  }

  if (!is.numeric(max_snps_per_gene) ||
    length(max_snps_per_gene) != 1L ||
    is.na(max_snps_per_gene) ||
    max_snps_per_gene < 1) {
    stop("'max_snps_per_gene' must be a positive number or Inf.")
  }

  if (!is.numeric(psd_tolerance) ||
    length(psd_tolerance) != 1L ||
    is.na(psd_tolerance) ||
    psd_tolerance < 0) {
    stop("'psd_tolerance' must be a nonnegative number.")
  }

  if (save_by_chr && is.null(output_dir)) {
    if (!is.null(output_file)) {
      output_dir <- dirname(output_file)
    } else {
      stop(
        "'output_dir' must be provided when save_by_chr = TRUE ",
        "and output_file is NULL."
      )
    }
  }

  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(
      output_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }

  if (verbose) {
    THREAD_log(
      "Response",
      sprintf(
        "compute_response settings: match_mode = %s, LD symmetrization = always enabled, PSD repair = trace-preserving eigenvalue clipping, n_cores = %d, set_dt_threads = %s.",
        match_mode,
        n_cores,
        if (set_dt_threads) "TRUE" else "FALSE"
      )
    )
  }

  gene_snp_matches <- load_table_or_object(
    gene_snp_matches,
    object_name = "gene_snp_matches"
  )

  matches_dt <- standardize_matches_for_response(
    matches = gene_snp_matches,
    match_mode = match_mode,
    snp_col = snp_col,
    snp_pos_col = snp_pos_col,
    gene_col = gene_col,
    chr_col = chr_col
  )

  if (!is.null(genes)) {
    genes <- as.character(genes)
    before <- data.table::uniqueN(matches_dt$gene_name)
    matches_dt <- matches_dt[gene_name %in% genes]

    if (verbose) {
      THREAD_log(
        "Response",
        sprintf(
          "Restricting gene-SNP matches to requested genes: %d of %d gene(s).",
          data.table::uniqueN(matches_dt$gene_name),
          before
        )
      )
    }

    if (nrow(matches_dt) == 0L) {
      stop("None of 'genes' were found in gene-SNP matches.")
    }
  }

  matches_dt <- merge_gwas_effects_for_response(
    matches = matches_dt,
    gwas_clean = gwas_clean,
    match_mode = match_mode,
    gwas_snp_col = gwas_snp_col,
    gwas_pos_col = gwas_pos_col,
    gwas_chr_col = gwas_chr_col,
    beta_col = beta_col,
    se_col = se_col,
    verbose = verbose
  )

  if (nrow(matches_dt) == 0L) {
    stop("No valid SNP-gene rows remained after merging Beta/SE.")
  }

  if (is.null(chromosomes)) {
    chromosomes <- sort(unique(matches_dt$chrom))
  } else {
    chromosomes <- normalize_chr_response(chromosomes)
  }

  all_results <- list()

  if (verbose) {
    THREAD_log("Response", sprintf(
      "Computing THREAD heteroscedastic response over %d chromosome(s).",
      length(chromosomes)
    ))
  }

  for (chr in chromosomes) {
    if (verbose) {
      THREAD_log("Response", sprintf("Processing chromosome %s ...", chr))
    }

    chr_file <- NULL

    if (save_by_chr) {
      chr_file <- file.path(
        output_dir,
        paste0("chr", chr, "_response.rds")
      )

      if (resume && file.exists(chr_file)) {
        if (verbose) {
          THREAD_log("Response", "Using existing chromosome checkpoint: ", chr_file)
        }

        all_results[[chr]] <- data.table::as.data.table(
          readRDS(chr_file)
        )

        next
      }
    }

    ld_path <- make_ld_file_path(
      ld_dir,
      chr,
      ld_file_template
    )

    if (!file.exists(ld_path)) {
      warning(
        "LD file not found for chromosome ",
        chr,
        ": ",
        ld_path
      )

      next
    }

    ld_chr <- readRDS(ld_path)

    ld_chr <- standardize_ld_table_for_response(
      ld_chr = ld_chr,
      match_mode = match_mode,
      ld_snp_a_col = ld_snp_a_col,
      ld_snp_b_col = ld_snp_b_col,
      ld_bp_a_col = ld_bp_a_col,
      ld_bp_b_col = ld_bp_b_col,
      ld_r_col = ld_r_col,
      verbose = verbose
    )

    chr_matches <- matches_dt[chrom == chr]

    if (nrow(chr_matches) == 0L) {
      if (verbose) {
        THREAD_log("Response", sprintf("No matches found for chromosome %s.", chr))
      }

      rm(ld_chr)
      gc(verbose = FALSE)
      next
    }

    genes_list <- split(
      chr_matches,
      by = "gene_name",
      keep.by = TRUE
    )

    if (!is.null(max_genes_per_chr) &&
      length(genes_list) > max_genes_per_chr) {
      genes_list <- genes_list[seq_len(max_genes_per_chr)]

      if (verbose) {
        THREAD_log(
          "Response",
          "Restricted to ",
          max_genes_per_chr,
          " gene(s) for chromosome ",
          chr,
          "."
        )
      }
    }

    worker_fun <- function(g_data) {
      compute_gene_hetero_metrics(
        g_data = g_data,
        ld_chr = ld_chr,
        max_snps_per_gene = max_snps_per_gene,
        psd_tolerance = psd_tolerance,
        include_diagnostics = include_diagnostics
      )
    }

    if (n_cores > 1L &&
      .Platform$OS.type != "windows") {
      chr_res_list <- parallel::mclapply(
        genes_list,
        worker_fun,
        mc.cores = n_cores
      )
    } else {
      chr_res_list <- lapply(
        genes_list,
        worker_fun
      )
    }

    chr_res_list <- Filter(
      Negate(is.null),
      chr_res_list
    )

    chr_res <- data.table::rbindlist(
      chr_res_list,
      fill = TRUE
    )

    if (nrow(chr_res) > 0L &&
      "status" %in% names(chr_res)) {
      if (verbose) {
        THREAD_log("Response", sprintf("Status summary for chromosome %s:", chr))
        print(table(chr_res$status))
      }
    }

    if (verbose &&
      include_diagnostics &&
      nrow(chr_res) > 0L) {
    ok_rows <- chr_res[status == "ok"]

    if (nrow(ok_rows) > 0L) {
      THREAD_log(
        "Response",
        message(
          sprintf(
            paste0(
              "[THREAD:Response] LD diagnostics for chromosome %s: ",
              "median missing-both fraction = %.4f; ",
              "median negative-eigen mass ratio = %.4f."
            ),
            chr,
            stats::median(
              ok_rows$ld_missing_both_fraction,
              na.rm = TRUE
            ),
            stats::median(
              ok_rows$negative_eigen_mass_ratio,
              na.rm = TRUE
            )
          )
        )
      )
    }
    }

    if (verbose) {
      n_ok <- if ("status" %in% names(chr_res)) {
        sum(chr_res$status == "ok")
      } else {
        nrow(chr_res)
      }

      THREAD_log(
        "Response",
        sprintf(
          "Chromosome %s complete: %d successful gene(s), %d returned record(s).",
          chr, n_ok, nrow(chr_res)
        )
      )
    }

    if (save_by_chr) {
      saveRDS(chr_res, chr_file)

      if (verbose) {
        THREAD_log("Response", "Saved chromosome response to: ", chr_file)
      }
    }

    all_results[[chr]] <- chr_res

    rm(
      ld_chr,
      chr_matches,
      genes_list,
      chr_res_list,
      chr_res
    )

    gc(verbose = FALSE)
  }

  if (length(all_results) == 0L) {
    stop("No chromosome produced response results.")
  }

  result <- data.table::rbindlist(
    all_results,
    fill = TRUE
  )

  if (drop_failed && "status" %in% names(result)) {
    result <- result[status == "ok"]
  }

  data.table::setorder(result, gene_name)

  required_output <- c(
    "gene_name",
    "n_snps",
    "d0",
    "d1",
    "d2",
    "final_y"
  )

  validate_compute_response_columns(
    result,
    required_output,
    "computed response"
  )

  if (verbose) {
    THREAD_log(
      "Response",
      sprintf(
        "Computed response for %d gene(s). Valid finite final_y: %d.",
        nrow(result),
        sum(is.finite(result$final_y))
      )
    )
  }

  if (!is.null(output_file)) {
    output_dir_final <- dirname(output_file)

    if (!dir.exists(output_dir_final)) {
      dir.create(
        output_dir_final,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }

    if (!dir.exists(output_dir_final)) {
      stop(
        "Failed to create output directory: ",
        output_dir_final
      )
    }

    saveRDS(result, file = output_file)

    if (verbose) {
      THREAD_log("Response", "Saved computed response to: ", output_file)
    }
  }

  as.data.frame(result)
}

#' Combine chromosome-level response files
#'
#' @param output_dir Directory containing chr*_response.rds files.
#' @param pattern File pattern for chromosome-level response files.
#' @param output_file Optional merged output RDS path.
#'
#' @return A data.frame containing the merged response table.
#'
#' @export
combine_response_by_chr <- function(output_dir,
                                    pattern = "^chr.*_response\\.rds$",
                                    output_file = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  if (!dir.exists(output_dir)) {
    stop("'output_dir' does not exist: ", output_dir)
  }

  files <- list.files(
    output_dir,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0L) {
    stop("No chromosome response files found in: ", output_dir)
  }

  res <- lapply(files, readRDS)
  res <- data.table::rbindlist(res, fill = TRUE)
  data.table::setorder(res, gene_name)

  if (!is.null(output_file)) {
    dir.create(
      dirname(output_file),
      recursive = TRUE,
      showWarnings = FALSE
    )

    saveRDS(res, output_file)
  }

  as.data.frame(res)
}
