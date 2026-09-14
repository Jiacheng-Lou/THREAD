# gwas.R
#   preprocess_gwas(): raw GWAS QC and column harmonisation

# gene_snp_mapping.R
#   match_snps_to_genes(): SNP-to-gene window mapping

# compute_response.R
#   compute_response(): LD-corrected gene response calculation
#
# Each step is independent: a user who already holds a clean sumstats table, a
# SNP-gene map, or pre-computed responses can enter the chain at any point.
# -----------------------------------------------------------------------------

# Internal helper: strip a leading 'chr' (any case) and coerce to integer.
normalize_chrom <- function(x) {
  suppressWarnings(as.integer(gsub("^chr", "", as.character(x), ignore.case = TRUE)))
}


#' Quality-control and harmonise GWAS summary statistics
#'
#' @description
#' Applies a standard set of GWAS quality-control filters and returns a clean,
#' uniformly named table that the downstream THREAD functions consume. The input
#' may use any column names; supply \code{col_map} to declare which of your
#' columns correspond to the required fields. A user who already holds a clean
#' table with the standard column names can pass it straight to
#' \code{\link{match_snps_to_genes}} and \code{\link{compute_response}} and skip
#' this step.
#'
#' @param sumstats A data.frame / data.table of summary statistics, or a single
#'   string giving the path to a delimited file (read with
#'   \code{data.table::fread}).
#' @param col_map Optional named character vector mapping standard field names
#'   to the column names used in \code{sumstats}. Recognised standard names are
#'   \code{chrom}, \code{pos}, \code{rsid}, \code{beta}, \code{se} (required)
#'   and \code{maf}, \code{freq}, \code{info}, \code{pval}, \code{n} (optional;
#'   each enables its corresponding filter). Example:
#'   \code{c(chrom = "CHR", pos = "POS", rsid = "SNP", beta = "BETA",
#'   se = "SE")}. If \code{NULL}, columns are assumed to already use the
#'   standard names.
#' @param remove_sex_chrom Logical; keep only autosomes 1-22. Default \code{TRUE}.
#' @param remove_mhc Logical; drop the MHC region on chromosome 6. Default
#'   \code{TRUE}.
#' @param mhc_range Integer length-2 giving the MHC start/end in base pairs
#'   (build hg19/GRCh37). Default \code{c(25e6, 34e6)}.
#' @param maf_min Minimum minor-allele frequency; applied only when a \code{maf}
#'   or \code{freq} column is available. Default \code{0.01}.
#' @param info_min Minimum imputation INFO score; applied only when an
#'   \code{info} column is available. Default \code{0.8}.
#' @param n_min Minimum per-SNP sample size; applied only when an \code{n}
#'   column is available and \code{n_min > 0}. Default \code{0} (no filter).
#' @param beta_max Maximum absolute effect size; SNPs with \code{abs(beta)}
#'   above this are dropped as likely artefacts. Default \code{10}. Set to
#'   \code{Inf} to disable.
#' @param verbose Logical; print a per-step row-count report. Default
#'   \code{TRUE}.
#'
#' @return A \code{data.table} with columns \code{chrom}, \code{pos},
#'   \code{rsid}, \code{beta}, \code{se}, and \code{maf} when available, sorted
#'   by \code{chrom} then \code{pos}.
#'
#' @examples
#' \dontrun{
#' clean <- preprocess_gwas(
#'   "PGC_BIP.tsv",
#'   col_map = c(
#'     chrom = "CHROM", pos = "POS", rsid = "ID",
#'     beta = "BETA", se = "SE", info = "IMPINFO"
#'   )
#' )
#' }
#'
#' @importFrom data.table fread as.data.table setnames setorder :=
#' @export
preprocess_gwas <- function(sumstats,
                            col_map = NULL,
                            remove_sex_chrom = TRUE,
                            remove_mhc = TRUE,
                            mhc_range = c(25e6, 34e6),
                            maf_min = 0.01,
                            info_min = 0.8,
                            n_min = 0,
                            beta_max = 10,
                            verbose = TRUE) {
  # ---- load ------------------------------------------------------------------
  if (is.character(sumstats) && length(sumstats) == 1L) {
    if (!file.exists(sumstats)) stop("File not found: ", sumstats)
    if (verbose) THREAD_log("GWAS", "Reading summary statistics from:", sumstats)
    dt <- data.table::fread(sumstats, na.strings = c("NA", "", "."))
  } else {
    dt <- data.table::as.data.table(sumstats)
  }
  if (nrow(dt) == 0L) stop("Input summary statistics are empty.")

  # ---- harmonise column names ------------------------------------------------
  if (!is.null(col_map)) {
    if (is.null(names(col_map))) {
      stop("'col_map' must be a *named* character vector, e.g. c(chrom = 'CHR').")
    }
    for (std in names(col_map)) {
      user_col <- col_map[[std]]
      if (!user_col %in% names(dt)) {
        stop(sprintf(
          "Column '%s' (mapped to '%s') not found in the input.",
          user_col, std
        ))
      }
      if (user_col != std) data.table::setnames(dt, user_col, std)
    }
  }

  required <- c("chrom", "pos", "rsid", "beta", "se")
  missing_req <- setdiff(required, names(dt))
  if (length(missing_req) > 0L) {
    stop(
      "Missing required column(s): ", paste(missing_req, collapse = ", "),
      ". Provide them via 'col_map'."
    )
  }

  # ---- derive maf from an allele-frequency column when needed ----------------
  if (!"maf" %in% names(dt) && "freq" %in% names(dt)) {
    dt[, maf := pmin(
      suppressWarnings(as.numeric(freq)),
      1 - suppressWarnings(as.numeric(freq))
    )]
  }

  # ---- type coercion ---------------------------------------------------------
  dt[, chrom := normalize_chrom(chrom)]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt[, rsid := as.character(rsid)]
  dt[, beta := suppressWarnings(as.numeric(beta))]
  dt[, se := suppressWarnings(as.numeric(se))]
  if ("maf" %in% names(dt)) dt[, maf := suppressWarnings(as.numeric(maf))]
  if ("info" %in% names(dt)) dt[, info := suppressWarnings(as.numeric(info))]
  if ("n" %in% names(dt)) dt[, n := suppressWarnings(as.numeric(n))]
  if ("pval" %in% names(dt)) dt[, pval := suppressWarnings(as.numeric(pval))]

  if (all(is.na(dt$chrom))) {
    stop(
      "All chromosome values are NA after coercion. ",
      "Check the chromosome coding (the 'chr' prefix is stripped ",
      "automatically, but values such as 'X'/'Y'/'MT' become NA)."
    )
  }

  n0 <- nrow(dt)
  if (verbose) THREAD_log("GWAS", sprintf("Loaded %d SNPs.", n0))

  log_step <- function(before, after, label) {
    if (verbose) {
      THREAD_log("GWAS", sprintf(
        "%-32s %9d -> %9d  (removed %d)",
        label, before, after, before - after
      ))
    }
  }

  # ---- filters ---------------------------------------------------------------
  b <- nrow(dt)
  dt <- dt[!is.na(chrom) & !is.na(pos) & !is.na(rsid) & !is.na(beta) & !is.na(se)]
  log_step(b, nrow(dt), "drop NA in required fields")

  if (remove_sex_chrom) {
    b <- nrow(dt)
    dt <- dt[chrom %in% 1:22]
    log_step(b, nrow(dt), "keep autosomes 1-22")
  }

  if (remove_mhc) {
    b <- nrow(dt)
    dt <- dt[!(chrom == 6 & pos >= mhc_range[1] & pos <= mhc_range[2])]
    log_step(b, nrow(dt), "remove MHC region (chr6)")
  }

  b <- nrow(dt)
  dt <- dt[se > 0]
  log_step(b, nrow(dt), "drop SE <= 0")

  if (is.finite(beta_max)) {
    b <- nrow(dt)
    dt <- dt[abs(beta) <= beta_max]
    log_step(b, nrow(dt), sprintf("drop |beta| > %g", beta_max))
  }

  if ("maf" %in% names(dt)) {
    b <- nrow(dt)
    dt <- dt[!is.na(maf) & maf >= maf_min]
    log_step(b, nrow(dt), sprintf("drop MAF < %g", maf_min))
  } else if (verbose) {
    THREAD_log("GWAS", "(no MAF/freq column - MAF filter skipped)")
  }

  if ("info" %in% names(dt)) {
    b <- nrow(dt)
    dt <- dt[!is.na(info) & info >= info_min]
    log_step(b, nrow(dt), sprintf("drop INFO < %g", info_min))
  }

  if ("n" %in% names(dt) && n_min > 0) {
    b <- nrow(dt)
    dt <- dt[!is.na(n) & n >= n_min]
    log_step(b, nrow(dt), sprintf("drop N < %g", n_min))
  }

  if ("pval" %in% names(dt)) {
    b <- nrow(dt)
    dt <- dt[is.finite(pval) & pval > 0 & pval <= 1]
    log_step(b, nrow(dt), "drop invalid P-values")
  }

  # ---- select, sort, return --------------------------------------------------
  out_cols <- c("chrom", "pos", "rsid", "beta", "se")
  if ("maf" %in% names(dt)) out_cols <- c(out_cols, "maf")
  dt <- dt[, ..out_cols]
  data.table::setorder(dt, chrom, pos)

  if (nrow(dt) == 0L) warning("No SNPs passed quality control.")
  if (verbose) THREAD_log("GWAS", sprintf("Done. %d of %d SNPs passed QC.", nrow(dt), n0))
  dt[]
}
