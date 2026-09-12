#' THREAD: Trait Heterogeneity Regression through Expression Annotation and Dirichlet processes
#'
#' THREAD integrates gene-level GWAS response statistics with gene-by-label
#' single-cell expression profiles to infer latent gene modules and
#' module-specific cellular coefficients.
#'
#' @keywords internal
#' @useDynLib THREAD, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

# Declare data.table columns referenced by non-standard evaluation so that
# R CMD check does not raise "no visible binding for global variable" notes.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    ".", "..out_cols", ".I", ".N",
    "A", "ambiguous_key", "B", "Beta", "beta", "BETA",
    "chrom", "end", "freq", "gene_name",
    "i.beta", "i.pair_id", "i.se", "info",
    "key_value", "maf", "n", "N", "pair_id", "pos", "pval",
    "R", "rsid", "se", "SE", "SNP", "SNP_A", "SNP_B", "snp_name", "snp_pos",
    "start", "status", "x.R"
  ))
}

# Internal helper for consistent package logging.
DPM_log <- function(module, ..., sep = " ") {
  message(sprintf("[DPM:%s] %s", module, paste(..., sep = sep)))
}
