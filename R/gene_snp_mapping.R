# Utility functions for GWAS SNP-to-gene matching.
# These functions reproduce the legacy DPM matching convention:
# gene body +/- window, SNP position as a one-base interval, and data.table::foverlaps().

resolve_column <- function(data, candidates, user_col = NULL, arg_name = "column") {
  if (!is.null(user_col)) {
    if (!is.character(user_col) || length(user_col) != 1L || is.na(user_col)) {
      stop("'", arg_name, "' must be a single column name.")
    }
    if (!user_col %in% names(data)) {
      stop("Column '", user_col, "' specified by '", arg_name, "' was not found.")
    }
    return(user_col)
  }

  lower_names <- tolower(names(data))
  lower_candidates <- tolower(candidates)

  idx <- match(lower_candidates, lower_names)
  idx <- idx[!is.na(idx)]

  if (length(idx) == 0L) {
    stop(
      "Could not infer ", arg_name, ". Tried candidate columns: ",
      paste(candidates, collapse = ", "), "."
    )
  }

  names(data)[idx[1L]]
}

normalize_chrom_for_matching <- function(x, keep_integer = TRUE) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  x <- toupper(x)

  if (keep_integer) {
    x_num <- suppressWarnings(as.integer(x))
    return(x_num)
  }

  x
}

#' Standardize a gene annotation table for SNP-to-gene matching
#'
#' Standardize a gene annotation table using the legacy DPM convention.
#' The output contains gene_name, chrom, start and end. The start/end columns
#' are not extended here; window extension is applied in match_snps_to_genes().
#'
#' @param annotation A data.frame or data.table containing gene annotation.
#' @param gene_col Optional gene name column.
#' @param chr_col Optional chromosome column.
#' @param start_col Optional gene start column.
#' @param end_col Optional gene end column.
#' @param keep_autosomes Logical. Whether to keep chromosomes 1--22 only.
#' @param collapse_genes Logical. Whether to collapse duplicated gene records.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A data.frame with columns gene_name, chrom, start and end.
#'
#' @export
standardize_gene_annotation_for_matching <- function(annotation,
                                                     gene_col = NULL,
                                                     chr_col = NULL,
                                                     start_col = NULL,
                                                     end_col = NULL,
                                                     keep_autosomes = TRUE,
                                                     collapse_genes = FALSE,
                                                     verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  if (!is.data.frame(annotation)) {
    stop("'annotation' must be a data.frame or data.table.")
  }

  dt <- data.table::as.data.table(annotation)

  gene_col <- resolve_column(
    dt,
    candidates = c(
      "gene_name", "gene", "symbol", "hgnc_symbol",
      "external_gene_name", "gene_symbol", "GENE"
    ),
    user_col = gene_col,
    arg_name = "gene_col"
  )

  chr_col <- resolve_column(
    dt,
    candidates = c(
      "seqname", "chrom", "chr", "chromosome",
      "chromosome_name", "gene_chr", "CHR"
    ),
    user_col = chr_col,
    arg_name = "chr_col"
  )

  start_col <- resolve_column(
    dt,
    candidates = c(
      "start", "gene_start", "start_position",
      "tx_start", "bp_start", "START"
    ),
    user_col = start_col,
    arg_name = "start_col"
  )

  end_col <- resolve_column(
    dt,
    candidates = c(
      "end", "gene_end", "stop", "end_position",
      "tx_end", "bp_end", "STOP"
    ),
    user_col = end_col,
    arg_name = "end_col"
  )

  out <- data.table::data.table(
    gene_name = as.character(dt[[gene_col]]),
    chrom = normalize_chrom_for_matching(dt[[chr_col]], keep_integer = TRUE),
    start = suppressWarnings(as.numeric(dt[[start_col]])),
    end = suppressWarnings(as.numeric(dt[[end_col]]))
  )

  swap_idx <- !is.na(out$start) & !is.na(out$end) & out$start > out$end
  if (any(swap_idx)) {
    tmp <- out$start[swap_idx]
    out$start[swap_idx] <- out$end[swap_idx]
    out$end[swap_idx] <- tmp
  }

  keep <- !is.na(out$gene_name) &
    out$gene_name != "" &
    !is.na(out$chrom) &
    is.finite(out$start) &
    is.finite(out$end) &
    out$start > 0 &
    out$end > 0

  out <- out[keep]

  if (keep_autosomes) {
    out <- out[chrom %in% seq_len(22)]
  }

  if (collapse_genes) {
    out <- out[
      ,
      .(
        start = min(start, na.rm = TRUE),
        end = max(end, na.rm = TRUE)
      ),
      by = .(gene_name, chrom)
    ]
  }

  data.table::setorder(out, chrom, start, end, gene_name)

  if (verbose) {
    THREAD_log(
      "Mapping",
      sprintf(
        "Standardized gene annotation for matching: %d gene record(s) across %d chromosome(s).",
        nrow(out), length(unique(out$chrom))
      )
    )
  }

  as.data.frame(out)
}

#' Standardize GWAS summary statistics for SNP-to-gene matching
#'
#' Standardize GWAS summary statistics using the legacy DPM convention.
#' The output preserves original GWAS columns and adds standardized columns:
#' SNP, chrom, pos, start and end.
#'
#' @param gwas A data.frame or data.table containing GWAS summary statistics.
#' @param snp_col Optional SNP ID column.
#' @param chr_col Optional chromosome column.
#' @param pos_col Optional base-pair position column.
#' @param keep_autosomes Logical. Whether to keep chromosomes 1--22 only.
#' @param drop_duplicates Logical. Whether to drop duplicated SNP records.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A data.frame with standardized GWAS columns.
#'
#' @export
standardize_gwas_for_matching <- function(gwas,
                                          snp_col = NULL,
                                          chr_col = NULL,
                                          pos_col = NULL,
                                          keep_autosomes = TRUE,
                                          drop_duplicates = TRUE,
                                          verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  if (!is.data.frame(gwas)) {
    stop("'gwas' must be a data.frame or data.table.")
  }

  dt <- data.table::as.data.table(gwas)

  snp_col <- resolve_column(
    dt,
    candidates = c(
      "SNP", "snp", "snp_name", "snp_id", "rsid",
      "rs_id", "variant_id", "markername", "ID"
    ),
    user_col = snp_col,
    arg_name = "snp_col"
  )

  chr_col <- resolve_column(
    dt,
    candidates = c("CHR", "chrom", "chr", "chromosome", "snp_chr"),
    user_col = chr_col,
    arg_name = "chr_col"
  )

  pos_col <- resolve_column(
    dt,
    candidates = c(
      "POS", "BP", "bp", "pos", "position",
      "base_pair_location", "snp_pos"
    ),
    user_col = pos_col,
    arg_name = "pos_col"
  )

  out <- data.table::copy(dt)

  out[, SNP := as.character(dt[[snp_col]])]
  out[, chrom := normalize_chrom_for_matching(dt[[chr_col]], keep_integer = TRUE)]
  out[, pos := suppressWarnings(as.numeric(dt[[pos_col]]))]
  out[, `:=`(start = pos, end = pos)]

  keep <- !is.na(out$SNP) &
    out$SNP != "" &
    !is.na(out$chrom) &
    is.finite(out$pos) &
    out$pos > 0

  out <- out[keep]

  if (keep_autosomes) {
    out <- out[chrom %in% seq_len(22)]
  }

  if (drop_duplicates) {
    out <- unique(out, by = c("SNP", "chrom", "pos"))
  }

  data.table::setcolorder(
    out,
    c("SNP", "chrom", "pos", "start", "end", setdiff(names(out), c("SNP", "chrom", "pos", "start", "end")))
  )

  if (verbose) {
    THREAD_log(
      "Mapping",
      "Standardized GWAS summary statistics for matching: ",
      nrow(out), " SNP record(s) across ",
      length(unique(out$chrom)), " chromosome(s)."
    )
  }


  as.data.frame(out)
}

#' Match GWAS SNPs to gene windows
#'
#' Match GWAS SNPs to genes using the legacy DPM convention:
#' gene body plus a symmetric window and data.table::foverlaps(type = "within").
#' The output is compatible with the historical gene_snp_matches.rds object:
#' it contains snp_name, snp_pos, gene_name and chrom.
#'
#' @param gwas A GWAS summary statistics data.frame.
#' @param gene_annotation A gene annotation data.frame.
#' @param window Numeric. Symmetric gene-body window in base pairs.
#' @param snp_col Optional SNP ID column in gwas.
#' @param gwas_chr_col Optional chromosome column in gwas.
#' @param pos_col Optional position column in gwas.
#' @param gene_col Optional gene column in gene_annotation.
#' @param gene_chr_col Optional chromosome column in gene_annotation.
#' @param gene_start_col Optional gene start column.
#' @param gene_end_col Optional gene end column.
#' @param keep_autosomes Logical. Whether to keep chromosomes 1--22 only.
#' @param collapse_genes Logical. Whether to collapse duplicated gene records before matching.
#' @param clean_output Logical. Whether to remove columns not used downstream.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A data.frame of SNP-to-gene matches.
#'
#' @export
match_snps_to_genes <- function(gwas,
                                gene_annotation,
                                window = 10000,
                                snp_col = NULL,
                                gwas_chr_col = NULL,
                                pos_col = NULL,
                                gene_col = NULL,
                                gene_chr_col = NULL,
                                gene_start_col = NULL,
                                gene_end_col = NULL,
                                keep_autosomes = TRUE,
                                collapse_genes = FALSE,
                                clean_output = TRUE,
                                verbose = TRUE) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Please install it first.")
  }

  if (!is.numeric(window) || length(window) != 1L || is.na(window) || window < 0) {
    stop("'window' must be a single non-negative numeric value.")
  }

  genes <- standardize_gene_annotation_for_matching(
    annotation = gene_annotation,
    gene_col = gene_col,
    chr_col = gene_chr_col,
    start_col = gene_start_col,
    end_col = gene_end_col,
    keep_autosomes = keep_autosomes,
    collapse_genes = collapse_genes,
    verbose = verbose
  )

  snps <- standardize_gwas_for_matching(
    gwas = gwas,
    snp_col = snp_col,
    chr_col = gwas_chr_col,
    pos_col = pos_col,
    keep_autosomes = keep_autosomes,
    drop_duplicates = TRUE,
    verbose = verbose
  )

  genes_dt <- data.table::as.data.table(genes)
  snps_dt <- data.table::as.data.table(snps)

  genes_dt[, `:=`(
    start = pmax(1, start - window),
    end = end + window
  )]

  data.table::setkey(genes_dt, chrom, start, end)
  data.table::setkey(snps_dt, chrom, start, end)

  if (verbose) {
    THREAD_log("Mapping", "Performing SNP-to-gene overlap matching with foverlaps(type = 'within')...")
  }

  matches <- data.table::foverlaps(
    x = snps_dt,
    y = genes_dt,
    type = "within",
    nomatch = NULL
  )

  if (nrow(matches) == 0L) {
    warning("No SNP-to-gene matches found. Check chromosome naming, genome build and position columns.")
    return(as.data.frame(matches))
  }

  columns_to_remove <- c(
    "source", "feature", "score", "strand", "frame",
    "attribute", "SE", "i.start", "i.end", "P", "AF",
    "SI", "N", "LP"
  )

  if (clean_output) {
    existing_columns <- intersect(columns_to_remove, names(matches))
    if (length(existing_columns) > 0L) {
      matches[, (existing_columns) := NULL]
    }
  }

  if ("pos" %in% names(matches)) {
    data.table::setnames(matches, "pos", "snp_pos")
  }

  if ("SNP" %in% names(matches)) {
    data.table::setnames(matches, "SNP", "snp_name")
  }

  if ("i.pos" %in% names(matches) && !"snp_pos" %in% names(matches)) {
    data.table::setnames(matches, "i.pos", "snp_pos")
  }

  if ("i.SNP" %in% names(matches) && !"snp_name" %in% names(matches)) {
    data.table::setnames(matches, "i.SNP", "snp_name")
  }

  required_cols <- c("gene_name", "snp_name", "snp_pos", "chrom")
  missing_required <- setdiff(required_cols, names(matches))
  if (length(missing_required) > 0L) {
    stop(
      "Internal matching output is missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }

  first_cols <- intersect(
    c("gene_name", "snp_name", "chrom", "start", "end", "snp_pos"),
    names(matches)
  )

  other_cols <- setdiff(names(matches), first_cols)
  matches <- matches[, c(first_cols, other_cols), with = FALSE]

  data.table::setorderv(
    matches,
    intersect(c("gene_name", "chrom", "snp_pos", "snp_name"), names(matches))
  )

  if (verbose) {
    THREAD_log(
      "Mapping",
      sprintf(
        "Matched %d SNP(s) to %d gene(s), producing %d SNP-gene pair(s).",
        length(unique(matches$snp_name)),
        length(unique(matches$gene_name)),
        nrow(matches)
      )
    )
  }

  as.data.frame(matches)
}
