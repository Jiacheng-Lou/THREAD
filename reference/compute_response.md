# Compute LD-based heteroscedastic gene-level THREAD response

Compute LD-based heteroscedastic gene-level THREAD response

## Usage

``` r
compute_response(
  gene_snp_matches = NULL,
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
  psd_tolerance = 1e-08,
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
  verbose = TRUE
)
```

## Arguments

- gene_snp_matches:

  A data.frame/data.table or file path to gene-SNP matches.

- matches:

  Deprecated alias for gene_snp_matches.

- gwas_clean:

  Optional data.frame/data.table or file path containing beta and SE.

- ld_dir:

  Directory containing chromosome-level LD RDS files.

- match_mode:

  Variant-key mode. Use "snp_name" (default and recommended) for
  rsID-based matching, or "snp_pos" for chromosome-position matching.
  The selected mode controls both GWAS matching and LD lookup.

- genes:

  Optional character vector restricting computation to selected genes.

- chromosomes:

  Optional chromosomes to process.

- n_cores:

  Number of forked workers for gene-level parallelism.

- set_dt_threads:

  Logical. If TRUE, data.table threads are set to 1 during computation.

- max_genes_per_chr:

  Optional maximum number of genes per chromosome for smoke tests.

- max_snps_per_gene:

  Maximum SNPs per gene. Larger genes are skipped.

- psd_tolerance:

  Eigenvalues below -psd_tolerance are counted as negative.

- include_diagnostics:

  Logical. Include LD completeness, symmetry and PSD diagnostic columns
  in the returned table.

- output_file:

  Optional final merged RDS output path.

- output_dir:

  Optional directory for chromosome-level checkpoint files.

- save_by_chr:

  Logical. Whether to save one response file per chromosome.

- resume:

  Logical. Whether to reuse existing chromosome-level files when
  save_by_chr is TRUE.

- snp_col:

  SNP-ID column in matches.

- snp_pos_col:

  SNP-position column in matches.

- gene_col:

  Gene column in matches.

- chr_col:

  Chromosome column in matches.

- gwas_snp_col:

  SNP-ID column in gwas_clean.

- gwas_pos_col:

  Position column in gwas_clean. If NULL, inferred for match_mode =
  "snp_pos".

- gwas_chr_col:

  Chromosome column in gwas_clean. Required for position matching.

- beta_col:

  Beta column in gwas_clean.

- se_col:

  Standard-error column in gwas_clean.

- ld_file_template:

  Character string defining chromosome-specific LD files under `ld_dir`.
  The token `chr`, enclosed in braces, is replaced by the chromosome
  value. The default corresponds to one RDS file named by chromosome.

- ld_snp_a_col:

  SNP_A column in SNP-name LD files.

- ld_snp_b_col:

  SNP_B column in SNP-name LD files.

- ld_bp_a_col:

  BP_A column in position LD files.

- ld_bp_b_col:

  BP_B column in position LD files.

- ld_r_col:

  LD correlation column.

- drop_failed:

  Logical. Whether to remove genes with non-ok status.

- verbose:

  Logical. Whether to print progress messages.

## Value

A data.frame containing gene_name, n_snps, d0, d1, d2, final_y, status
and, optionally, LD diagnostic columns.
