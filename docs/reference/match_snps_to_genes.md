# Match GWAS SNPs to gene windows

Match GWAS SNPs to genes using the legacy DPM convention: gene body plus
a symmetric window and data.table::foverlaps(type = "within"). The
output is compatible with the historical gene_snp_matches.rds object: it
contains snp_name, snp_pos, gene_name and chrom.

## Usage

``` r
match_snps_to_genes(
  gwas,
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
  verbose = TRUE
)
```

## Arguments

- gwas:

  A GWAS summary statistics data.frame.

- gene_annotation:

  A gene annotation data.frame.

- window:

  Numeric. Symmetric gene-body window in base pairs.

- snp_col:

  Optional SNP ID column in gwas.

- gwas_chr_col:

  Optional chromosome column in gwas.

- pos_col:

  Optional position column in gwas.

- gene_col:

  Optional gene column in gene_annotation.

- gene_chr_col:

  Optional chromosome column in gene_annotation.

- gene_start_col:

  Optional gene start column.

- gene_end_col:

  Optional gene end column.

- keep_autosomes:

  Logical. Whether to keep chromosomes 1–22 only.

- collapse_genes:

  Logical. Whether to collapse duplicated gene records before matching.

- clean_output:

  Logical. Whether to remove columns not used downstream.

- verbose:

  Logical. Whether to print progress messages.

## Value

A data.frame of SNP-to-gene matches.
