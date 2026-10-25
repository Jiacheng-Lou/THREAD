# Standardize GWAS summary statistics for SNP-to-gene matching

Standardize GWAS summary statistics using the legacy DPM convention. The
output preserves original GWAS columns and adds standardized columns:
SNP, chrom, pos, start and end.

## Usage

``` r
standardize_gwas_for_matching(
  gwas,
  snp_col = NULL,
  chr_col = NULL,
  pos_col = NULL,
  keep_autosomes = TRUE,
  drop_duplicates = TRUE,
  verbose = TRUE
)
```

## Arguments

- gwas:

  A data.frame or data.table containing GWAS summary statistics.

- snp_col:

  Optional SNP ID column.

- chr_col:

  Optional chromosome column.

- pos_col:

  Optional base-pair position column.

- keep_autosomes:

  Logical. Whether to keep chromosomes 1–22 only.

- drop_duplicates:

  Logical. Whether to drop duplicated SNP records.

- verbose:

  Logical. Whether to print progress messages.

## Value

A data.frame with standardized GWAS columns.
