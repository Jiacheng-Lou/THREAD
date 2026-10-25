# Standardize a gene annotation table for SNP-to-gene matching

Standardize a gene annotation table using the legacy DPM convention. The
output contains gene_name, chrom, start and end. The start/end columns
are not extended here; window extension is applied in
match_snps_to_genes().

## Usage

``` r
standardize_gene_annotation_for_matching(
  annotation,
  gene_col = NULL,
  chr_col = NULL,
  start_col = NULL,
  end_col = NULL,
  keep_autosomes = TRUE,
  collapse_genes = FALSE,
  verbose = TRUE
)
```

## Arguments

- annotation:

  A data.frame or data.table containing gene annotation.

- gene_col:

  Optional gene name column.

- chr_col:

  Optional chromosome column.

- start_col:

  Optional gene start column.

- end_col:

  Optional gene end column.

- keep_autosomes:

  Logical. Whether to keep chromosomes 1–22 only.

- collapse_genes:

  Logical. Whether to collapse duplicated gene records.

- verbose:

  Logical. Whether to print progress messages.

## Value

A data.frame with columns gene_name, chrom, start and end.
