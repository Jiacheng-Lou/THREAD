# Read a gene annotation table

Read a gene annotation table from a local text file. This is a
lightweight reader only; column standardization is performed by
standardize_gene_annotation().

## Usage

``` r
read_gene_annotation(path, ...)
```

## Arguments

- path:

  Character string. Path to a gene annotation file.

- ...:

  Additional arguments passed to data.table::fread().

## Value

A data.frame containing the raw gene annotation table.
