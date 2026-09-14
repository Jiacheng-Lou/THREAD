# Read a gene2vec embedding file

Read an external gene embedding file for gene2vec-based initialization
in THREAD. The file should contain one gene identifier column followed
by numeric embedding dimensions. A Word2Vec-style metadata header, such
as "25442 200", is automatically detected and skipped.

## Usage

``` r
read_gene2vec(path, n_max = Inf, gene_col = 1L, skip = NULL, verbose = TRUE)
```

## Arguments

- path:

  Character string. Path to the gene2vec embedding file.

- n_max:

  Integer or Inf. Maximum number of rows to read after skipping a
  possible metadata header. Use Inf to read the full file.

- gene_col:

  Integer. Column index containing gene identifiers.

- skip:

  Integer or NULL. Number of rows to skip. If NULL, the function
  automatically detects a Word2Vec-style metadata header.

- verbose:

  Logical. Whether to print progress messages.

## Value

A numeric matrix with genes as row names and embedding dimensions as
columns.
