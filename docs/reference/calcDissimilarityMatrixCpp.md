# Posterior dissimilarity matrix from MCMC label draws

For every pair of genes, returns one minus the fraction of MCMC draws in
which the two genes share a cluster label (i.e. one minus their
posterior co-assignment probability). The diagonal is zero.

## Usage

``` r
calcDissimilarityMatrixCpp(zmat)
```

## Arguments

- zmat:

  An integer matrix of label draws with one row per retained draw and
  one column per gene.

## Value

A symmetric numeric matrix of dissimilarities (genes by genes).
