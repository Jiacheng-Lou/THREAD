# Recommend empirical-Bayes prior hyperparameters

Estimates a rough coefficient scale through a weighted ridge
approximation. The response is kept on its raw LD-corrected scale; no
log transform is applied to y.

## Usage

``` r
recommend_priors(
  X,
  response,
  z_init = NULL,
  gene_col = NULL,
  response_col = "final_y",
  alpha0 = 2,
  r = 1,
  delta = 1,
  min_beta0 = 1e-16,
  min_var = 1e-12,
  ridge_scale = 1,
  min_cluster_genes = NULL,
  min_genes = 10,
  verbose = TRUE
)
```

## Arguments

- X:

  Gene-by-subtype matrix.

- response:

  Response data frame with gene, final_y, d0, d1 and d2.

- z_init:

  Optional initial labels, vector or data frame.

- gene_col:

  Optional gene column name.

- response_col:

  Response column name. Default final_y.

- alpha0:

  Default inverse-gamma shape for sigma2. Default 2.

- r:

  Default gamma shape for lambda2. Default 1.

- delta:

  Default gamma rate for lambda2. Default 1.

- min_beta0:

  Lower bound for beta0. Default 1e-16.

- min_var:

  Lower bound for approximate variance components. Default 1e-12.

- ridge_scale:

  Multiplicative scale for the ridge penalty. Default 1.

- min_cluster_genes:

  Minimum genes required to use a cluster-specific ridge fit. Defaults
  to max(p + 1, 10).

- min_genes:

  Minimum overlap required between X and response.

- verbose:

  Logical.

## Value

List with alpha0, beta0, r, delta and diagnostic fields.
