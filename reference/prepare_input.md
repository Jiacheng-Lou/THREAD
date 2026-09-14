# Prepare THREAD model input

Aligns a gene-by-subtype expression matrix with a gene-level GWAS
response, creates or validates initial cluster labels, recommends prior
hyperparameters when they are not supplied, and returns a `thread_data`
object consumed by
[`run_thread()`](https://jiacheng-lou.github.io/THREAD/reference/run_thread.md).

The response is kept on the raw LD-corrected scale. No log transform or
truncation is applied to y, because y, d0, d1 and d2 must remain on the
same scale as the heteroscedastic likelihood.

## Usage

``` r
prepare_input(
  X,
  response,
  z_init = NULL,
  priors = NULL,
  gene2vec = NULL,
  K_init = 10,
  k_neighbors = 100,
  gene_col = NULL,
  response_col = "final_y",
  x_transform = c("none", "log2p1", "log1p"),
  center_x = FALSE,
  scale_x = FALSE,
  min_genes = 50,
  seed = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- X:

  Gene-by-subtype expression matrix. Rows must be genes.

- response:

  Data frame with gene, final_y, d0, d1 and d2 columns.

- z_init:

  Optional initial labels. Can be a vector aligned to genes, a named
  vector, or a data frame with gene and cluster columns.

- priors:

  Optional list with alpha0, beta0, r and delta. When NULL,
  [`recommend_priors()`](https://jiacheng-lou.github.io/THREAD/reference/recommend_priors.md)
  is called.

- gene2vec:

  Optional gene2vec matrix/data frame/path used only when `z_init` is
  NULL.

- K_init:

  Initial number of clusters when z_init is not supplied.

- k_neighbors:

  Number of neighbours for gene2vec-missing genes.

- gene_col:

  Optional gene column name in response and z_init data frames.

- response_col:

  Response column name. Default final_y.

- x_transform:

  Expression transform: none, log2p1 or log1p. Default none.

- center_x:

  Center expression columns before model fitting. Default FALSE.

- scale_x:

  Scale expression columns before model fitting. Default FALSE.

- min_genes:

  Minimum number of overlapping valid genes. Default 50.

- seed:

  Optional seed used by init_clusters.

- verbose:

  Logical.

- ...:

  Reserved for future extensions.

## Value

A list of class `thread_data` with x_train, y_train, d0, d1, d2, z_init,
gene_names and priors.
