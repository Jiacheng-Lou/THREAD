# Generate initial gene-cluster labels

Uses gene2vec embeddings when supplied. Genes present in gene2vec are
first clustered by k-means in embedding space. Genes missing from
gene2vec are assigned by KNN voting in the expression-profile space. If
gene2vec is not supplied, the function falls back to k-means on cbind(y,
X).

## Usage

``` r
init_clusters(
  X,
  response,
  gene2vec = NULL,
  K_init = 10,
  k_neighbors = 100,
  gene_col = NULL,
  response_col = "final_y",
  seed = NULL,
  nstart = 25,
  iter.max = 1000,
  tie_method = c("first", "random"),
  min_genes = 10,
  verbose = TRUE
)
```

## Arguments

- X:

  Gene-by-subtype matrix.

- response:

  Response data frame with gene, final_y, d0, d1 and d2.

- gene2vec:

  Optional matrix/data frame/path with genes as rows.

- K_init:

  Initial number of clusters. Default 10.

- k_neighbors:

  Number of neighbours for KNN voting. Default 100.

- gene_col:

  Optional gene column name in response.

- response_col:

  Response column name. Default final_y.

- seed:

  Optional random seed.

- nstart:

  Number of k-means starts. Default 25.

- iter.max:

  Maximum k-means iterations. Default 1000.

- tie_method:

  Tie handling for KNN voting, "first" or "random".

- min_genes:

  Minimum overlap required between X and response.

- verbose:

  Logical.

## Value

Named integer vector of length N, aligned to the intersected genes.
