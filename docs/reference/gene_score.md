# Score and rank genes within each module

Scores every gene by an expression-weighted combination of its module's
coefficients, \\h\_{jk} = \sum_p \max(\hat\gamma\_{kp}, 0)\\ x\_{jp}\\,
restricting the sum to positive coefficients so the score reflects only
positive contributions to genetic variance, and ranks genes within each
module in descending order.

This score is an expression-weighted prioritisation, not a direct
measure of disease relevance: because it scales with expression
magnitude, broadly and highly expressed genes (including housekeeping
genes) tend to rank highly. The top-ranked genes should therefore be
treated as a permissive candidate set for set-level analysis rather than
as individually validated targets, and no conclusion should be drawn
from the single highest-ranked gene. Restricting the set to
protein-coding genes is left to the user (see the package
documentation), as it requires an external gene-type annotation.

## Usage

``` r
gene_score(
  partition,
  coefficients,
  x_train,
  top_n = 200,
  out_dir = NULL,
  prefix = "THREAD",
  verbose = TRUE
)
```

## Arguments

- partition:

  A partition from
  [`get_partition`](https://jiacheng-lou.github.io/THREAD/reference/get_partition.md)
  (uses `$assignments`).

- coefficients:

  The per-module coefficient list returned by
  [`test_significance`](https://jiacheng-lou.github.io/THREAD/reference/test_significance.md)
  as `$coefficients` (named by module id, each a vector over cell
  subtypes).

- x_train:

  The gene-by-subtype matrix used to fit the model
  (`thread_data$x_train`); row names are gene identifiers and column
  order must match the coefficient vectors.

- top_n:

  Number of top genes to keep per module. Default `200`.

- out_dir:

  Optional directory; when supplied, one CSV of the top genes is written
  per module. When `NULL` (default) no file is written and the data are
  only returned.

- prefix:

  File-name prefix used when `out_dir` is supplied. Default `"THREAD"`.

- verbose:

  Logical; print a short summary. Default `TRUE`.

## Value

A list with

- scores:

  named numeric vector of the score of every gene (`NA` for genes in a
  module that could not be fit).

- top_genes:

  data frame with `cluster`, `rank`, `gene` and `score`, holding the top
  `top_n` genes per module.

## Examples

``` r
if (FALSE) { # \dontrun{
sig <- test_significance(fit, partition, thread_data)
gs <- gene_score(partition, sig$coefficients, thread_data$x_train, top_n = 200)
head(gs$top_genes)
} # }
```
