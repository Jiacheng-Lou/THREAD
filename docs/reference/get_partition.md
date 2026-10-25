# Point estimate of the gene partition from a fitted THREAD model

Turns posterior label draws into a single point-estimate partition.
Because mixture labels are exchangeable, the function first constructs
the posterior dissimilarity matrix `1 - PSM`. It then computes the
trajectory mode of the number of major clusters, ignoring transient
micro-clusters with at most `min_cluster_size` genes. The final K is
selected by PAM silhouette within a local window around this trajectory
mode.

## Usage

``` r
get_partition(
  fit,
  min_cluster_size = 10,
  mode_window = 3,
  min_silhouette = NULL,
  return_dissimilarity = FALSE,
  verbose = TRUE
)
```

## Arguments

- fit:

  A `thread_fit` object from
  [`run_thread`](https://jiacheng-lou.github.io/THREAD/reference/run_thread.md).
  For backward compatibility, objects with `$results` instead of
  `$samples` are also accepted.

- min_cluster_size:

  Clusters with at most this many genes are excluded when counting the
  number of major modules along the chain. Default `10`.

- mode_window:

  Non-negative integer. Candidate K values are restricted to
  `mode_k - mode_window` through `mode_k + mode_window`, clipped to the
  valid range. Default `3`.

- min_silhouette:

  Optional minimum average silhouette width. If supplied and the best
  silhouette is below this threshold, the result collapses to `K = 1`.
  Default `NULL`.

- return_dissimilarity:

  Logical; whether to include the full gene-by-gene dissimilarity matrix
  in the output. This can be large, so the default is `FALSE`.

- verbose:

  Logical; print trajectory and selection diagnostics. Default `TRUE`.

## Value

A list with

- assignments:

  Integer vector of length N giving each gene's module.

- K:

  Number of modules in the point estimate.

- silhouette:

  List containing candidate K values, silhouette widths, selected K and
  maximum width.

- mode_k:

  Trajectory mode of the major-cluster count.

- major_cluster_count:

  Major-cluster count for every retained draw.

- dissimilarity:

  The full `1 - PSM` matrix, only when `return_dissimilarity = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
partition <- get_partition(fit, min_cluster_size = 10, mode_window = 3)
partition$K
table(partition$assignments)
} # }
```
