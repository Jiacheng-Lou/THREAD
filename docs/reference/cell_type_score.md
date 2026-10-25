# Integrated cell-subtype importance score (PLACEHOLDER)

Summarises the overall importance of each cell subtype across the whole
model as the total absolute coefficient across modules, \\S_p = \sum_k
\|\hat\gamma\_{kp}\|\\. Absolute values are used because, unlike the
within-module gene score, a negative coefficient still encodes a genuine
association and therefore contributes to a subtype's overall effect. The
score is a magnitude summary only; statistical prioritisation is carried
by the dual-significance criterion, so a subtype is flagged as
significant only when it is dual-significant in at least one module,
graded by its strongest evidence (smallest FDR among its
dual-significant modules).

NOTE: this scoring scheme is a placeholder following the manuscript and
is not yet finalised; it may be revised. Only this function need change.

## Usage

``` r
cell_type_score(sig_result)
```

## Arguments

- sig_result:

  The list returned by
  [`test_significance`](https://jiacheng-lou.github.io/THREAD/reference/test_significance.md)
  (uses `$table`).

## Value

A data frame ordered by `S` (descending), with columns `subtype`, `S`
(integrated score), `significant` (dual in at least one module),
`min_fdr` (smallest FDR among dual-significant modules, `NA` if none)
and `grade` (0/1/2/3 for `q >= 0.05`, `< 0.05`, `< 0.01`, `< 0.001`).

## Examples

``` r
if (FALSE) { # \dontrun{
sig <- test_significance(fit, partition, thread_data)
cts <- cell_type_score(sig)
head(cts)
} # }
```
