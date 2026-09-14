# Dual-significance testing of module / cell-subtype associations

Evaluates each module / cell-subtype association under two independent
criteria on the point-estimate partition and reports those supported by
both. The frequentist criterion refits each module by feasible GLS
(reweighting genes by the precision of their own estimate, with sandwich
standard errors), runs a one-sided test for positive enrichment, and
controls the FDR by Benjamini-Hochberg. The Bayesian criterion aligns
the posterior coefficient draws to the modules and computes the local
false sign rate. An association is dual-significant when `fdr < fdr_thr`
and `lfsr < lfsr_thr`. The reported effect size is the FGLS coefficient,
since the Bayesian posterior mean is shrunk by the Lasso prior.

## Usage

``` r
test_significance(
  fit,
  partition,
  thread_data,
  method = "FGLS",
  fdr_thr = 0.05,
  lfsr_thr = 0.05,
  cred_mass = 0.95,
  verbose = TRUE
)
```

## Arguments

- fit:

  A `thread_fit` from
  [`run_thread`](https://jiacheng-lou.github.io/THREAD/reference/run_thread.md);
  each draw carries `$z` and `$gamma`.

- partition:

  A partition from
  [`get_partition`](https://jiacheng-lou.github.io/THREAD/reference/get_partition.md)
  (uses `$assignments`).

- thread_data:

  The `thread_data` used to fit the model; supplies `x_train`, `y_train`
  and `d0`, `d1`, `d2` on the same scale the sampler used.

- method:

  Frequentist engine, "FGLS" (default) or "OLS".

- fdr_thr:

  FDR threshold. Default `0.05`.

- lfsr_thr:

  lFSR threshold. Default `0.05`.

- cred_mass:

  Credible mass for HPD intervals. Default `0.95`.

- verbose:

  Logical; print progress. Default `TRUE`.

## Value

A list with the following elements:

- table:

  A data frame with one row per module and cellular subtype. It contains
  module dimensions, FGLS estimates, standard errors, one-sided P
  values, FDR values, posterior summaries, local false sign rates, and
  the dual-significance indicator.

- coefficients:

  A named list of per-module FGLS coefficient vectors. These
  coefficients are used by
  [`gene_score`](https://jiacheng-lou.github.io/THREAD/reference/gene_score.md).

- module_diagnostics:

  A data frame with one row per module, reporting the number of genes,
  number of predictors, design rank, and whether frequentist inference
  was estimable.

- method:

  The frequentist inference method used.

- thresholds:

  A list containing the FDR threshold, lFSR threshold, and credible
  mass.

## Examples

``` r
if (FALSE) { # \dontrun{
sig <- test_significance(fit, partition, thread_data)
subset(sig$table, dual_significant)
} # }
```
