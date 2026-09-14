# Fit the THREAD model by MCMC

Runs the Dirichlet-process-mixture sampler on a prepared `thread_data`
object and returns a `thread_fit` object holding the posterior draws,
diagnostics and settings. Cluster labels are updated by Neal's Algorithm
8 (auxiliary components proposed from the prior); the coefficient
vectors are updated by stepping-out slice sampling; the Bayesian-Lasso
scale parameters have closed-form conjugate updates.

All hyperparameters are explicit arguments or are read from
`thread_data$priors`; the function holds no hidden state. Posterior
draws are stored separately from the input data because they can be
large; the downstream functions therefore take both the `thread_fit` and
the `thread_data`.

## Usage

``` r
run_thread(
  thread_data,
  alpha = 0.1,
  m = 3,
  K_init = 10,
  n_iter = 30000,
  burnin = 25000,
  thin = 1,
  priors = NULL,
  seed = NULL,
  adaptive_ridge = TRUE,
  adaptive_w = TRUE,
  save_scales = FALSE,
  verbose = TRUE
)
```

## Arguments

- thread_data:

  A list produced by
  [`prepare_input()`](https://jiacheng-lou.github.io/THREAD/reference/prepare_input.md)
  containing at least `x_train` (gene-by-subtype matrix), `y_train`
  (response vector), `d0`, `d1`, `d2` (heteroscedastic-variance
  coefficients) and, optionally, `z_init` (starting labels),
  `gene_names` and `priors`.

- alpha:

  Dirichlet-process concentration parameter. Default `0.1`.

- m:

  Number of auxiliary components in Neal's Algorithm 8. Default `3`.

- K_init:

  Number of starting clusters when `z_init` is not supplied. Default
  `10`. When `z_init` is supplied this is taken from the labels.

- n_iter:

  Total number of MCMC iterations. Default `30000`.

- burnin:

  Number of initial iterations discarded. Default `25000`.

- thin:

  Keep one draw every `thin` post-burn-in iterations. Default set to `1`
  to keep every draw (larger output).

- priors:

  Optional list overriding `thread_data$priors`; must contain `alpha0`,
  `beta0`, `r`, `delta`.

- seed:

  Optional integer random seed for reproducibility. Default `NULL`.

- adaptive_ridge:

  Logical; passed to the warm-start initialisation. Default `TRUE`.

- adaptive_w:

  Logical; passed to the slice sampler. Default `TRUE`.

- save_scales:

  Logical. Whether to retain tau2 and lambda2 in every posterior draw.
  Default FALSE because these objects can substantially increase memory
  and storage use.

- verbose:

  Logical; show a progress bar and status messages. Default `TRUE`.

## Value

An object of class `thread_fit`: a list with

- samples:

  list with one element per retained draw, each a list of `z` (labels),
  `gamma` (K-by-p coefficient matrix), `sigma2` and `K`.

- diagnostics:

  list with the per-draw cluster count (`cluster_counts`) and the
  wall-clock `runtime`.

- settings:

  the hyperparameters and run controls used.

- dims:

  `n_genes`, `n_subtypes`, `gene_names`, `subtype_names`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- run_thread(thread_data, n_iter = 30000, burnin = 25000, thin = 1, seed = 123)
table(fit$diagnostics$cluster_counts)
} # }
```
