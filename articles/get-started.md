# Getting Started with THREAD

## Overview

THREAD integrates an LD-corrected gene-level GWAS response with a
gene-by-cellular-label expression matrix. It uses a heteroscedastic
Dirichlet process mixture to infer latent gene modules and
module-specific cellular coefficients.

The shortest workflow starts from:

1.  a gene-by-cellular-label expression matrix, `X`;
2.  a response table containing `gene_name`, `final_y`, `d0`, `d1`, and
    `d2`.

The installed example data are deliberately small. They are intended to
demonstrate the software interface and must not be used for biological
interpretation.

## Load the installed example data

``` r

library(THREAD)

example_dir <- system.file(
  "extdata",
  "example",
  package = "THREAD",
  mustWork = TRUE
)

example_X <- readRDS(
  file.path(example_dir, "example_x.rds")
)

example_response <- readRDS(
  file.path(example_dir, "example_response.rds")
)

dim(example_X)
#> [1] 200  31
head(example_response)
#>   gene_name n_snps           d0           d1        d2      final_y status
#> 1     ACAP3    123 2.676094e-08 2.327622e-04 0.6352512 0.0006409818     ok
#> 2   ADPRHL2     58 1.427176e-07 5.672954e-04 0.6411754 0.0031773968     ok
#> 3      AGO1    120 2.326618e-08 2.693607e-04 0.8702472 0.0003077566     ok
#> 4      AGO4    110 9.465677e-09 1.599227e-04 0.7608213 0.0011992010     ok
#> 5      AGRN    168 1.834882e-09 7.207823e-05 0.7747462 0.0001823946     ok
#> 6     AHDC1    114 1.699830e-08 2.402827e-04 0.9454711 0.0002825845     ok
#>   ld_forward_coverage ld_reverse_coverage ld_both_present_fraction
#> 1          0.04804745          0.04804745                        0
#> 2          0.05353902          0.05353902                        0
#> 3          0.06001401          0.06001401                        0
#> 4          0.07180984          0.07180984                        0
#> 5          0.06750784          0.06750784                        0
#> 6          0.04595560          0.04595560                        0
#>   ld_missing_both_fraction ld_max_reciprocal_difference ld_max_asymmetry
#> 1                0.9039051                           NA                0
#> 2                0.8929220                           NA                0
#> 3                0.8799720                           NA                0
#> 4                0.8563803                           NA                0
#> 5                0.8649843                           NA                0
#> 6                0.9080888                           NA                0
#>   min_eigenvalue_raw n_negative_eigenvalues negative_eigen_mass_ratio
#> 1      -2.273143e+00                     13              4.537463e-02
#> 2      -4.375298e-01                      2              1.134668e-02
#> 3      -5.888512e-01                      3              1.243750e-02
#> 4      -1.500979e+00                      7              3.124168e-02
#> 5      -1.469352e+00                     13              4.184783e-02
#> 6      -7.490211e-07                      2              7.303095e-09
#>   psd_repaired
#> 1         TRUE
#> 2         TRUE
#> 3         TRUE
#> 4         TRUE
#> 5         TRUE
#> 6         TRUE
table(example_response$status, useNA = "ifany")
#> 
#>  ok 
#> 200
```

The rows of `example_X` are genes and the columns are cellular labels.
The response is retained on its original LD-corrected scale. Negative
`final_y` values, when present, are not truncated or log transformed.

## Prepare model input

[`prepare_input()`](https://jiacheng-lou.github.io/THREAD/reference/prepare_input.md)
aligns genes between the expression matrix and the response, retains
response records with `status == "ok"`, generates initial cluster
labels, and recommends prior hyperparameters when priors are not
supplied.

``` r

thread_data <- prepare_input(
  X = example_X,
  response = example_response,
  K_init = 3L,
  min_genes = 10L,
  seed = 2028L,
  verbose = TRUE
)
#> [THREAD:Input] Using 200 response row(s) with status == 'ok'; excluded 0 non-ok row(s).
#> [THREAD:Input] Aligned 200 genes across X (200 total) and response (200 total).
#> [THREAD:Input] Initialisation: k-means fallback on cbind(y, X), K = 3.
#> [THREAD:Input] Prior recommendation: cluster-wise weighted ridge over K = 3 initial clusters.
#> [THREAD:Input] Recommended priors: alpha0 = 2, beta0 = 1.277495e-07, r = 1, delta = 1, prior_sd = 2.918327e-04.
#> [THREAD:Input] Prepared thread_data: 200 genes x 31 subtypes, K_init = 3.

dim(thread_data$x_train)
#> [1] 200  31
table(thread_data$z_init)
#> 
#>   1   2   3 
#>  69  22 109
thread_data$priors
#> $alpha0
#> [1] 2
#> 
#> $beta0
#> [1] 1.277495e-07
#> 
#> $r
#> [1] 1
#> 
#> $delta
#> [1] 1
#> 
#> $beta_variance
#> [1] 2.55499e-07
#> 
#> $base_lambda
#> [1] 1213.694
#> 
#> $mode
#> [1] "cluster-wise"
#> 
#> $prior_sd
#> [1] 0.0002918327
```

Gene2vec initialization is optional. When no embedding is supplied,
THREAD uses a k-means fallback based on the aligned response and
expression matrix.

## Fit the THREAD model

The settings below are intentionally short so that the vignette builds
quickly. A production analysis should use substantially longer chains
and should evaluate stability across seeds.

``` r

fit <- run_thread(
  thread_data = thread_data,
  alpha = 0.1,
  m = 2L,
  n_iter = 60L,
  burnin = 30L,
  thin = 1L,
  seed = 2028L,
  save_scales = FALSE,
  verbose = TRUE
)
#> [THREAD:MCMC] Using K_init = 3 from supplied z_init (ignoring K_init = 10).
#> [THREAD:MCMC] 200 genes x 31 subtypes | alpha = 0.1, m = 2, K_init = 3
#> [THREAD:MCMC] iterations = 60, burnin = 30, thin = 1, expected retained draws = 30
#> [THREAD:MCMC] auxiliary-cluster proposal prior_sd = 0.000291833
#> [THREAD:MCMC] save_scales =  FALSE
#>   |                                                                              |                                                                      |   0%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |====                                                                  |   5%  |                                                                              |=====                                                                 |   7%  |                                                                              |======                                                                |   8%  |                                                                              |=======                                                               |  10%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |==========                                                            |  15%  |                                                                              |============                                                          |  17%  |                                                                              |=============                                                         |  18%  |                                                                              |==============                                                        |  20%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  23%  |                                                                              |==================                                                    |  25%  |                                                                              |===================                                                   |  27%  |                                                                              |====================                                                  |  28%  |                                                                              |=====================                                                 |  30%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |========================                                              |  35%  |                                                                              |==========================                                            |  37%  |                                                                              |===========================                                           |  38%  |                                                                              |============================                                          |  40%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |================================                                      |  45%  |                                                                              |=================================                                     |  47%  |                                                                              |==================================                                    |  48%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |======================================                                |  55%  |                                                                              |========================================                              |  57%  |                                                                              |=========================================                             |  58%  |                                                                              |==========================================                            |  60%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  63%  |                                                                              |==============================================                        |  65%  |                                                                              |===============================================                       |  67%  |                                                                              |================================================                      |  68%  |                                                                              |=================================================                     |  70%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================                |  77%  |                                                                              |=======================================================               |  78%  |                                                                              |========================================================              |  80%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |============================================================          |  85%  |                                                                              |=============================================================         |  87%  |                                                                              |==============================================================        |  88%  |                                                                              |===============================================================       |  90%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |==================================================================    |  95%  |                                                                              |====================================================================  |  97%  |                                                                              |===================================================================== |  98%  |                                                                              |======================================================================| 100%
#> [THREAD:MCMC] MCMC complete in 0.8 secs. Retained 30 draws.
#> [THREAD:MCMC] Occupied modules across retained draws: min = 3, median = 4.0, max = 4.

length(fit$samples)
#> [1] 30
summary(fit$diagnostics$cluster_counts)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>     3.0     3.0     4.0     3.6     4.0     4.0
```

## Obtain a point-estimate partition

THREAD summarizes posterior label draws through a posterior similarity
matrix. The point-estimate partition is selected by PAM silhouette
within a local window around the posterior-supported number of major
modules.

``` r

partition <- get_partition(
  fit = fit,
  min_cluster_size = 1L,
  mode_window = 1L,
  return_dissimilarity = FALSE,
  verbose = TRUE
)
#> [THREAD:Partition] Trajectory mode of major clusters: 4 (overall 50.0%, last 30 draws 50.0%).
#> [THREAD:Partition] Silhouette K search window: 3, 4, 5 (mode_k = 4, mode_window = 1).
#> [THREAD:Partition] Selected K = 3 by maximum average silhouette width within the mode-local window (0.538).

partition$K
#> [1] 3
table(partition$assignments)
#> 
#>   1   2   3 
#> 141  39  20
```

## Test module-cellular-label associations

[`test_significance()`](https://jiacheng-lou.github.io/THREAD/reference/test_significance.md)
combines:

- feasible generalized least squares with BH FDR control;
- Bayesian local false sign rate;
- a dual-significance rule requiring both criteria.

``` r

significance <- test_significance(
  fit = fit,
  partition = partition,
  thread_data = thread_data,
  fdr_thr = 0.05,
  lfsr_thr = 0.05,
  verbose = TRUE
)
#> [THREAD:Inference] Frequentist inference (FGLS) over 3 module(s).
#> [THREAD:Inference] Module 3 has 20 genes and 31 predictors; frequentist inference is underdetermined. Coefficients are retained, but SE, P value and FDR are set to NA.
#> [THREAD:Inference] Frequentist inference was estimable for 2 of 3 module(s).
#> [THREAD:Inference] Bayesian inference (lFSR): aligning posterior draws to the partition.
#> [THREAD:Inference] Dual-significant (FDR < 0.05 and lFSR < 0.05): 3 module-subtype association(s).

head(significance$table)
#>   cluster subtype module_n_genes n_predictors design_rank fgls_estimable
#> 1       1 Astro-1            141           31          31           TRUE
#> 2       1 Astro-2            141           31          31           TRUE
#> 3       1 Astro-3            141           31          31           TRUE
#> 4       1 Astro-4            141           31          31           TRUE
#> 5       1 Astro-5            141           31          31           TRUE
#> 6       1 Astro-6            141           31          31           TRUE
#>      gamma_fgls      se_fgls  pval_fgls       fdr sig_fdr gamma_post_mean
#> 1  6.594922e-05 4.433371e-05 0.06986349 0.4640932   FALSE    2.448297e-05
#> 2 -1.926641e-05 4.401078e-05 0.66879278 1.0000000   FALSE    1.050650e-05
#> 3 -8.716004e-06 4.104011e-05 0.58389729 1.0000000   FALSE    3.321492e-05
#> 4  1.740479e-05 2.287822e-05 0.22421459 0.8019983   FALSE   -1.674448e-05
#> 5  5.287563e-05 5.299265e-05 0.16028513 0.7569698   FALSE    6.122777e-06
#> 6  1.732439e-05 3.200594e-05 0.29470228 0.9788326   FALSE    4.317505e-06
#>       hpd_lower     hpd_upper      lfsr sig_lfsr dual_significant
#> 1  1.992231e-05  3.258222e-05 0.0000000     TRUE            FALSE
#> 2  4.973565e-06  1.574071e-05 0.0000000     TRUE            FALSE
#> 3  2.641793e-05  4.039899e-05 0.0000000     TRUE            FALSE
#> 4 -2.713681e-05 -5.422519e-06 1.0000000    FALSE            FALSE
#> 5 -1.785176e-06  1.141155e-05 0.1333333    FALSE            FALSE
#> 6 -1.207750e-05  1.547528e-05 0.3000000    FALSE            FALSE
significance$module_diagnostics
#>   cluster module_n_genes n_predictors design_rank fgls_estimable
#> 1       1            141           31          31           TRUE
#> 2       2             39           31          31           TRUE
#> 3       3             20           31          20          FALSE
subset(significance$table, dual_significant)
#>    cluster subtype module_n_genes n_predictors design_rank fgls_estimable
#> 13       1 Micro-3            141           31          31           TRUE
#> 32       2 Astro-1             39           31          31           TRUE
#> 55       2 Oligo-3             39           31          31           TRUE
#>      gamma_fgls      se_fgls    pval_fgls          fdr sig_fdr gamma_post_mean
#> 13 8.537816e-05 1.852199e-05 5.469186e-06 0.0005086343    TRUE    5.540312e-05
#> 32 6.908162e-04 1.404908e-04 5.840909e-04 0.0181068184    TRUE    6.219437e-05
#> 55 1.546473e-03 3.368359e-04 8.879157e-04 0.0206440398    TRUE    9.571439e-05
#>       hpd_lower    hpd_upper       lfsr sig_lfsr dual_significant
#> 13 3.672685e-05 7.224654e-05 0.00000000     TRUE             TRUE
#> 32 6.664876e-06 2.152448e-04 0.03333333     TRUE             TRUE
#> 55 1.305387e-05 1.304463e-04 0.00000000     TRUE             TRUE
```

The example contains only 24 genes but 31 cellular predictors.
Consequently, some or all example modules can be underdetermined for
frequentist inference. For such modules, THREAD retains descriptive
coefficients but sets the frequentist standard error, P value, and FDR
to `NA`. They cannot be declared dual-significant. This behavior is
intentional and allows the example to test the package interface without
presenting the miniature dataset as a valid biological analysis.

## Rank genes and cellular labels

``` r

gene_results <- gene_score(
  partition = partition,
  coefficients = significance$coefficients,
  x_train = thread_data$x_train,
  top_n = 10L,
  verbose = TRUE
)

head(gene_results$top_genes)
#>        cluster rank   gene       score
#> RPL11        1    1  RPL11 0.006323269
#> SFPQ         1    2   SFPQ 0.006246336
#> HP1BP3       1    3 HP1BP3 0.006234229
#> RPS8         1    4   RPS8 0.006167004
#> ARID1A       1    5 ARID1A 0.005951408
#> PRDX1        1    6  PRDX1 0.005866957

cellular_results <- cell_type_score(significance)
head(cellular_results)
#>    subtype            S significant    min_fdr grade
#> 1  Oligo-3 0.0016725200        TRUE 0.02064404     1
#> 2 Neuron-5 0.0014970239       FALSE         NA     0
#> 3    OPC-1 0.0013618312       FALSE         NA     0
#> 4  Oligo-5 0.0012350096       FALSE         NA     0
#> 5  Oligo-6 0.0010775965       FALSE         NA     0
#> 6  Astro-1 0.0007567654        TRUE 0.01810682     1
```

The gene score is an expression-weighted prioritization score rather
than an independent significance test. Biological conclusions should
combine module enrichment, transcription-factor analysis, external
validation, and other downstream evidence.
