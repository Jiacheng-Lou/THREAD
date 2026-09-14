# Computing the LD-Corrected Gene-Level Response

## Overview

THREAD models a gene-level response derived from GWAS marginal effect
estimates and local linkage disequilibrium. For each gene,
[`compute_response()`](https://jiacheng-lou.github.io/THREAD/reference/compute_response.md)
returns:

- `final_y`, the LD-corrected response;
- `d0`, `d1`, and `d2`, the heteroscedastic variance coefficients;
- `n_snps`, the number of variants used;
- `status`, the computation status;
- optional LD completeness and PSD-repair diagnostics.

## Locate the installed example inputs

``` r

library(THREAD)

example_dir <- system.file(
  "extdata",
  "example",
  package = "THREAD",
  mustWork = TRUE
)

matches_path <- file.path(
  example_dir,
  "example_gene_snp_matches.rds"
)

gwas_path <- file.path(
  example_dir,
  "example_gwas.tsv.gz"
)

example_X <- readRDS(
  file.path(example_dir, "example_x.rds")
)
```

The example LD file is stored as:

``` text
example_ld_chr1.rds
```

The `{chr}` placeholder in `ld_file_template` allows
[`compute_response()`](https://jiacheng-lou.github.io/THREAD/reference/compute_response.md)
to locate chromosome-specific files.

## Compute the response

``` r

observed_response <- compute_response(
  gene_snp_matches = matches_path,
  gwas_clean = gwas_path,
  ld_dir = example_dir,
  ld_file_template = "example_ld_chr{chr}.rds",
  match_mode = "snp_name",
  genes = rownames(example_X),
  chromosomes = "1",
  n_cores = 1L,
  max_genes_per_chr = NULL,
  max_snps_per_gene = Inf,
  include_diagnostics = TRUE,
  drop_failed = FALSE,
  verbose = TRUE
)
#> [THREAD:Response] compute_response settings: match_mode = snp_name, LD symmetrization = always enabled, PSD repair = trace-preserving eigenvalue clipping, n_cores = 1, set_dt_threads = TRUE.
#> [THREAD:Response] Restricting gene-SNP matches to requested genes: 200 of 200 gene(s).
#> [THREAD:Response] Merged GWAS Beta/SE into matches:  23461  matched SNP-gene row(s).
#> [THREAD:Response] Computing THREAD heteroscedastic response over 1 chromosome(s).
#> [THREAD:Response] Processing chromosome 1 ...
#> [THREAD:Response] Status summary for chromosome 1:
#> 
#>  ok 
#> 200
#> [THREAD:Response] LD diagnostics for chromosome 1: median missing-both fraction = 0.8726; median negative-eigen mass ratio = 0.0156.
#> 
#> [THREAD:Response] Chromosome 1 complete: 200 successful gene(s), 200 returned record(s).
#> [THREAD:Response] Computed response for 200 gene(s). Valid finite final_y: 200.

table(observed_response$status, useNA = "ifany")
#> 
#>  ok 
#> 200
head(observed_response)
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
```

With `drop_failed = FALSE`, the returned object retains both successful
and unsuccessful records. Downstream model preparation uses records with
`status == "ok"`.

## Compare with the frozen example response

The installed example includes a response generated by the validated
implementation. It is used as a regression reference for package
testing.

``` r

expected_response <- readRDS(
  file.path(example_dir, "example_response.rds")
)

observed_response <- observed_response[
  order(observed_response$gene_name),
  ,
  drop = FALSE
]

expected_response <- expected_response[
  order(expected_response$gene_name),
  ,
  drop = FALSE
]

all.equal(
  observed_response[, c(
    "gene_name",
    "n_snps",
    "d0",
    "d1",
    "d2",
    "final_y",
    "status"
  )],
  expected_response[, c(
    "gene_name",
    "n_snps",
    "d0",
    "d1",
    "d2",
    "final_y",
    "status"
  )],
  tolerance = 1e-8,
  check.attributes = FALSE
)
#> [1] TRUE
```

## Variant matching modes

The recommended mode is:

``` r

match_mode = "snp_name"
```

It requires consistent SNP identifiers across the GWAS table, the
SNP-to-gene mapping, and the LD reference.

The fallback mode is:

``` r

match_mode = "snp_pos"
```

Position matching uses chromosome and position jointly. The GWAS summary
statistics, SNP-to-gene mapping, and LD reference must use the same
genome build.

## LD matrix reconstruction

Chromosome-level LD files may store an off-diagonal SNP pair in only one
direction. THREAD recovers the reciprocal entry and constructs a
symmetric gene-level LD matrix.

Pairs absent in both directions are treated as zero and their fraction
is retained as a diagnostic. When sparse pairwise storage produces an
indefinite symmetric matrix, negative eigenvalues are clipped to zero
and the retained spectrum is rescaled to preserve the trace.

## Full-size analyses

A complete analysis typically processes all autosomes and restricts
response calculation to genes present in the final expression matrix.

``` r

checkpoint_dir <- file.path(tempdir(), "THREAD_response_checkpoints")

response <- compute_response(
  gene_snp_matches = gene_snp_matches,
  gwas_clean = gwas_clean,
  ld_dir = ld_directory,
  match_mode = "snp_name",
  genes = rownames(X),
  chromosomes = 1:22,
  n_cores = 20L,
  max_genes_per_chr = NULL,
  max_snps_per_gene = Inf,
  output_dir = checkpoint_dir,
  save_by_chr = TRUE,
  resume = TRUE,
  include_diagnostics = TRUE,
  drop_failed = FALSE,
  verbose = TRUE
)
```

The installed example uses only a small chromosome 1 subset and must not
be interpreted as a complete disease analysis.
