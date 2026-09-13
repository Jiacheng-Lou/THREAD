
<!-- README.md is generated from README.Rmd. Edit README.Rmd, then run devtools::build_readme(). -->

# THREAD

[![R-CMD-check](https://github.com/Jiacheng-Lou/DPM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jiacheng-Lou/DPM/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Jiacheng-Lou/DPM/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Jiacheng-Lou/DPM/actions/workflows/pkgdown.yaml)

**THREAD** (*Trait Heterogeneity Regression through Expression
Annotation and Dirichlet processes*) is an R package for integrating
LD-corrected gene-level GWAS response statistics with
gene-by-cellular-label expression profiles. THREAD infers latent gene
modules and estimates module-specific cellular coefficients under a
heteroscedastic Bayesian nonparametric model.

The package supports a modular workflow. Users can start from raw GWAS
and single-cell inputs, from an externally prepared SNP-to-gene map,
from an existing gene-level response, or directly from a prepared THREAD
model object.

## Main capabilities

THREAD provides:

- preprocessing of GWAS summary statistics;
- aggregation of single-cell expression into a gene-by-cellular-label
  matrix;
- SNP-to-gene mapping using the annotated gene body plus a configurable
  genomic window;
- LD-aware construction of the gene-level response and heteroscedastic
  variance terms;
- SNP-name and chromosome-position matching modes;
- automatic reciprocal completion and symmetrization of sparse LD
  records;
- trace-preserving PSD repair by eigenvalue clipping;
- gene2vec-based or k-means-based initialization;
- Bayesian Dirichlet process mixture fitting;
- posterior-similarity-based gene partitioning;
- frequentist and Bayesian module-cellular-label inference;
- downstream gene and cellular-label prioritization scores.

## Installation

Install the development version from GitHub:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github(
  "Jiacheng-Lou/DPM",
  build_vignettes = TRUE
)

Load the package:

```r
library(THREAD)
```

## Quick start with the installed example data

THREAD includes a compact AD-derived example containing 200 genes and 31
cellular labels. The example is intended only for software testing and
interface demonstration; it is not a complete disease analysis and
should not be interpreted biologically.

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
head(example_response)
```

Prepare aligned model input:

``` r
thread_data <- prepare_input(
  X = example_X,
  response = example_response,
  K_init = 3L,
  min_genes = 10L,
  seed = 2028L,
  verbose = TRUE
)
```

Run THREAD:

``` r
fit <- run_thread(
  thread_data = thread_data,
  alpha = 0.1,
  m = 3L,
  n_iter = 30000L,
  burnin = 25000L,
  thin = 1L,
  seed = 2028L,
  save_scales = FALSE,
  verbose = TRUE
)
```

The chain above represents a production-style template. Short chains may
be used only for installation checks and software tests, not for
scientific inference.

Obtain a point-estimate partition:

``` r
partition <- get_partition(
  fit = fit,
  min_cluster_size = 10L,
  mode_window = 3L,
  verbose = TRUE
)

partition$K
table(partition$assignments)
```

Test module-cellular-label associations:

``` r
significance <- test_significance(
  fit = fit,
  partition = partition,
  thread_data = thread_data,
  fdr_thr = 0.05,
  lfsr_thr = 0.05,
  cred_mass = 0.95,
  verbose = TRUE
)

head(significance$table)
subset(significance$table, dual_significant)
```

Calculate downstream prioritization scores:

``` r
gene_results <- gene_score(
  partition = partition,
  coefficients = significance$coefficients,
  x_train = thread_data$x_train,
  top_n = 20L,
  verbose = TRUE
)

cellular_results <- cell_type_score(significance)

head(gene_results$top_genes)
head(cellular_results)
```

## Complete GWAS-to-THREAD workflow

### 1. Construct a gene-by-cellular-label expression matrix

Users who already have an aligned numeric matrix can skip this step. The
required orientation is genes in rows and cellular labels in columns.

``` r
X <- preprocess_scrna(
  object = single_cell_object,
  subtype_col = "cell_subtype",
  already_normalized = FALSE
)
```

`preprocess_scrna()` is an aggregation helper. Upstream quality control,
batch correction, integration, and cell annotation remain the user’s
responsibility.

### 2. Preprocess GWAS summary statistics

``` r
gwas_clean <- preprocess_gwas(
  sumstats = "path/to/gwas.tsv.gz",
  col_map = c(
    chrom = "CHR",
    pos = "BP",
    rsid = "SNP",
    beta = "BETA",
    se = "SE"
  )
)
```

The standardized core fields are:

``` text
chrom, pos, rsid, beta, se
```

### 3. Map SNPs to genes

``` r
gene_annotation <- read_gene_annotation(
  "path/to/gene_annotation.tsv.gz"
)

gene_snp_matches <- match_snps_to_genes(
  gwas = gwas_clean,
  gene_annotation = gene_annotation,
  window = 10000,
  verbose = TRUE
)
```

The standard mapping output contains:

``` text
gene_name, snp_name, snp_pos, chrom
```

A SNP can map to multiple genes when extended gene windows overlap.

### 4. Compute the LD-corrected gene-level response

``` r
checkpoint_dir <- file.path(tempdir(), "THREAD_response_checkpoints")

response <- compute_response(
  gene_snp_matches = gene_snp_matches,
  gwas_clean = gwas_clean,
  ld_dir = "path/to/ld_directory",
  match_mode = "snp_name",
  genes = rownames(X),
  chromosomes = 1:22,
  n_cores = 20L,
  max_genes_per_chr = NULL,
  max_snps_per_gene = Inf,
  include_diagnostics = TRUE,
  output_dir = checkpoint_dir,
  save_by_chr = TRUE,
  resume = TRUE,
  drop_failed = FALSE,
  verbose = TRUE
)
```

Each successful gene record contains:

``` text
gene_name, n_snps, final_y, d0, d1, d2, status
```

Additional LD completeness and PSD-repair diagnostics are retained when
`include_diagnostics = TRUE`.

THREAD retains negative values of `final_y`. They are not truncated or
log transformed.

### 5. Prepare, fit, partition, and infer

``` r
thread_data <- prepare_input(
  X = X,
  response = response,
  gene2vec = NULL,
  K_init = 10L,
  seed = 2028L
)

fit <- run_thread(
  thread_data = thread_data,
  alpha = 0.1,
  m = 3L,
  n_iter = 30000L,
  burnin = 25000L,
  thin = 1L,
  seed = 2028L
)

partition <- get_partition(fit)

significance <- test_significance(
  fit = fit,
  partition = partition,
  thread_data = thread_data
)
```

## LD reference data

Full genome-wide LD files are not bundled with the package. Users must
supply a compatible, ancestry-matched external LD reference in the
documented input format; see the *Computing the LD-Corrected Gene-Level
Response* vignette.

For SNP-name matching, chromosome-specific LD files should contain:

``` text
SNP_A, SNP_B, R
```

For position matching, they should additionally contain:

``` text
BP_A, BP_B
```

A typical directory contains one file per chromosome:

``` text
ld_directory/
├── 1.rds
├── 2.rds
├── ...
└── 22.rds
```

Sparse LD storage may contain only one orientation of an off-diagonal
pair. THREAD automatically restores the reciprocal value and constructs
a symmetric gene-level matrix. Pairs absent in both directions are
treated as zero and reported through diagnostics.

All GWAS, gene annotation, SNP-to-gene mapping, and LD resources used in
position mode must use the same genome build.

## Supported workflow entry points

| Available inputs | First THREAD function |
|----|----|
| Raw single-cell data and raw GWAS | `preprocess_scrna()` and `preprocess_gwas()` |
| Gene-by-label matrix and raw GWAS | `preprocess_gwas()` |
| Clean GWAS and gene annotation | `match_snps_to_genes()` |
| SNP-gene map, clean GWAS, and LD | `compute_response()` |
| Expression matrix and gene-level response | `prepare_input()` |
| Prepared `thread_data` object | `run_thread()` |
| Fitted `thread_fit` object | `get_partition()` |
| Fit, partition, and model data | `test_significance()` |
| Significance output | `gene_score()` and `cell_type_score()` |

## Statistical interpretation

THREAD uses a heteroscedastic module-specific regression model of the
form

$$
y_i \mid c_i = k \sim \mathcal{N}\left(x_i^\top \gamma_k,\;
d_{0i} + d_{1i}\mu_{ik} + d_{2i}\mu_{ik}^{2}\right),
\qquad
\mu_{ik} = x_i^\top \gamma_k.
$$

Gene modules are inferred through a Dirichlet process mixture.
Module-cellular-label associations are evaluated using both:

- a one-sided FGLS-based test with BH FDR control;
- a Bayesian posterior sign criterion reported through the package lFSR
  output.

An association is marked `dual_significant` only when both criteria pass
their configured thresholds.

For a module with no more genes than predictors, THREAD retains
descriptive coefficients but sets frequentist standard errors, P values,
and FDR values to `NA`. Such a module cannot be dual-significant.

## Documentation

After installation:

``` r
browseVignettes("THREAD")
```

The package includes:

- **Getting Started with THREAD**
- **Computing the LD-Corrected Gene-Level Response**
- **THREAD Input Formats and Workflow Entry Points**

The pkgdown site is available at:

``` text
https://jiacheng-lou.github.io/DPM/
```

## Reproducibility

For formal analyses, record at least:

- THREAD package version;
- Git commit or GitHub release tag;
- random seed;
- GWAS and expression input versions;
- genome build and LD ancestry;
- response diagnostics;
- `sessionInfo()`.

The included example data are checked by automated regression and
end-to-end tests.

## Citation

To obtain the current package citation, run:

``` r
citation("THREAD")
```

## License

THREAD is released under GPL-3.

## Issues and contact

Please report software problems through GitHub Issues:

``` text
https://github.com/Jiacheng-Lou/DPM/issues
```

For scientific questions, contact:

``` text
jiachenglou82@gmail.com
```
