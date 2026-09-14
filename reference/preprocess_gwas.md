# Quality-control and harmonise GWAS summary statistics

Applies a standard set of GWAS quality-control filters and returns a
clean, uniformly named table that the downstream THREAD functions
consume. The input may use any column names; supply `col_map` to declare
which of your columns correspond to the required fields. A user who
already holds a clean table with the standard column names can pass it
straight to
[`match_snps_to_genes`](https://jiacheng-lou.github.io/THREAD/reference/match_snps_to_genes.md)
and
[`compute_response`](https://jiacheng-lou.github.io/THREAD/reference/compute_response.md)
and skip this step.

## Usage

``` r
preprocess_gwas(
  sumstats,
  col_map = NULL,
  remove_sex_chrom = TRUE,
  remove_mhc = TRUE,
  mhc_range = c(2.5e+07, 3.4e+07),
  maf_min = 0.01,
  info_min = 0.8,
  n_min = 0,
  beta_max = 10,
  verbose = TRUE
)
```

## Arguments

- sumstats:

  A data.frame / data.table of summary statistics, or a single string
  giving the path to a delimited file (read with
  [`data.table::fread`](https://rdrr.io/pkg/data.table/man/fread.html)).

- col_map:

  Optional named character vector mapping standard field names to the
  column names used in `sumstats`. Recognised standard names are
  `chrom`, `pos`, `rsid`, `beta`, `se` (required) and `maf`, `freq`,
  `info`, `pval`, `n` (optional; each enables its corresponding filter).
  Example:
  `c(chrom = "CHR", pos = "POS", rsid = "SNP", beta = "BETA", se = "SE")`.
  If `NULL`, columns are assumed to already use the standard names.

- remove_sex_chrom:

  Logical; keep only autosomes 1-22. Default `TRUE`.

- remove_mhc:

  Logical; drop the MHC region on chromosome 6. Default `TRUE`.

- mhc_range:

  Integer length-2 giving the MHC start/end in base pairs (build
  hg19/GRCh37). Default `c(25e6, 34e6)`.

- maf_min:

  Minimum minor-allele frequency; applied only when a `maf` or `freq`
  column is available. Default `0.01`.

- info_min:

  Minimum imputation INFO score; applied only when an `info` column is
  available. Default `0.8`.

- n_min:

  Minimum per-SNP sample size; applied only when an `n` column is
  available and `n_min > 0`. Default `0` (no filter).

- beta_max:

  Maximum absolute effect size; SNPs with `abs(beta)` above this are
  dropped as likely artefacts. Default `10`. Set to `Inf` to disable.

- verbose:

  Logical; print a per-step row-count report. Default `TRUE`.

## Value

A `data.table` with columns `chrom`, `pos`, `rsid`, `beta`, `se`, and
`maf` when available, sorted by `chrom` then `pos`.

## Examples

``` r
if (FALSE) { # \dontrun{
clean <- preprocess_gwas(
  "PGC_BIP.tsv",
  col_map = c(
    chrom = "CHROM", pos = "POS", rsid = "ID",
    beta = "BETA", se = "SE", info = "IMPINFO"
  )
)
} # }
```
