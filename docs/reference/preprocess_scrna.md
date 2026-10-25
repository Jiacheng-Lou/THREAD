# Preprocess single-cell expression into a THREAD expression matrix

Converts a Seurat object or a gene-by-cell expression matrix into a
gene-by-label pseudobulk matrix. Labels can be cell types, subtypes,
states, clusters, domains, or any user-defined grouping variable. The
returned matrix has genes in rows and labels in columns.

The default workflow is:

1.  obtain normalized expression. If `already_normalized = FALSE`,
    Seurat's `NormalizeData()` is used with `LogNormalize`; if
    `already_normalized = TRUE`, the input expression is used directly;

2.  average expression over all cells sharing the same label;

3.  optionally apply a second log transform to the pseudobulk means.

This function is intentionally not a full scRNA-seq preprocessing
pipeline. Users should perform cell filtering, doublet removal, batch
correction and cell annotation upstream when needed. Users who already
have a gene-by-label matrix can pass it directly to
[`prepare_input()`](https://jiacheng-lou.github.io/THREAD/reference/prepare_input.md)
and skip this function.

## Usage

``` r
preprocess_scrna(
  object,
  meta = NULL,
  subtype_col,
  cell_col = NULL,
  assay = NULL,
  normalized_layer = "data",
  counts_layer = "counts",
  already_normalized = FALSE,
  normalization_method = "LogNormalize",
  scale_factor = 10000,
  second_log = TRUE,
  second_log_base = 2,
  pseudocount = 1,
  min_cells_per_subtype = 1,
  min_cells_per_gene = 0,
  drop_na_subtype = TRUE,
  return_sparse = FALSE,
  seurat_verbose = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object or a numeric matrix-like object with genes in rows and
  cells in columns. Sparse matrices from the Matrix package are
  supported for matrix input.

- meta:

  A data frame containing cell metadata. Required for matrix input;
  ignored for Seurat input unless supplied, in which case it overrides
  the object's metadata. Row names should match `colnames(object)`,
  unless `cell_col` is supplied.

- subtype_col:

  A single column name in the metadata containing the label to aggregate
  by, or a vector of labels with length equal to the number of cells.

- cell_col:

  Optional column name in `meta` containing cell IDs.

- assay:

  Seurat assay to use. If `NULL`, the object's default assay is used.
  Ignored for matrix input.

- normalized_layer:

  Assay layer/slot containing normalized expression when
  `already_normalized = TRUE`, and the layer/slot read after Seurat
  normalization when `already_normalized = FALSE`. Default `"data"`.

- counts_layer:

  Assay layer/slot containing raw counts for Seurat input when
  `already_normalized = FALSE`. Default `"counts"`.

- already_normalized:

  Logical. If `FALSE`, expression is treated as raw counts and
  normalized with Seurat's `NormalizeData()`. If `TRUE`, expression is
  treated as already normalized non-negative expression and is
  aggregated directly.

- normalization_method:

  Method passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html).
  Default `"LogNormalize"`.

- scale_factor:

  Scale factor passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html).
  Default `1e4`.

- second_log:

  Logical. Whether to apply a second log transform after label-level
  averaging. Default `TRUE`.

- second_log_base:

  Logarithm base for the second log transform. Supported values are `2`
  and `exp(1)`. Default `2`.

- pseudocount:

  Non-negative pseudocount used in the second log transform. Default
  `1`.

- min_cells_per_subtype:

  Minimum number of cells required for a label to be retained. Default
  `1`.

- min_cells_per_gene:

  Optional minimum number of cells in which a gene must have non-zero
  expression before aggregation. Default `0` disables this filter.

- drop_na_subtype:

  Logical. Whether to remove cells with missing labels. Default `TRUE`.

- return_sparse:

  Logical. Whether to return a sparse matrix when possible. The default
  `FALSE` returns a dense matrix because the number of label columns is
  usually small.

- seurat_verbose:

  Logical. Verbosity passed to Seurat normalization. Default `FALSE`.

- verbose:

  Logical. Whether to print concise progress messages. Default `TRUE`.

## Value

A gene-by-label numeric matrix suitable for
[`prepare_input()`](https://jiacheng-lou.github.io/THREAD/reference/prepare_input.md).

## Examples

``` r
if (FALSE) { # \dontrun{
X <- preprocess_scrna(
  object = seurat_obj,
  subtype_col = "cell_subtype",
  already_normalized = FALSE
)
} # }
```
