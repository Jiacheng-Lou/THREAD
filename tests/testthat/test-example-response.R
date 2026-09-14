test_that("example LD response is reproducible", {
  example_dir <- system.file(
    "extdata",
    "example",
    package = "THREAD",
    mustWork = TRUE
  )

  X <- readRDS(
    file.path(
      example_dir,
      "example_x.rds"
    )
  )

  expected <- readRDS(
    file.path(
      example_dir,
      "example_response.rds"
    )
  )

  observed <- compute_response(
    gene_snp_matches = file.path(
      example_dir,
      "example_gene_snp_matches.rds"
    ),
    gwas_clean = file.path(
      example_dir,
      "example_gwas.tsv.gz"
    ),
    ld_dir = example_dir,
    ld_file_template = "example_ld_chr{chr}.rds",
    match_mode = "snp_name",
    genes = rownames(X),
    chromosomes = "1",
    n_cores = 1L,
    include_diagnostics = TRUE,
    drop_failed = FALSE,
    verbose = FALSE
  )

  observed <- observed[
    order(observed$gene_name), ,
    drop = FALSE
  ]

  expected <- expected[
    order(expected$gene_name), ,
    drop = FALSE
  ]

  expect_equal(
    observed$gene_name,
    expected$gene_name
  )

  expect_equal(
    observed$n_snps,
    expected$n_snps
  )

  expect_equal(
    observed$status,
    expected$status
  )

  expect_equal(
    observed[, c(
      "d0",
      "d1",
      "d2",
      "final_y"
    )],
    expected[, c(
      "d0",
      "d1",
      "d2",
      "final_y"
    )],
    tolerance = 1e-8
  )

  expect_true(all(observed$status == "ok"))
  expect_true(all(is.finite(observed$final_y)))
  expect_true(all(observed$d0 > 0))
  expect_true(all(observed$d1 >= 0))
  expect_true(all(observed$d2 >= 0))

  if ("ld_max_asymmetry" %in% names(observed)) {
    expect_true(
      all(
        observed$ld_max_asymmetry < 1e-12,
        na.rm = TRUE
      )
    )
  }
})