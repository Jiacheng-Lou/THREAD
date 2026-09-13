test_that("installed THREAD example files are complete and unchanged", {
  example_dir <- system.file(
    "extdata",
    "example",
    package = "THREAD",
    mustWork = TRUE
  )

  expected_files <- c(
    "example_x.rds",
    "example_response.rds",
    "example_gwas.tsv.gz",
    "example_gene_snp_matches.rds",
    "example_ld_chr1.rds",
    "example_data_manifest.csv",
    "example_data_README.txt"
  )

  expected_paths <- file.path(
    example_dir,
    expected_files
  )

  expect_true(all(file.exists(expected_paths)))

  manifest_path <- file.path(
    example_dir,
    "example_data_manifest.csv"
  )

  manifest <- data.table::fread(
    manifest_path,
    showProgress = FALSE
  )

  expect_true(
    all(
      c(
        "file",
        "size_bytes",
        "md5"
      ) %in% names(manifest)
    )
  )

  manifest_paths <- file.path(
    example_dir,
    manifest$file
  )

  expect_true(all(file.exists(manifest_paths)))

  observed_size <- unname(
    file.info(manifest_paths)$size
  )

  observed_md5 <- unname(
    tools::md5sum(manifest_paths)
  )

  expect_equal(
    observed_size,
    manifest$size_bytes
  )

  expect_equal(
    as.character(observed_md5),
    as.character(manifest$md5)
  )
})