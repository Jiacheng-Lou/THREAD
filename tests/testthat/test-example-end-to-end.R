test_that("example data complete the THREAD workflow", {
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

  response <- readRDS(
    file.path(
      example_dir,
      "example_response.rds"
    )
  )

  thread_data <- prepare_input(
    X = X,
    response = response,
    K_init = 3L,
    min_genes = min(10L, nrow(X)),
    seed = 2028L,
    verbose = FALSE
  )

  expect_s3_class(
    thread_data,
    "thread_data"
  )

  expect_equal(
    nrow(thread_data$x_train),
    nrow(X)
  )

  expect_equal(
    ncol(thread_data$x_train),
    ncol(X)
  )

  fit <- run_thread(
    thread_data = thread_data,
    alpha = 0.1,
    m = 2L,
    n_iter = 40L,
    burnin = 20L,
    thin = 1L,
    seed = 2028L,
    save_scales = FALSE,
    verbose = FALSE
  )

  expect_s3_class(
    fit,
    "thread_fit"
  )

  expect_length(
    fit$samples,
    20L
  )

  partition <- get_partition(
    fit = fit,
    min_cluster_size = 1L,
    mode_window = 1L,
    return_dissimilarity = FALSE,
    verbose = FALSE
  )

  expect_length(
    partition$assignments,
    nrow(X)
  )

  expect_gte(
    partition$K,
    1L
  )

  sig <- test_significance(
    fit = fit,
    partition = partition,
    thread_data = thread_data,
    verbose = FALSE
  )

  expect_true(
    all(
      c(
        "fgls_estimable",
        "gamma_fgls",
        "se_fgls",
        "pval_fgls",
        "fdr",
        "sig_fdr",
        "lfsr",
        "sig_lfsr",
        "dual_significant"
      ) %in% names(sig$table)
    )
  )

  non_estimable <- !sig$table$fgls_estimable

  if (any(non_estimable, na.rm = TRUE)) {
    expect_true(
      all(
        is.na(
          sig$table$se_fgls[
            non_estimable
          ]
        )
      )
    )

    expect_true(
      all(
        is.na(
          sig$table$pval_fgls[
            non_estimable
          ]
        )
      )
    )

    expect_true(
      all(
        is.na(
          sig$table$fdr[
            non_estimable
          ]
        )
      )
    )

    expect_true(
      all(
        !sig$table$sig_fdr[
          non_estimable
        ]
      )
    )

    expect_true(
      all(
        !sig$table$dual_significant[
          non_estimable
        ]
      )
    )
  }

  gene_results <- gene_score(
    partition = partition,
    coefficients = sig$coefficients,
    x_train = thread_data$x_train,
    top_n = 10L,
    verbose = FALSE
  )

  expect_true(
    all(
      c(
        "scores",
        "top_genes"
      ) %in% names(gene_results)
    )
  )

  expect_length(
    gene_results$scores,
    nrow(X)
  )

  cellular_results <- cell_type_score(sig)

  expect_s3_class(
    cellular_results,
    "data.frame"
  )
})