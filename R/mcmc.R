# mcmc.R
# -----------------------------------------------------------------------------
# Core MCMC sampler for the DPM model. The whole sampler lives in one file so
# that the main loop and every conditional update it calls can be read
# together.
#
#   run_thread()             public entry point; consumes a `thread_data` object
#                            and returns a `thread_fit` object
#
# Internal conditional updates (not exported):
#   initialize_parameters()  starting state (external labels or k-means)
#   getLogPostgamma_p()       conditional log-posterior of one gamma coordinate
#   update_gamma_slice()      stepping-out slice sampler for gamma (non-conjugate)
#   update_sigma2_gibbs()    Gibbs update of the global Lasso scale sigma^2
#   update_tau2()            update of the local Lasso scales tau^2
#   update_lambda2()         update of the Lasso rate lambda^2
#
# The cluster-label update (Neal's Algorithm 8) is implemented in C++ and lives
# in src/update_clusters.cpp; it is exposed to R as update_clusters_cpp_hetero()
# and called directly from run_thread().
# -----------------------------------------------------------------------------


#' Initialise the sampler state
#'
#' @param K_init Integer number of starting clusters.
#' @param p Integer number of cell subtypes (columns of \code{x_train}).
#' @param x_train Numeric gene-by-subtype matrix.
#' @param y_train Numeric response vector.
#' @param alpha0 Inverse-Gamma hyperparameter for \code{sigma^2} (shape)
#' @param beta0 Inverse-Gamma hyperparameter for \code{sigma^2} (rate)
#' @param r Gamma hyperparameter (shape) for \code{lambda^2}
#' @param delta Gamma hyperparameter (rate) for \code{lambda^2}
#' @param init_method Either "external" (use \code{z_init}) or "kmeans".
#' @param z_init Optional integer starting labels (required for "external").
#' @param adaptive_ridge Logical; scale the ridge penalty of the per-cluster
#'   warm-start regression to the data. Default \code{TRUE}.
#'
#' @return A list with elements \code{z}, \code{gamma}, \code{sigma2},
#'   \code{tau2}, \code{lambda2}.
#'
#' @importFrom stats kmeans rgamma rexp
#' @noRd
initialize_parameters <- function(K_init, p, x_train, y_train,
                                  alpha0, beta0, r, delta,
                                  init_method = "kmeans",
                                  z_init = NULL,
                                  adaptive_ridge = TRUE) {
  if (init_method == "kmeans") {
    km <- stats::kmeans(cbind(y_train, x_train), centers = K_init, nstart = 20)
    z <- km$cluster
  } else if (init_method == "external") {
    if (is.null(z_init)) stop("'z_init' must be provided when init_method = 'external'.")
    z <- as.integer(z_init)
    if (any(z < 1 | z > K_init)) {
      stop("Labels in 'z_init' must lie between 1 and K_init.")
    }
  } else {
    stop("Unknown init_method: choose 'kmeans' or 'external'.")
  }

  lambda2 <- stats::rgamma(K_init, shape = r, rate = delta)
  sigma2 <- 1 / stats::rgamma(K_init, shape = alpha0, rate = beta0)
  tau2 <- matrix(stats::rexp(K_init * p, rate = 0.5), nrow = K_init, ncol = p)

  gamma <- matrix(0, nrow = K_init, ncol = p)

  if (adaptive_ridge) {
    base_lambda_ridge <- sum(diag(t(x_train) %*% x_train)) / nrow(x_train)
    base_lambda_ridge <- max(base_lambda_ridge, 1e-4)
  } else {
    base_lambda_ridge <- 1.0
  }

  for (k in 1:K_init) {
    idx <- which(z == k)
    x_k <- x_train[idx, , drop = FALSE]
    y_k <- y_train[idx]

    if (length(y_k) > p) {
      gamma[k, ] <- tryCatch(
        {
          lambda_ridge <- if (adaptive_ridge) {
            base_lambda_ridge * (length(y_k) / nrow(x_train))
          } else {
            base_lambda_ridge
          }
          p_dim <- ncol(x_k)
          as.vector(solve(t(x_k) %*% x_k + lambda_ridge * diag(p_dim)) %*% t(x_k) %*% y_k)
        },
        error = function(e) rep(0, p)
      )
    } else {
      gamma[k, ] <- rep(0, p)
    }
  }

  rownames(gamma) <- as.character(1:K_init)
  rownames(tau2) <- as.character(1:K_init)
  names(sigma2) <- as.character(1:K_init)
  names(lambda2) <- as.character(1:K_init)

  list(z = z, gamma = gamma, sigma2 = sigma2, tau2 = tau2, lambda2 = lambda2)
}


#' Conditional log-posterior of a single gamma coordinate
#'
#' @param gamma_p_val Candidate value for coordinate \code{p}.
#' @param p Index of the coordinate being updated.
#' @param gamma_k_current Current coefficient vector of the cluster.
#' @param x_k,y_k Cluster design matrix and response.
#' @param d0_k,d1_k,d2_k Cluster heteroscedastic-variance coefficients.
#' @param sigma2_k,tau2_k Cluster Lasso scales.
#'
#' @return Scalar log-posterior (log-likelihood + log-prior of coordinate p).
#'
#' @importFrom stats dnorm
#' @noRd
getLogPostgamma_p <- function(gamma_p_val, p, gamma_k_current, x_k, y_k,
                             d0_k, d1_k, d2_k, sigma2_k, tau2_k) {
  gamma_temp <- gamma_k_current
  gamma_temp[p] <- gamma_p_val

  mu_k <- as.vector(x_k %*% gamma_temp)

  v_k <- d0_k + d1_k * mu_k + d2_k * (mu_k^2)
  v_k[v_k <= 1e-10] <- 1e-10 # strict positivity guard

  logLike <- sum(stats::dnorm(y_k, mean = mu_k, sd = sqrt(v_k), log = TRUE))
  logPrior <- -(gamma_p_val^2) / (2 * sigma2_k * tau2_k[p])

  logLike + logPrior
}


#' Stepping-out slice sampler for a cluster's coefficient vector
#'
#' @description
#' Updates \code{gamma_k} one coordinate at a time. The kernel is non-conjugate
#' because the heteroscedastic variance depends on \code{gamma_k} through both
#' the mean and the variance of the Gaussian likelihood.
#'
#' @param z Current label vector.
#' @param k Cluster being updated.
#' @param x_train,y_train Full design matrix and response.
#' @param d0,d1,d2 Full heteroscedastic-variance coefficient vectors.
#' @param gamma_k_current Current coefficient vector of cluster \code{k}.
#' @param sigma2_k,tau2_k Cluster Lasso scales.
#' @param adaptive_w Logical; set the initial slice width from the prior scale
#'   rather than using a fixed width. Default \code{TRUE}.
#'
#' @return Updated coefficient vector for cluster \code{k}.
#'
#' @importFrom stats rnorm rexp runif
#' @noRd
update_gamma_slice <- function(z, k, x_train, y_train, d0, d1, d2,
                              gamma_k_current, sigma2_k, tau2_k,
                              adaptive_w = TRUE) {
  is_k <- (z == k)
  x_k <- x_train[is_k, , drop = FALSE]
  y_k <- y_train[is_k]
  d0_k <- d0[is_k]
  d1_k <- d1[is_k]
  d2_k <- d2[is_k]

  p_dim <- ncol(x_k)

  # empty cluster: draw from the prior
  if (nrow(x_k) == 0) {
    return(stats::rnorm(p_dim, mean = 0, sd = sqrt(sigma2_k * tau2_k)))
  }

  gamma_new <- gamma_k_current
  E <- 40 # maximum number of stepping-out expansions

  # initial slice width: prior scale when adaptive, else a small fixed width
  w <- if (adaptive_w) max(sqrt(sigma2_k * mean(tau2_k)), 1e-6) else 0.001

  for (p in 1:p_dim) {
    x_current <- gamma_new[p]

    # slice height
    z_slice <- getLogPostgamma_p(
      x_current, p, gamma_new, x_k, y_k,
      d0_k, d1_k, d2_k, sigma2_k, tau2_k
    ) -
      stats::rexp(1, rate = 1)

    u <- stats::runif(1)
    L <- x_current - w * u
    R <- L + w

    # stepping out
    J <- floor(E * stats::runif(1))
    K_exp <- (E - 1) - J

    while (J > 0) {
      L <- L - w
      if (getLogPostgamma_p(
        L, p, gamma_new, x_k, y_k,
        d0_k, d1_k, d2_k, sigma2_k, tau2_k
      ) <= z_slice) {
        break
      }
      J <- J - 1
    }
    while (K_exp > 0) {
      R <- R + w
      if (getLogPostgamma_p(
        R, p, gamma_new, x_k, y_k,
        d0_k, d1_k, d2_k, sigma2_k, tau2_k
      ) <= z_slice) {
        break
      }
      K_exp <- K_exp - 1
    }

    # propose and shrink
    repeat {
      newX <- L + stats::runif(1) * (R - L)
      currentLogPost <- getLogPostgamma_p(
        newX, p, gamma_new, x_k, y_k,
        d0_k, d1_k, d2_k, sigma2_k, tau2_k
      )
      if (currentLogPost >= z_slice) break
      if (newX < x_current) L <- newX else R <- newX
    }

    gamma_new[p] <- newX
  }

  gamma_new
}


#' Gibbs update of the global Lasso scale sigma^2 (Inverse-Gamma)
#'
#' @param gamma_k Cluster coefficient vector.
#' @param tau2_k Cluster local scales.
#' @param alpha0 Inverse-Gamma hyperparameter for \code{sigma^2} (shape)
#' @param beta0 Inverse-Gamma hyperparameter for \code{sigma^2} (rate)
#'
#' @return Updated \code{sigma^2} (scalar, lower-bounded for stability).
#'
#' @importFrom stats rgamma
#' @noRd
update_sigma2_gibbs <- function(gamma_k, tau2_k, alpha0, beta0) {
  p <- length(gamma_k)
  shape_post <- alpha0 + p / 2
  rate_post <- beta0 + 0.5 * sum((gamma_k^2) / tau2_k)
  sigma2_new <- 1 / stats::rgamma(1, shape = shape_post, rate = rate_post)
  max(sigma2_new, 1e-8)
}


#' Update the local Lasso scales tau^2 (Inverse-Gaussian for 1/tau^2)
#'
#' @param gamma_k Cluster coefficient vector.
#' @param lambda2_k Cluster Lasso rate.
#' @param sigma2_k Cluster global scale.
#'
#' @return Updated vector of \code{tau^2} (lower-bounded for stability).
#'
#' @importFrom statmod rinvgauss
#' @noRd
update_tau2 <- function(gamma_k, lambda2_k, sigma2_k) {
  p <- length(gamma_k)
  # 1/tau^2 ~ InvGauss(mean = sqrt(lambda^2 sigma^2 / gamma^2), shape = lambda^2)
  mu_k <- sqrt(lambda2_k * sigma2_k / (gamma_k^2 + 1e-10))
  shape_k <- lambda2_k

  inv_tau2 <- tryCatch(
    statmod::rinvgauss(n = p, mean = mu_k, shape = shape_k),
    # when gamma is tiny, mu_k -> Inf and the sampler may fail: 1/tau^2 -> large
    error = function(e) rep(1e10, p)
  )

  pmax(1 / inv_tau2, 1e-10)
}


#' Gibbs update of the Lasso rate lambda^2 (Gamma)
#'
#' @param tau2_k Cluster local scales.
#' @param r,delta Gamma hyperparameters.
#'
#' @return Updated \code{lambda^2} (scalar).
#'
#' @importFrom stats rgamma
#' @noRd
update_lambda2 <- function(tau2_k, r, delta) {
  p <- length(tau2_k)
  shape <- max(r + p, 1e-6)
  rate <- max(delta + sum(tau2_k) / 2, 1e-10)
  stats::rgamma(1, shape = shape, rate = rate)
}


#' Fit the THREAD model by MCMC
#'
#' @description
#' Runs the Dirichlet-process-mixture sampler on a prepared \code{thread_data}
#' object and returns a \code{thread_fit} object holding the posterior draws,
#' diagnostics and settings. Cluster labels are updated by Neal's Algorithm 8
#' (auxiliary components proposed from the prior); the coefficient vectors are
#' updated by stepping-out slice sampling; the Bayesian-Lasso scale parameters
#' have closed-form conjugate updates.
#'
#' All hyperparameters are explicit arguments or are read from
#' \code{thread_data$priors}; the function holds no hidden state. Posterior draws
#' are stored separately from the input data because they can be large; the
#' downstream functions therefore take both the \code{thread_fit} and the
#' \code{thread_data}.
#'
#' @param thread_data A list produced by \code{prepare_input()} containing at least
#'   \code{x_train} (gene-by-subtype matrix), \code{y_train} (response vector),
#'   \code{d0}, \code{d1}, \code{d2} (heteroscedastic-variance coefficients) and,
#'   optionally, \code{z_init} (starting labels), \code{gene_names} and
#'   \code{priors}.
#' @param alpha Dirichlet-process concentration parameter. Default \code{0.1}.
#' @param m Number of auxiliary components in Neal's Algorithm 8. Default
#'   \code{3}.
#' @param K_init Number of starting clusters when \code{z_init} is not supplied.
#'   Default \code{10}. When \code{z_init} is supplied this is taken from the
#'   labels.
#' @param n_iter Total number of MCMC iterations. Default \code{30000}.
#' @param burnin Number of initial iterations discarded. Default \code{25000}.
#' @param thin Keep one draw every \code{thin} post-burn-in iterations. Default
#'   set to \code{1} to keep every draw (larger output).
#' @param priors Optional list overriding \code{thread_data$priors}; must contain
#'   \code{alpha0}, \code{beta0}, \code{r}, \code{delta}.
#' @param seed Optional integer random seed for reproducibility. Default
#'   \code{NULL}.
#' @param adaptive_ridge Logical; passed to the warm-start initialisation.
#'   Default \code{TRUE}.
#' @param adaptive_w Logical; passed to the slice sampler. Default \code{TRUE}.
#' @param verbose Logical; show a progress bar and status messages. Default
#'   \code{TRUE}.
#' @param save_scales Logical. Whether to retain tau2 and lambda2 in every
#'   posterior draw. Default FALSE because these objects can substantially
#'   increase memory and storage use.
#'
#' @return An object of class \code{thread_fit}: a list with
#'   \describe{
#'     \item{samples}{list with one element per retained draw, each a list of
#'       \code{z} (labels), \code{gamma} (K-by-p coefficient matrix),
#'       \code{sigma2} and \code{K}.}
#'     \item{diagnostics}{list with the per-draw cluster count
#'       (\code{cluster_counts}) and the wall-clock \code{runtime}.}
#'     \item{settings}{the hyperparameters and run controls used.}
#'     \item{dims}{\code{n_genes}, \code{n_subtypes}, \code{gene_names},
#'       \code{subtype_names}.}
#'   }
#'
#' @examples
#' \dontrun{
#' fit <- run_thread(thread_data, n_iter = 30000, burnin = 25000, thin = 1, seed = 123)
#' table(fit$diagnostics$cluster_counts)
#' }
#'
#' @importFrom stats rgamma rnorm rexp
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
run_thread <- function(thread_data,
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
                    verbose = TRUE) {
  # ---- local validation helpers ---------------------------------------------
  validate_integer_control <- function(x,
                                       name,
                                       minimum = 0L) {
    if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x != floor(x) ||
      x < minimum) {
      stop(
        "'", name,
        "' must be an integer >= ",
        minimum,
        "."
      )
    }

    as.integer(x)
  }

  validate_scalar_logical <- function(x, name) {
    if (!is.logical(x) ||
      length(x) != 1L ||
      is.na(x)) {
      stop(
        "'", name,
        "' must be TRUE or FALSE."
      )
    }

    invisible(TRUE)
  }

  validate_positive_scalar <- function(x, name) {
    if (!is.numeric(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !is.finite(x) ||
      x <= 0) {
      stop(
        "'", name,
        "' must be a finite, strictly positive scalar."
      )
    }

    as.numeric(x)
  }

  # ---- validate control arguments ------------------------------------------
  alpha <- validate_positive_scalar(alpha, "alpha")
  m <- validate_integer_control(m, "m", minimum = 1L)
  K_init <- validate_integer_control(
    K_init,
    "K_init",
    minimum = 1L
  )
  n_iter <- validate_integer_control(
    n_iter,
    "n_iter",
    minimum = 1L
  )
  burnin <- validate_integer_control(
    burnin,
    "burnin",
    minimum = 0L
  )
  thin <- validate_integer_control(
    thin,
    "thin",
    minimum = 1L
  )

  validate_scalar_logical(
    adaptive_ridge,
    "adaptive_ridge"
  )
  validate_scalar_logical(
    adaptive_w,
    "adaptive_w"
  )
  validate_scalar_logical(
    save_scales,
    "save_scales"
  )
  validate_scalar_logical(
    verbose,
    "verbose"
  )

  if (!is.null(seed)) {
    seed <- validate_integer_control(
      seed,
      "seed",
      minimum = 0L
    )
  }

  if (burnin >= n_iter) {
    stop(
      "'burnin' must be strictly less than 'n_iter'."
    )
  }

  if ((n_iter - burnin) < thin) {
    stop(
      "The requested burnin/thin settings retain no posterior draws."
    )
  }

  # ---- validate thread_data -------------------------------------------------
  if (!is.list(thread_data)) {
    stop(
      "'thread_data' must be a list produced by prepare_input()."
    )
  }

  required_fields <- c(
    "x_train",
    "y_train",
    "d0",
    "d1",
    "d2"
  )

  missing_fields <- setdiff(
    required_fields,
    names(thread_data)
  )

  if (length(missing_fields) > 0L) {
    stop(
      "'thread_data' is missing field(s): ",
      paste(missing_fields, collapse = ", "),
      ". Build it with prepare_input()."
    )
  }

  x_train <- as.matrix(thread_data$x_train)
  storage.mode(x_train) <- "double"

  y_train <- as.numeric(thread_data$y_train)
  d0 <- as.numeric(thread_data$d0)
  d1 <- as.numeric(thread_data$d1)
  d2 <- as.numeric(thread_data$d2)

  if (nrow(x_train) == 0L ||
    ncol(x_train) == 0L) {
    stop(
      "'x_train' must have at least one row and one column."
    )
  }

  N <- nrow(x_train)
  p_dim <- ncol(x_train)

  response_lengths <- c(
    y_train = length(y_train),
    d0 = length(d0),
    d1 = length(d1),
    d2 = length(d2)
  )

  if (any(response_lengths != N)) {
    stop(
      "Dimension mismatch: x_train has ",
      N,
      " rows, but y_train/d0/d1/d2 have lengths ",
      paste(response_lengths, collapse = "/"),
      "."
    )
  }

  if (any(!is.finite(x_train))) {
    stop(
      "'x_train' contains non-finite values."
    )
  }

  if (any(!is.finite(y_train))) {
    stop(
      "'y_train' contains non-finite values."
    )
  }

  if (any(!is.finite(d0)) ||
    any(!is.finite(d1)) ||
    any(!is.finite(d2))) {
    stop(
      "'d0', 'd1' and 'd2' must contain only finite values."
    )
  }

  if (any(d0 <= 0)) {
    stop(
      "'d0' must be strictly positive for every gene."
    )
  }

  if (any(d1 < 0) || any(d2 < 0)) {
    stop(
      "'d1' and 'd2' must be non-negative for every gene."
    )
  }

  # ---- validate priors ------------------------------------------------------
  if (is.null(priors)) {
    priors <- thread_data$priors
  }

  if (is.null(priors)) {
    stop(
      "No priors found: pass 'priors' or include them in ",
      "'thread_data$priors'."
    )
  }

  needed_priors <- c(
    "alpha0",
    "beta0",
    "r",
    "delta"
  )

  missing_priors <- setdiff(
    needed_priors,
    names(priors)
  )

  if (length(missing_priors) > 0L) {
    stop(
      "'priors' is missing: ",
      paste(missing_priors, collapse = ", ")
    )
  }

  alpha0 <- validate_positive_scalar(
    priors$alpha0,
    "priors$alpha0"
  )
  beta0 <- validate_positive_scalar(
    priors$beta0,
    "priors$beta0"
  )
  r <- validate_positive_scalar(
    priors$r,
    "priors$r"
  )
  delta <- validate_positive_scalar(
    priors$delta,
    "priors$delta"
  )

  # ---- starting labels ------------------------------------------------------
  z_init <- thread_data$z_init

  if (!is.null(z_init)) {
    if (length(z_init) != N) {
      stop(
        "'z_init' length (",
        length(z_init),
        ") does not match the number of genes (",
        N,
        ")."
      )
    }

    if (anyNA(z_init) ||
      any(!is.finite(as.numeric(z_init)))) {
      stop(
        "'z_init' contains missing or non-finite labels."
      )
    }

    z_init <- match(
      z_init,
      sort(unique(z_init))
    )
    z_init <- as.integer(z_init)

    K_init_used <- max(z_init)

    if (verbose && K_init_used != K_init) {
      THREAD_log(
        "MCMC",
        sprintf(
          "Using K_init = %d from supplied z_init (ignoring K_init = %d).",
          K_init_used,
          K_init
        )
      )
    }

    K_init <- K_init_used
    init_method <- "external"
  } else {
    if (K_init > N) {
      stop(
        "'K_init' cannot exceed the number of genes when ",
        "z_init is not supplied."
      )
    }

    init_method <- "kmeans"
  }

  # ---- metadata -------------------------------------------------------------
  gene_names <- thread_data$gene_names

  if (is.null(gene_names)) {
    gene_names <- rownames(x_train)
  }

  if (is.null(gene_names)) {
    gene_names <- paste0(
      "Gene",
      seq_len(N)
    )
  }

  if (length(gene_names) != N) {
    stop(
      "The number of gene names does not match nrow(x_train)."
    )
  }

  subtype_names <- colnames(x_train)

  if (is.null(subtype_names)) {
    subtype_names <- paste0(
      "Subtype",
      seq_len(p_dim)
    )
    colnames(x_train) <- subtype_names
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ---- auxiliary proposal scale --------------------------------------------
  sigma2_prior_mode <- beta0 / (alpha0 + 1)
  tau2_prior_mean <- 2 / (r / delta)

  prior_sd_proposal <- sqrt(
    sigma2_prior_mode * tau2_prior_mean
  )

  if (!is.finite(prior_sd_proposal) ||
    prior_sd_proposal <= 0) {
    stop(
      "The auxiliary-cluster proposal scale is not finite and positive."
    )
  }

  total_save <- floor(
    (n_iter - burnin) / thin
  )

  if (total_save < 1L) {
    stop(
      "No posterior draws would be retained with the supplied controls."
    )
  }

  if (verbose) {
    THREAD_log(
      "MCMC",
      sprintf(
        "%d genes x %d subtypes | alpha = %g, m = %d, K_init = %d",
        N,
        p_dim,
        alpha,
        m,
        K_init
      )
    )

    THREAD_log(
      "MCMC",
      sprintf(
        "iterations = %d, burnin = %d, thin = %d, expected retained draws = %d",
        n_iter,
        burnin,
        thin,
        total_save
      )
    )

    THREAD_log(
      "MCMC",
      sprintf("auxiliary-cluster proposal prior_sd = %.6g", prior_sd_proposal)
    )

    THREAD_log("MCMC", "save_scales = ", save_scales)
  }

  # ---- initialise sampler ---------------------------------------------------
  current_state <- initialize_parameters(
    K_init = K_init,
    p = p_dim,
    x_train = x_train,
    y_train = y_train,
    alpha0 = alpha0,
    beta0 = beta0,
    r = r,
    delta = delta,
    init_method = init_method,
    z_init = z_init,
    adaptive_ridge = adaptive_ridge
  )

  if (any(!is.finite(current_state$gamma)) ||
    any(!is.finite(current_state$sigma2)) ||
    any(!is.finite(current_state$tau2)) ||
    any(!is.finite(current_state$lambda2))) {
    stop(
      "Sampler initialisation produced non-finite parameter values."
    )
  }

  samples <- vector(
    "list",
    total_save
  )
  cluster_counts <- integer(
    total_save
  )
  save_idx <- 1L

  pb <- NULL

  if (verbose) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = n_iter,
      style = 3
    )

    on.exit(
      {
        if (!is.null(pb)) {
          try(close(pb), silent = TRUE)
        }
      },
      add = TRUE
    )
  }

  t_start <- Sys.time()

  # ---- MCMC loop ------------------------------------------------------------
  for (iter in seq_len(n_iter)) {
    # 1. Update cluster labels using Neal's Algorithm 8.
    res <- update_clusters_cpp_hetero(
      x_train = x_train,
      y_train = y_train,
      d0 = d0,
      d1 = d1,
      d2 = d2,
      z = current_state$z,
      gamma = current_state$gamma,
      alpha = alpha,
      m = m,
      prior_sd = prior_sd_proposal
    )

    z_new <- as.integer(res$z)
    gamma_new <- as.matrix(res$gamma)
    cluster_map <- as.integer(
      res$old_cluster_map
    )

    K_new <- nrow(gamma_new)

    if (K_new < 1L ||
      ncol(gamma_new) != p_dim) {
      stop(
        "Internal error: invalid gamma dimensions returned ",
        "by update_clusters_cpp_hetero()."
      )
    }

    if (length(z_new) != N ||
      anyNA(z_new) ||
      any(z_new < 1L) ||
      any(z_new > K_new)) {
      stop(
        "Internal error: invalid labels returned by ",
        "update_clusters_cpp_hetero()."
      )
    }

    if (length(cluster_map) != K_new) {
      stop(
        "Internal error: length(old_cluster_map) does not ",
        "match nrow(gamma_new)."
      )
    }

    if (any(!is.finite(gamma_new))) {
      stop(
        "Non-finite gamma values were returned by the cluster update."
      )
    }

    old_K <- nrow(current_state$gamma)

    if (any(cluster_map > old_K)) {
      stop(
        "Internal error: old_cluster_map contains an invalid cluster index."
      )
    }

    # 2. Reconstruct Lasso scales after relabelling.
    sigma2_new <- numeric(K_new)
    lambda2_new <- numeric(K_new)
    tau2_new <- matrix(
      0,
      nrow = K_new,
      ncol = p_dim
    )

    for (k in seq_len(K_new)) {
      old_idx <- cluster_map[k]

      if (old_idx > 0L) {
        old_name <- as.character(old_idx)

        sigma2_new[k] <-
          current_state$sigma2[old_name]
        lambda2_new[k] <-
          current_state$lambda2[old_name]
        tau2_new[k, ] <-
          current_state$tau2[old_name, ]
      } else {
        sigma2_new[k] <- 1 / stats::rgamma(
          1L,
          shape = alpha0,
          rate = beta0
        )

        lambda2_new[k] <- stats::rgamma(
          1L,
          shape = r,
          rate = delta
        )

        tau_rate <- max(
          lambda2_new[k] / 2,
          1e-10
        )

        tau2_new[k, ] <- stats::rexp(
          p_dim,
          rate = tau_rate
        )

        tau2_new[k, ] <- pmax(
          tau2_new[k, ],
          1e-10
        )
      }
    }

    if (any(!is.finite(sigma2_new)) ||
      any(!is.finite(lambda2_new)) ||
      any(!is.finite(tau2_new))) {
      stop(
        "Non-finite Lasso scale values were produced during ",
        "cluster reconstruction."
      )
    }

    rownames(gamma_new) <-
      as.character(seq_len(K_new))
    rownames(tau2_new) <-
      as.character(seq_len(K_new))
    names(sigma2_new) <-
      as.character(seq_len(K_new))
    names(lambda2_new) <-
      as.character(seq_len(K_new))

    current_state$z <- z_new
    current_state$gamma <- gamma_new
    current_state$sigma2 <- sigma2_new
    current_state$tau2 <- tau2_new
    current_state$lambda2 <- lambda2_new

    # 3. Update module-specific parameters.
    K <- nrow(current_state$gamma)

    for (k in seq_len(K)) {
      current_state$lambda2[k] <- update_lambda2(
        current_state$tau2[k, ],
        r,
        delta
      )

      current_state$tau2[k, ] <- update_tau2(
        current_state$gamma[k, ],
        current_state$lambda2[k],
        current_state$sigma2[k]
      )

      current_state$sigma2[k] <- update_sigma2_gibbs(
        current_state$gamma[k, ],
        current_state$tau2[k, ],
        alpha0,
        beta0
      )

      current_state$gamma[k, ] <- update_gamma_slice(
        current_state$z,
        k,
        x_train,
        y_train,
        d0,
        d1,
        d2,
        current_state$gamma[k, ],
        current_state$sigma2[k],
        current_state$tau2[k, ],
        adaptive_w = adaptive_w
      )
    }

    if (any(!is.finite(current_state$gamma)) ||
      any(!is.finite(current_state$sigma2)) ||
      any(!is.finite(current_state$tau2)) ||
      any(!is.finite(current_state$lambda2))) {
      stop(
        "Sampler produced non-finite parameter values at iteration ",
        iter,
        "."
      )
    }

    # 4. Store a thinned posterior draw.
    if (iter > burnin &&
      ((iter - burnin) %% thin == 0L) &&
      save_idx <= total_save) {
      draw <- list(
        z = current_state$z,
        gamma = current_state$gamma,
        sigma2 = current_state$sigma2,
        K = K
      )

      if (save_scales) {
        draw$tau2 <- current_state$tau2
        draw$lambda2 <- current_state$lambda2
      }

      samples[[save_idx]] <- draw
      cluster_counts[save_idx] <- K
      save_idx <- save_idx + 1L
    }

    if (verbose) {
      utils::setTxtProgressBar(
        pb,
        iter
      )
    }
  }

  runtime <- Sys.time() - t_start

  if (verbose && !is.null(pb)) {
    close(pb)
    pb <- NULL
  }

  n_saved <- save_idx - 1L

  if (n_saved < 1L) {
    stop(
      "MCMC completed but no posterior draws were retained."
    )
  }

  if (n_saved != total_save) {
    samples <- samples[seq_len(n_saved)]
    cluster_counts <- cluster_counts[seq_len(n_saved)]

    warning(
      "Expected ", total_save,
      " retained draws but stored ", n_saved,
      "."
    )
  }

  if (verbose) {
    THREAD_log(
      "MCMC",
      sprintf(
        "MCMC complete in %s. Retained %d draws.",
        format(round(runtime, 2)),
        n_saved
      )
    )

    THREAD_log(
      "MCMC",
      sprintf(
        "Occupied modules across retained draws: min = %d, median = %.1f, max = %d.",
        min(cluster_counts),
        stats::median(cluster_counts),
        max(cluster_counts)
      )
    )
  }

  structure(
    list(
      samples = samples,
      diagnostics = list(
        cluster_counts = cluster_counts,
        runtime = runtime,
        n_saved = n_saved
      ),
      settings = list(
        alpha0 = alpha0,
        beta0 = beta0,
        r = r,
        delta = delta,
        alpha = alpha,
        m = m,
        K_init = K_init,
        n_iter = n_iter,
        burnin = burnin,
        thin = thin,
        seed = seed,
        adaptive_ridge = adaptive_ridge,
        adaptive_w = adaptive_w,
        save_scales = save_scales
      ),
      dims = list(
        n_genes = N,
        n_subtypes = p_dim,
        gene_names = gene_names,
        subtype_names = subtype_names
      )
    ),
    class = "thread_fit"
  )
}
