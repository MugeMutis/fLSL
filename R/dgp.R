# ================================================================
# Data-generating processes of the Monte Carlo experiments (Section 4)
#
# design = "main" (Experiments 1-3, eq. 4.1):
#   X_i(t) = sum_{k=1}^{nphi} sqrt(2) k^{-3/2} {a_ik cos(k pi t) + b_ik sin(k pi t)}
#   beta(t) = beta_scale * sin(2 pi t)            (paper: beta_scale = 2)
#   m_i = alpha + int X_i(t) beta(t) dt           (trapezoidal rule)
# design = "exact" (Experiment 2, exact design):
#   X_i(t) = sum_{k=1}^4 lambda_k^{1/2} zeta_ik phi_k(t), phi_k = sqrt(2) sin(k pi t),
#   beta = c (phi_1 - phi_2 + phi_3 - phi_4), c such that sd(m_i) = 1.5,
#   m_i computed exactly from the scores.
# Both:
#   locations uniform on [0,1]^2, W row-standardized inverse-distance
#   weights over the k_nn = 4 nearest neighbours,
#   p_i = Lambda( (S m)_i / s_i ),  S = (I - rho W)^{-1},  s_i^2 = (S S')_ii,
#   Y_i | X, W ~ Bernoulli(p_i), conditionally independent.
# dependence_range > 0 (Experiment 2, Part C): the responses are instead
#   Y_i = 1{Phi(e_i) <= p_i}, e ~ N(0, C), C_ij = exp(-d_ij / range),
#   so that P(Y_i = 1 | X, W) = p_i while nearby responses are dependent.
#
# The random-number streams follow the experiment scripts of the paper, so
# that a given seed reproduces the same data set.
# W_format = "sparse" returns W as a sparse Matrix (used for large n); the
# k-nearest-neighbour search is done in row blocks so that the full
# n x n distance matrix is never stored. When rho = 0 no n x n inverse is
# computed.
# ================================================================

.trap_weights_dgp <- function(gp) {
  h <- diff(gp)
  c(h[1L] / 2, if (length(h) > 1L) (utils::head(h, -1L) + utils::tail(h, -1L)) / 2 else numeric(0),
    utils::tail(h, 1L) / 2)
}

.knn_W_dgp <- function(coords, k_nn = 4L, power = 1, sparse = FALSE, chunk = 1000L) {
  n <- nrow(coords)
  if (length(k_nn) != 1L || !is.finite(k_nn) || k_nn != as.integer(k_nn) ||
      k_nn < 1L || k_nn >= n)
    stop("'k_nn' must be an integer between 1 and n-1.")
  if (length(power) != 1L || !is.finite(power) || power < 0)
    stop("'power' must be a finite nonnegative number.")
  k_nn <- as.integer(k_nn)
  ii <- integer(n * k_nn); jj <- integer(n * k_nn); xx <- numeric(n * k_nn)
  for (start in seq(1L, n, by = chunk)) {
    rows <- start:min(n, start + chunk - 1L)
    d <- sqrt(outer(coords[rows, 1L], coords[, 1L], "-")^2 +
              outer(coords[rows, 2L], coords[, 2L], "-")^2)
    d[cbind(seq_along(rows), rows)] <- Inf
    for (r in seq_along(rows)) {
      i <- rows[r]
      nb <- order(d[r, ], seq_len(n))[seq_len(k_nn)]
      dn <- d[r, nb]
      # Equal weights at coincident locations avoid dividing by zero.
      wi <- if (any(dn == 0)) as.numeric(dn == 0) else
        if (power == 0) rep(1, k_nn) else dn^(-power)
      idx <- (i - 1L) * k_nn + seq_len(k_nn)
      ii[idx] <- i; jj[idx] <- nb; xx[idx] <- wi / sum(wi)
    }
  }
  if (sparse) {
    W <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
    return(Matrix::drop0(W))
  }
  W <- matrix(0, n, n)
  W[cbind(ii, jj)] <- xx
  W
}

#' k-nearest-neighbour weight matrix
#'
#' Row-standardized weights over the \code{k} nearest neighbours, equal
#' (\code{power = 0}) or inverse-distance (\code{power > 0}), computed in row
#' blocks so that the full distance matrix is never stored.
#'
#' @param coords n x 2 matrix of coordinates.
#' @param k Number of neighbours.
#' @param power Inverse-distance power; 0 gives equal weights.
#' @param sparse Logical; return a sparse \code{Matrix}.
#' @return Row-standardized n x n weight matrix.
#' @export
knn_weights <- function(coords, k = 4L, power = 0, sparse = TRUE) {
  coords <- as.matrix(coords)
  if (!is.numeric(coords) || ncol(coords) != 2L || any(!is.finite(coords)))
    stop("'coords' must be a finite n x 2 numeric matrix.", call. = FALSE)
  .knn_W_dgp(coords, k_nn = k, power = power, sparse = sparse)
}

#' Simulate data from the functional spatial-index logit model
#'
#' Generates the designs of the Monte Carlo experiments of the paper: the main
#' design (Experiments 1-3), the exact design (Experiment 2) and spatially
#' dependent responses with the same conditional mean (Experiment 2, Part C).
#'
#' @param n Number of units.
#' @param p Number of grid points on [0, 1] (default 101).
#' @param rho Spatial parameter, \eqn{|\rho| < 1}.
#' @param alpha Intercept (default 0).
#' @param nphi Number of Fourier pairs of the main design (default 10).
#' @param beta_scale Amplitude of \eqn{\beta(t) = }\code{beta_scale}
#'   \eqn{\sin(2\pi t)} in the main design (default 2, as in the paper).
#' @param k_nn Number of nearest neighbours of W (default 4).
#' @param W_type \code{"knn"} (default) or \code{"index"}.
#' @param power Inverse-distance power of the k-NN weights (default 1).
#' @param coords Optional n x 2 coordinates.
#' @param W Optional weight matrix (overrides \code{W_type}).
#' @param center_X Logical; center the Fourier coefficients (main design).
#' @param W_format \code{"dense"} or \code{"sparse"}.
#' @param design \code{"main"} (default) or \code{"exact"}.
#' @param exact_lambda,exact_sign,exact_sd_m Eigenvalues, signs of the
#'   coefficient function and standard deviation of the index in the exact
#'   design.
#' @param dependence_range Range of the exponential correlation of the
#'   Gaussian copula; 0 (default) gives conditionally independent responses.
#' @param response_model Only \code{"standardized_spatial_index_logit"}.
#' @return A list of class \code{"flsl_dgp"} with the responses \code{y}, the
#'   curves \code{x}, the weight matrix \code{w}, the grid, the true
#'   coefficient function, the coordinates, the probabilities, the true effect
#'   multipliers and diagnostics.
#' @export
dgp <- function(n = 200L, p = 101L, rho = 0.3, alpha = 0,
                nphi = 10L, beta_scale = 2, k_nn = 4L,
                W_type = c("knn", "index"), power = 1,
                coords = NULL, W = NULL, center_X = FALSE,
                W_format = c("dense", "sparse"),
                design = c("main", "exact"),
                exact_lambda = c(1, 0.7, 0.5, 0.35), exact_sign = c(1, -1, 1, -1),
                exact_sd_m = 1.5, dependence_range = 0,
                response_model = "standardized_spatial_index_logit") {
  W_type <- match.arg(W_type)
  W_format <- match.arg(W_format)
  design <- match.arg(design)
  if (!identical(response_model, "standardized_spatial_index_logit"))
    stop("Only response_model='standardized_spatial_index_logit' is implemented.")
  for (nm in c("n", "p", "nphi")) {
    v <- get(nm)
    if (length(v) != 1L || !is.finite(v) || v != as.integer(v) || v < 1L)
      stop(sprintf("'%s' must be a positive integer.", nm))
  }
  n <- as.integer(n); p <- as.integer(p); nphi <- as.integer(nphi)
  if (n < 3L || p < 5L) stop("Require n>=3 and p>=5.")
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1)
    stop("For a nonnegative row-standardized W, require abs(rho)<1.")
  if (length(alpha) != 1L || !is.finite(alpha) ||
      length(beta_scale) != 1L || !is.finite(beta_scale))
    stop("'alpha' and 'beta_scale' must be finite scalars.")
  if (length(dependence_range) != 1L || !is.finite(dependence_range) || dependence_range < 0)
    stop("'dependence_range' must be a nonnegative number.")
  if (design == "exact" && (length(exact_lambda) != length(exact_sign) ||
                            any(exact_lambda <= 0) || exact_sd_m <= 0))
    stop("Invalid exact-design settings.")

  ## Functional predictors (drawn first: same RNG order as the experiments)
  gp <- seq(0, 1, length.out = p)
  if (design == "main") {
    k <- seq_len(nphi)
    phi_cos <- outer(k, gp, function(kk, t) sqrt(2) * kk^(-1.5) * cos(kk * pi * t))
    phi_sin <- outer(k, gp, function(kk, t) sqrt(2) * kk^(-1.5) * sin(kk * pi * t))
    coeff <- matrix(stats::rnorm(n * 2L * nphi), nrow = n, ncol = 2L * nphi)
    if (center_X) coeff <- sweep(coeff, 2L, colMeans(coeff))
    x <- coeff %*% rbind(phi_cos, phi_sin)
    beta <- beta_scale * sin(2 * pi * gp)
    functional_effect <- as.vector(x %*% (beta * .trap_weights_dgp(gp)))
    beta_label <- "beta_scale * sin(2*pi*t)"
  } else {
    Kx <- length(exact_lambda)
    Phi <- sapply(seq_len(Kx), function(kk) sqrt(2) * sin(kk * pi * gp))     # p x K
    b_true <- exact_sign * exact_sd_m / sqrt(sum(exact_lambda * exact_sign^2))
    xi <- matrix(stats::rnorm(n * Kx), n, Kx) %*% diag(sqrt(exact_lambda), Kx)
    x <- xi %*% t(Phi)
    beta <- as.numeric(Phi %*% b_true)
    functional_effect <- as.numeric(xi %*% b_true)                        # exact index
    beta_label <- "c * sum_k sign_k * sqrt(2) sin(k pi t)"
  }

  ## Spatial weights
  W_supplied <- !is.null(W)
  if (!W_supplied) {
    if (W_type == "knn") {
      if (is.null(coords)) coords <- cbind(stats::runif(n), stats::runif(n))
      coords <- as.matrix(coords)
      if (!is.numeric(coords) || !identical(dim(coords), c(n, 2L)) || any(!is.finite(coords)))
        stop("'coords' must be a finite n x 2 numeric matrix.")
      W <- .knn_W_dgp(coords, k_nn = k_nn, power = power,
                      sparse = (W_format == "sparse"))
    } else {
      idx <- seq_len(n)
      W <- 1 / (1 + abs(outer(idx, idx, `-`)))
      diag(W) <- 0
      W <- W / rowSums(W)
      if (W_format == "sparse") W <- Matrix::Matrix(W, sparse = TRUE)
    }
  } else {
    if (!inherits(W, "sparseMatrix")) W <- as.matrix(W)
    if (nrow(W) != n || ncol(W) != n) stop("'W' must be n x n.")
    W <- normalize_weights(W)
    if (W_format == "dense") W <- as.matrix(W)
    if (W_format == "sparse" && !inherits(W, "sparseMatrix"))
      W <- Matrix::Matrix(W, sparse = TRUE)
    if (!is.null(coords)) coords <- as.matrix(coords)
  }

  ## Probabilities and true effect multipliers
  m <- alpha + functional_effect
  rcond_spatial <- NA_real_
  if (rho == 0) {
    sigma <- rep(1, n)
    eta <- m
    rcond_spatial <- 1
    probability <- stats::plogis(eta)
    c_direct <- mean(probability * (1 - probability)); c_total <- c_direct
  } else {
    Wd <- as.matrix(W)
    A <- -rho * Wd
    diag(A) <- diag(A) + 1
    if (n <= 2000L) rcond_spatial <- rcond(A)
    S <- solve(A)
    sigma <- sqrt(rowSums(S * S))
    eta <- as.vector(S %*% m) / sigma
    probability <- stats::plogis(eta)
    weight <- probability * (1 - probability) / sigma
    c_direct <- mean(weight * diag(S)); c_total <- mean(weight * rowSums(S))
    rm(S, A, Wd)
  }
  y <- stats::rbinom(n, size = 1L, prob = probability)

  ## Spatially dependent responses with the same conditional mean
  if (dependence_range > 0) {
    if (is.null(coords)) stop("'coords' are required for dependent responses.")
    D <- as.matrix(stats::dist(coords))
    e <- as.numeric(crossprod(chol(exp(-D / dependence_range)), stats::rnorm(n)))
    y <- as.integer(stats::pnorm(e) <= probability)
  }

  rs <- if (inherits(W, "sparseMatrix")) Matrix::rowSums(W) else rowSums(W)
  out <- list(
    y = as.integer(y), x = unname(x), w = W, gp = gp, beta = as.vector(beta),
    coords = coords, probability = probability, eta = eta, sigma = sigma,
    functional_effect = functional_effect, m = m,
    true_multipliers = c(direct = c_direct, indirect = c_total - c_direct, total = c_total),
    diagnostics = list(class_balance = mean(y), mean_probability = mean(probability),
                       probability_range = range(probability),
                       functional_effect_sd = stats::sd(functional_effect),
                       index_sd = stats::sd(eta),
                       row_sum_error = max(abs(rs - 1)),
                       rcond_spatial = rcond_spatial),
    settings = list(n = n, p = p, rho = rho, alpha = alpha, nphi = nphi,
                    beta_scale = beta_scale, k_nn = k_nn, design = design,
                    W_type = if (W_supplied) "supplied" else W_type,
                    W_format = W_format, power = power, center_X = center_X,
                    dependence_range = dependence_range,
                    response_model = response_model,
                    beta_function = beta_label)
  )
  class(out) <- c("flsl_dgp", "list")
  out
}
