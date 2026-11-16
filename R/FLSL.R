# ==============================================================
# Functional Linearized Spatial Logit (FLSL)
#
# FLSL()             : FPCA + FLSL_from_scores()
# FLSL_from_scores() : nonspatial logit start, one Gauss-Newton step for the
#                      projected GMM criterion (eq. 3.6-3.7 of the paper), and
#                      the inversion-free Wald test of H0: rho = 0
#                      (Corollary 3.1).
#
# No n x n matrix is inverted. W may be dense or a sparse Matrix; only
# products W %*% (n x K) are formed.
#
# 'train' restricts the moments (and the nonspatial start) to a subset of
# units while the curves and W of all units are kept. This is the fixed-
# location training/test design of the application (Section 5): withheld
# responses may be NA.
# Stage timings (seconds): fpca, logit, flsl, wald.
# ==============================================================


#' Resolve a training-unit specification
#' @keywords internal
.resolve_train <- function(train, n) {
  if (is.null(train)) return(seq_len(n))
  if (is.logical(train)) {
    if (length(train) != n || anyNA(train))
      stop("A logical 'train' must have length n and no missing values.", call. = FALSE)
    tr <- which(train)
  } else {
    if (!is.numeric(train) || anyNA(train) || any(train != floor(train)) ||
        any(train < 1) || any(train > n))
      stop("'train' must be a logical vector or integer indices in 1..n.", call. = FALSE)
    tr <- sort(unique(as.integer(train)))
  }
  if (length(tr) < 3L) stop("At least three training units are required.", call. = FALSE)
  tr
}


#' Functional linearized spatial logit model
#'
#' Computes the one-step functional linearized spatial logit (FLSL) estimator
#' and the inversion-free Wald test of no spatial propagation. The curves are
#' represented by FPCA on a cubic B-spline basis, the number of components is
#' chosen by the fraction-of-variance rule, and the instruments are
#' \eqn{Z = [1_n, \Xi, W\Xi, \ldots, W^L\Xi]}.
#'
#' @param y Binary response vector of length n (0/1). Responses of units not
#'   in \code{train} may be \code{NA}.
#' @param x Matrix of functional observations (rows = units, columns = grid).
#' @param W Spatial weight matrix (dense or sparse Matrix); it is
#'   row-standardized internally.
#' @param gp Optional observation grid; defaults to an equally spaced grid on
#'   [0, 1].
#' @param reduction Functional reduction method; only "FPCA".
#' @param nbasis Number of B-spline basis functions (default 15).
#' @param ncomp Optional fixed number of FPCA components.
#' @param threshold Fraction-of-variance threshold when \code{ncomp = NULL}
#'   (default 0.95).
#' @param eig_tol Numerical eigenvalue tolerance.
#' @param instrument_order Spatial lag order L of the instruments (default 1).
#' @param compute_wald Logical; compute the inversion-free Wald test.
#' @param train Optional training units (logical vector or indices). The
#'   moments use only these units; the FPCA, the curves and W use all units.
#' @return An object of class \code{"FLSL"}.
#' @export
FLSL <- function(y, x, W, gp = NULL, reduction = "FPCA", nbasis = 15,
                 ncomp = NULL, threshold = 0.95, eig_tol = 1e-10,
                 instrument_order = 1L, compute_wald = TRUE, train = NULL) {
  if (!identical(reduction, "FPCA")) stop("Only FPCA is supported.", call. = FALSE)
  if (!is.matrix(x)) x <- as.matrix(x)
  if (!is.numeric(x) || any(!is.finite(x)))
    stop("'x' must be a finite numeric matrix.", call. = FALSE)
  if (is.null(gp)) gp <- seq(0, 1, length.out = ncol(x))
  gp <- as.numeric(gp)
  if (length(gp) != ncol(x) || any(!is.finite(gp)))
    stop("'gp' must have length ncol(x).", call. = FALSE)

  t0 <- .elapsed()
  red <- getFPCA(data = x, nbasis = nbasis, gp = gp, ncomp = ncomp,
                 threshold = threshold, eig_tol = eig_tol)
  t_fpca <- .elapsed() - t0

  out <- FLSL_from_scores(y = y, Xi = red$PCAscore, W = W, reduction = red,
                          gp = gp, instrument_order = instrument_order,
                          compute_wald = compute_wald, train = train)
  out$call <- match.call()
  out$timing <- c(fpca = t_fpca, out$timing)
  out$threshold <- threshold
  out$eig_tol <- eig_tol
  out$nbasis <- nbasis
  out
}


#' FLSL from precomputed functional scores
#'
#' @inheritParams FLSL
#' @param Xi n x K matrix of functional scores.
#' @param reduction Optional \code{getFPCA()} object (needed for the
#'   coefficient function, the effect curves and prediction from curves).
#' @return An object of class \code{"FLSL"}.
#' @export
FLSL_from_scores <- function(y, Xi, W, reduction = NULL, gp = NULL,
                             instrument_order = 1L, compute_wald = TRUE,
                             train = NULL) {

  ## ------------------------------------------------------------
  ## 1. Input checks
  ## ------------------------------------------------------------
  if (length(instrument_order) != 1L || !is.numeric(instrument_order) ||
      !is.finite(instrument_order) || instrument_order < 1 ||
      instrument_order != floor(instrument_order))
    stop("'instrument_order' must be a single positive integer.", call. = FALSE)
  instrument_order <- as.integer(instrument_order)

  Xi <- as.matrix(Xi)
  if (!is.numeric(Xi) || any(!is.finite(Xi)))
    stop("'Xi' must be a finite numeric matrix.", call. = FALSE)
  n <- nrow(Xi); K <- ncol(Xi)
  if (K < 1L) stop("At least one functional component is required.", call. = FALSE)

  y <- as.numeric(y)
  if (length(y) != n) stop("length(y) must equal nrow(Xi).", call. = FALSE)
  tr <- .resolve_train(train, n)
  ytr <- y[tr]
  if (any(!is.finite(ytr)) || !all(ytr %in% c(0, 1)))
    stop("'y' must contain only 0 and 1 at the training units.", call. = FALSE)
  if (length(unique(ytr)) < 2L)
    stop("'y' must contain both response classes at the training units.", call. = FALSE)
  ntr <- length(tr)

  if (inherits(W, "Matrix") && !inherits(W, "sparseMatrix")) W <- as.matrix(W)
  if (!inherits(W, "sparseMatrix")) {
    W <- as.matrix(W)
    if (!is.numeric(W) || any(!is.finite(W)))
      stop("'W' must be a finite numeric matrix.", call. = FALSE)
  }
  if (nrow(W) != n || ncol(W) != n)
    stop("'W' must be an n x n matrix.", call. = FALSE)
  W <- normalize_weights(W)
  Wmul <- function(X) as.matrix(W %*% X)

  ## ------------------------------------------------------------
  ## 2. Nonspatial logit start (the FLR fit) on the training units
  ## ------------------------------------------------------------
  t0 <- .elapsed()
  X0 <- cbind(intercept = 1, Xi)
  colnames(X0) <- c("intercept", paste0("b", seq_len(K)))
  fit0 <- stats::glm.fit(x = X0[tr, , drop = FALSE], y = ytr,
                         family = stats::binomial(link = "logit"))
  t_logit <- .elapsed() - t0
  if (!isTRUE(fit0$converged))
    warning("The nonspatial logit start did not converge.", call. = FALSE)

  coef0 <- as.numeric(fit0$coefficients)
  if (length(coef0) != K + 1L || any(!is.finite(coef0)))
    stop("The nonspatial logistic regression produced non-finite estimates.",
         call. = FALSE)
  alpha0 <- coef0[1L]
  b0 <- stats::setNames(coef0[-1L], paste0("b", seq_len(K)))

  ## ------------------------------------------------------------
  ## 3. One Gauss-Newton step at theta0 = (delta0, 0)
  ##    G0 = R0 [B, W B delta0]   (eq. 3.5; no n x n inverse)
  ## ------------------------------------------------------------
  t0 <- .elapsed()
  m0 <- as.numeric(alpha0 + Xi %*% b0)
  P0 <- stats::plogis(m0)
  q0 <- P0 * (1 - P0)
  Wm0 <- as.numeric(Wmul(m0))
  G <- cbind(q0, Xi * q0, q0 * Wm0)
  colnames(G) <- c("alpha", paste0("b", seq_len(K)), "rho")

  if (1 + (instrument_order + 1) * K > ntr)
    stop("The instrument matrix has more columns than training units.", call. = FALSE)
  Z <- cbind(intercept = rep(1, n), Xi)
  colnames(Z) <- c("intercept", paste0("Xi", seq_len(K)))
  lagged_Xi <- vector("list", instrument_order)
  current_lag <- Xi
  for (ell in seq_len(instrument_order)) {
    current_lag <- Wmul(current_lag)
    colnames(current_lag) <- paste0("W", ell, "Xi", seq_len(K))
    Z <- cbind(Z, current_lag)
    lagged_Xi[[ell]] <- current_lag
  }
  Ztr <- Z[tr, , drop = FALSE]
  Gtr <- G[tr, , drop = FALSE]
  u0 <- ytr - P0[tr]

  qr_Z <- qr(Ztr)
  if (qr_Z$rank < ncol(Ztr))
    stop("The instrument matrix is rank deficient: rank = ", qr_Z$rank,
         ", required = ", ncol(Ztr), ".", call. = FALSE)
  Ghat <- qr.fitted(qr_Z, Gtr)
  qr_Ghat <- qr(Ghat)
  if (qr_Ghat$rank < ncol(Ghat))
    stop("The projected Jacobian is rank deficient: rank = ", qr_Ghat$rank,
         ", required = ", ncol(Ghat), ".", call. = FALSE)
  theta0 <- c(alpha = alpha0, b0, rho = 0)
  v <- as.numeric(u0 + Gtr %*% theta0)
  theta_hat <- as.numeric(qr.coef(qr_Ghat, v))
  names(theta_hat) <- colnames(G)
  if (any(!is.finite(theta_hat)))
    stop("The linearized IV step produced non-finite estimates.", call. = FALSE)
  t_flsl <- .elapsed() - t0

  ## ------------------------------------------------------------
  ## 4. Inversion-free Wald test of H0: rho = 0 (conditional independence)
  ##    V^(0) = sandwich (3.9) evaluated at theta0. For spatially dependent
  ##    responses use wald_rho(fit, type = "hac", ...).
  ## ------------------------------------------------------------
  wald <- NULL
  t_wald <- NA_real_
  if (isTRUE(compute_wald)) {
    t0 <- .elapsed()
    wald <- .flsl_wald(Z = Ztr, G0 = Gtr, q0 = q0[tr], rho_hat = theta_hat[["rho"]])
    t_wald <- .elapsed() - t0
  }

  ## ------------------------------------------------------------
  ## 5. Diagnostics (not timed) and beta_hat
  ## ------------------------------------------------------------
  kappa_Ghat <- tryCatch(kappa(Ghat, exact = TRUE), error = function(e) Inf)
  beta_hat <- NULL
  if (!is.null(reduction) && !is.null(gp)) {
    Phi <- fda::eval.fd(gp, reduction$PCAcoef)
    beta_hat <- as.numeric(Phi %*% theta_hat[paste0("b", seq_len(K))])
  }

  out <- list(
    method = "FLSL-FPCA", reduction_method = "FPCA",
    alpha_hat = unname(theta_hat["alpha"]),
    b_hat = theta_hat[paste0("b", seq_len(K))],
    rho_hat = unname(theta_hat["rho"]),
    theta_hat = theta_hat, beta_hat = beta_hat,
    admissible = abs(theta_hat[["rho"]]) < 1,
    converged = isTRUE(fit0$converged),
    logit_converged = isTRUE(fit0$converged),
    logit_boundary = isTRUE(fit0$boundary),
    y = y, train = tr, n = n, n_train = ntr,
    gp = gp, W = W, reduction = reduction, ncomp = K, Xi = Xi,
    WXi = lagged_Xi[[1L]], instrument_order = instrument_order,
    lagged_Xi = lagged_Xi,
    initial_fit = fit0, alpha0 = alpha0, b0 = b0, P0 = P0, q0 = q0, m0 = m0,
    u0 = u0, v = v, G = G, Ghat = Ghat, Z = Z, theta0 = theta0,
    rank_Z = qr_Z$rank, rank_Ghat = qr_Ghat$rank, n_parameters = K + 2L,
    kappa_Ghat = kappa_Ghat,
    wald = wald,
    timing = c(logit = t_logit, flsl = t_flsl, wald = t_wald)
  )
  class(out) <- "FLSL"
  out
}


#' Inversion-free Wald test of H0: rho = 0 (Corollary 3.1)
#'
#' Sandwich covariance (3.9) evaluated at the nonspatial logit fit, with
#' Omega computed under conditional independence from the fitted
#' probabilities or supplied (for example the spatial HAC estimate). Only
#' (K+2)- and q-dimensional systems are solved.
#' @keywords internal
.flsl_wald <- function(Z, G0, q0, rho_hat, level = 0.05, Omega = NULL) {
  n <- nrow(Z); d <- ncol(G0)
  A <- solve(crossprod(Z) / n)
  J <- -crossprod(Z, G0) / n
  Om <- if (is.null(Omega)) crossprod(Z * q0, Z) / n else Omega
  M <- crossprod(J, A %*% J)
  Mi <- solve(M)
  V <- Mi %*% t(J) %*% A %*% Om %*% A %*% J %*% Mi
  V <- (V + t(V)) / 2
  dimnames(V) <- list(colnames(G0), colnames(G0))
  se <- sqrt(pmax(diag(V), 0) / n)
  names(se) <- colnames(G0)
  se_rho <- se[d]
  stat <- if (is.finite(se_rho) && se_rho > 0) (rho_hat / se_rho)^2 else NA_real_
  list(statistic = unname(stat),
       p_value = if (is.finite(stat)) stats::pchisq(stat, df = 1, lower.tail = FALSE) else NA_real_,
       reject = if (is.finite(stat)) stat > stats::qchisq(1 - level, df = 1) else NA,
       level = level, se = se, se_rho = unname(se_rho), V = V)
}
