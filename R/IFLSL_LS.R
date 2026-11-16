# ==============================================================
# IFLSL-LS: bounded, Armijo-safeguarded Gauss-Newton iteration for the
# projected GMM criterion (Section 3.3 and Algorithm 1, Appendix D).
#
#   * start at the nonspatial logit fit, theta = (delta0, 0); the first unit
#     step reproduces FLSL;
#   * direction (H + lambda I)^{-1} r, falling back to the negative gradient;
#   * projected Armijo rule, |rho| <= rho_max - eps_b;
#   * stop when the projected-gradient mapping is below grad_tol.
# Line-search trials evaluate probabilities only; the Jacobian is evaluated
# once per accepted iterate and uses rowSums(T * S), so S S' is never formed.
# The moments use the training units of the FLSL fit.
# ==============================================================

#' Iterated functional linearized spatial logit (IFLSL-LS)
#'
#' Computes the projected nonlinear GMM estimator by the bounded,
#' Armijo-safeguarded Gauss-Newton iteration of Algorithm 1 (Appendix D of the
#' paper), starting from the nonspatial logit fit stored in \code{flsl_fit}.
#'
#' @param y Binary response vector of length n; defaults to the response
#'   stored in \code{flsl_fit}.
#' @param flsl_fit A fitted \code{"FLSL"} object, which supplies the FPCA
#'   scores, W, the instruments and the training units.
#' @param rho_bound Bound rho_max of the spatial parameter (default 0.99).
#' @param maxit Maximum number of iterations N_max (default 1500).
#' @param max_halving Maximum number of step halvings I_max (default 40).
#' @param armijo Armijo constant c (default 1e-4).
#' @param grad_tol Projected-gradient tolerance eps_g (default 1e-5).
#' @param boundary_tol Boundary margin eps_b (default 1e-7).
#' @param ridge Regularization lambda of the Gauss-Newton matrix
#'   (default 1e-10).
#' @return An object of class \code{"IFLSL_LS"}.
#' @export
IFLSL_LS <- function(y = flsl_fit$y, flsl_fit, rho_bound = 0.99,
                     maxit = 1500L, max_halving = 40L,
                     armijo = 1e-4, grad_tol = 1e-5,
                     boundary_tol = 1e-7, ridge = 1e-10) {
  if (!inherits(flsl_fit, "FLSL")) stop("'flsl_fit' must inherit from 'FLSL'.")
  Xi <- as.matrix(flsl_fit$Xi)
  W <- as.matrix(flsl_fit$W)
  Ws <- .as_sparse_W(flsl_fit$W)
  Z <- as.matrix(flsl_fit$Z)
  n <- nrow(Xi); K <- ncol(Xi); d <- K + 2L
  tr <- if (is.null(flsl_fit$train)) seq_len(n) else flsl_fit$train
  y <- as.numeric(y)
  if (length(y) != n) stop("'y' must have length nrow(Xi).")
  ytr <- y[tr]
  if (any(!is.finite(ytr)) || !all(ytr %in% c(0, 1)))
    stop("'y' must be binary at the training units.")
  if (!all(dim(W) == c(n, n)) || nrow(Z) != n ||
      any(!is.finite(Xi)) || any(!is.finite(W)) || any(!is.finite(Z)))
    stop("Invalid Xi, W or Z matrices.")
  check_scalar <- function(v, lo, hi, name, closed_lo = FALSE) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) ||
        (if (closed_lo) v < lo else v <= lo) || v >= hi)
      stop(sprintf("Invalid %s.", name))
  }
  check_scalar(rho_bound, 0, 1, "rho_bound")
  check_scalar(armijo, 0, 1, "armijo")
  check_scalar(grad_tol, 0, Inf, "grad_tol")
  check_scalar(boundary_tol, 0, rho_bound, "boundary_tol")
  check_scalar(ridge, 0, Inf, "ridge", closed_lo = TRUE)
  valid_int <- function(v, lower) is.numeric(v) && length(v) == 1L &&
    is.finite(v) && v == floor(v) && v >= lower
  if (!valid_int(maxit, 1) || !valid_int(max_halving, 0))
    stop("Invalid maxit or max_halving.")
  maxit <- as.integer(maxit); max_halving <- as.integer(max_halving)

  Ztr <- Z[tr, , drop = FALSE]; ntr <- length(tr)
  qz <- qr(Ztr)
  if (qz$rank < ncol(Ztr) || ncol(Ztr) < d)
    stop("Instrument matrix rank deficient or too small.")
  bound <- rho_bound - min(boundary_tol, rho_bound / 100)
  project <- function(x) { x[d] <- max(-bound, min(bound, x[d])); x }
  pg_norm <- function(th, g) sqrt(sum((as.numeric(project(th - as.numeric(g)) - th))^2))

  theta <- stats::setNames(c(flsl_fit$alpha0, as.numeric(flsl_fit$b0), 0),
                           c("alpha", paste0("b", seq_len(K)), "rho"))
  if (length(theta) != d || any(!is.finite(theta)))
    stop("Invalid FLSL initial coefficients.")

  ## Criterion only (used for line-search trials)
  evaluate_J <- function(th) {
    e <- .flsl_probability(th[1], th[1 + seq_len(K)], th[d], Xi, W)
    Pu <- qr.fitted(qz, ytr - e$P[tr])
    sum(Pu^2) / ntr
  }
  ## Criterion, gradient and projected Jacobian (accepted iterates only)
  evaluate <- function(th) {
    e <- .flsl_general_gradient(th[1], th[1 + seq_len(K)], th[d], Xi, W, Ws)
    u <- as.numeric(ytr - e$P[tr])
    Pu <- as.numeric(qr.fitted(qz, u))
    PG <- qr.fitted(qz, e$G[tr, , drop = FALSE])
    list(J = sum(Pu^2) / ntr, gradient = as.numeric(-2 * crossprod(PG, Pu) / ntr),
         PG = PG, Pu = Pu, probability = as.numeric(e$P), residuals = u)
  }

  t_start <- .elapsed()
  cur <- evaluate(theta); history <- list(); status <- "maxit"
  n_trials <- 0L
  for (it in seq_len(maxit)) {
    g <- cur$gradient
    if (pg_norm(theta, g) <= grad_tol) { status <- "converged"; break }
    H <- crossprod(cur$PG)
    rhs <- as.numeric(crossprod(cur$PG, cur$Pu))
    delta <- tryCatch(as.numeric(solve(H + diag(ridge, d), rhs)),
                      error = function(e) rep(NA_real_, d))
    directions <- list(gauss_newton = delta, projected_gradient = -as.numeric(g))
    accepted <- FALSE; chosen <- NA_character_
    for (nm in names(directions)) {
      direction <- directions[[nm]]
      if (any(!is.finite(direction))) next
      step <- 1
      for (j in 0:max_halving) {
        candidate <- project(theta + step * direction)
        displacement <- as.numeric(candidate - theta)
        slope_actual <- sum(g * displacement)
        if (sum(displacement^2) > 0 && is.finite(slope_actual) && slope_actual < 0) {
          n_trials <- n_trials + 1L
          trialJ <- tryCatch(evaluate_J(candidate), error = function(e) NA_real_)
          if (is.finite(trialJ) &&
              trialJ <= cur$J + armijo * slope_actual + 1e-14) {
            accepted <- TRUE; chosen <- nm; break
          }
        }
        step <- step / 2
      }
      if (accepted) break
    }
    if (!accepted) { status <- "line_search_failed"; break }
    nxt <- tryCatch(evaluate(candidate), error = function(e) NULL)
    if (is.null(nxt)) { status <- "evaluation_failed"; break }
    update <- max(abs(displacement))
    theta <- stats::setNames(candidate, names(theta)); cur <- nxt
    history[[length(history) + 1L]] <- data.frame(
      iteration = it, alpha = unname(theta[1]), rho = unname(theta[d]),
      J = cur$J, step_size = step, max_update = update,
      gradient_norm = sqrt(sum(cur$gradient^2)),
      projected_gradient_norm = pg_norm(theta, cur$gradient),
      direction = chosen, slope_actual = slope_actual)
    if (pg_norm(theta, cur$gradient) <= grad_tol) { status <- "converged"; break }
  }
  if (status == "maxit" && pg_norm(theta, cur$gradient) <= grad_tol)
    status <- "converged"
  elapsed <- .elapsed() - t_start

  b_hat <- stats::setNames(as.numeric(theta[1 + seq_len(K)]), paste0("b", seq_len(K)))
  beta_hat <- NULL
  if (!is.null(flsl_fit$reduction) && !is.null(flsl_fit$gp)) {
    Phi <- fda::eval.fd(flsl_fit$gp, flsl_fit$reduction$PCAcoef)
    beta_hat <- as.numeric(Phi %*% b_hat)
  }
  result <- list(call = match.call(), method = "IFLSL-LS",
    instrument_order = flsl_fit$instrument_order,
    theta_hat = theta, alpha_hat = unname(theta[1]), b_hat = b_hat,
    rho_hat = unname(theta[d]), beta_hat = beta_hat,
    objective = cur$J, J = cur$J, gradient_norm = sqrt(sum(cur$gradient^2)),
    projected_gradient_norm = pg_norm(theta, cur$gradient),
    boundary = abs(theta[d]) >= bound - boundary_tol,
    converged = identical(status, "converged"),
    failed = !identical(status, "converged"),
    status = status, iterations = length(history), line_search_trials = n_trials,
    history = if (length(history)) do.call(rbind, history) else data.frame(),
    probability = cur$probability, fitted.values = cur$probability,
    residuals = cur$residuals, y = y, train = tr, n = n, n_train = ntr,
    Xi = Xi, W = flsl_fit$W, Z = Z, gp = flsl_fit$gp,
    reduction = flsl_fit$reduction, ncomp = K,
    rho_bound = rho_bound, effective_rho_bound = bound,
    grad_tol = grad_tol, maxit = maxit, max_halving = max_halving,
    armijo = armijo, boundary_tol = boundary_tol, ridge = ridge,
    one_step_theta = flsl_fit$theta_hat, elapsed_seconds = elapsed,
    flsl_fit = flsl_fit)
  class(result) <- "IFLSL_LS"
  result
}

#' @rdname IFLSL_LS
#' @param ... Arguments passed to \code{IFLSL_LS()}.
#' @export
IFLSL_LS_from_fit <- function(y = flsl_fit$y, flsl_fit, ...) IFLSL_LS(y = y, flsl_fit = flsl_fit, ...)

#' Projected-gradient norm (diagnostic helper)
#' @param theta Parameter vector (alpha, b, rho).
#' @param gradient Gradient of the criterion at theta.
#' @param rho_bound,boundary_tol Bound settings as in \code{IFLSL_LS()}.
#' @return The norm of the projected-gradient mapping.
#' @keywords internal
projected_gradient_norm <- function(theta, gradient, rho_bound = 0.99,
                                    boundary_tol = 1e-7) {
  if (length(theta) != length(gradient)) stop("theta and gradient lengths differ.")
  bound <- rho_bound - min(boundary_tol, rho_bound / 100)
  tnew <- as.numeric(theta) - as.numeric(gradient)
  tnew[length(tnew)] <- max(-bound, min(bound, tnew[length(tnew)]))
  sqrt(sum((tnew - as.numeric(theta))^2))
}
