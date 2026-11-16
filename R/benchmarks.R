# ==============================================================
# Benchmark estimators of Section 3.1 / 4.1: Bernoulli QMLE and NLS for the
# same probability model, the same FPCA scores and the same W, with the
# bound |rho| <= rho_bound - 1e-7 and the analytic gradient supplied to
# L-BFGS-B. Both evaluate S_rho at every evaluation of the criterion.
# ==============================================================

#' Bernoulli QMLE or NLS benchmark for the functional spatial-index logit
#'
#' @param y Binary response vector of length n; defaults to the response
#'   stored in \code{flsl_fit}.
#' @param flsl_fit A fitted \code{"FLSL"} object (scores, W, training units
#'   and the nonspatial starting value).
#' @param method \code{"QMLE"} (default) or \code{"NLS"}.
#' @param rho_bound Bound of the spatial parameter (default 0.99).
#' @param maxit Maximum number of L-BFGS-B iterations.
#' @return An object of class \code{"fLSL_benchmark"}.
#' @export
spatial_benchmark <- function(y = flsl_fit$y, flsl_fit, method = c("QMLE", "NLS"),
                              rho_bound = 0.99, maxit = 2000L) {
  method <- match.arg(method)
  if (!inherits(flsl_fit, "FLSL")) stop("'flsl_fit' must inherit from 'FLSL'.", call. = FALSE)
  Xi <- as.matrix(flsl_fit$Xi); Wd <- as.matrix(flsl_fit$W); Ws <- .as_sparse_W(flsl_fit$W)
  n <- nrow(Xi); K <- ncol(Xi); d <- K + 2L
  tr <- if (is.null(flsl_fit$train)) seq_len(n) else flsl_fit$train
  y <- as.numeric(y)
  if (length(y) != n) stop("'y' must have length n.", call. = FALSE)
  ytr <- y[tr]; ntr <- length(tr)
  if (any(!is.finite(ytr)) || !all(ytr %in% c(0, 1)))
    stop("'y' must be binary at the training units.", call. = FALSE)

  prob <- function(th) {
    if (any(!is.finite(th)) || abs(th[d]) >= 1) return(NULL)
    .flsl_probability(th[1], th[1 + seq_len(K)], th[d], Xi, Wd)$P
  }
  fn <- function(th) {
    pp <- tryCatch(prob(th), error = function(e) NULL)
    if (is.null(pp) || any(!is.finite(pp))) return(1e15)
    pp <- pp[tr]
    if (method == "QMLE") {
      pc <- pmin(1 - 1e-12, pmax(1e-12, pp))
      return(-mean(ytr * log(pc) + (1 - ytr) * log1p(-pc)))
    }
    mean((ytr - pp)^2)
  }
  gr <- function(th) {
    o <- tryCatch(.flsl_general_gradient(th[1], th[1 + seq_len(K)], th[d], Xi, Wd, Ws),
                  error = function(e) NULL)
    if (is.null(o) || any(!is.finite(o$G))) return(rep(0, d))
    P <- o$P[tr]; G <- o$G[tr, , drop = FALSE]
    r <- if (method == "QMLE") (ytr - P) / pmax(P * (1 - P), 1e-12) else 2 * (ytr - P)
    -as.numeric(crossprod(G, r)) / ntr
  }
  start <- c(flsl_fit$alpha0, as.numeric(flsl_fit$b0), 0)
  lower <- c(rep(-Inf, K + 1L), -rho_bound + 1e-7)
  upper <- c(rep(Inf, K + 1L), rho_bound - 1e-7)
  t0 <- .elapsed()
  opt <- stats::optim(start, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(maxit = as.integer(maxit), factr = 1e7))
  elapsed <- .elapsed() - t0
  th <- stats::setNames(opt$par, c("alpha", paste0("b", seq_len(K)), "rho"))
  P <- tryCatch(prob(th), error = function(e) rep(NA_real_, n))
  beta_hat <- if (!is.null(flsl_fit$reduction) && !is.null(flsl_fit$gp))
    as.numeric(fda::eval.fd(flsl_fit$gp, flsl_fit$reduction$PCAcoef) %*% th[1 + seq_len(K)]) else NULL
  out <- list(method = method, theta_hat = th, alpha_hat = unname(th[1]),
              b_hat = th[1 + seq_len(K)], rho_hat = unname(th[d]), beta_hat = beta_hat,
              converged = (opt$convergence == 0L), failed = (opt$convergence != 0L),
              status = paste0(method, "_optim_", opt$convergence),
              iterations = as.numeric(opt$counts[["function"]]),
              criterion_value = opt$value, probability = P, fitted.values = P,
              y = y, train = tr, n = n, n_train = ntr, ncomp = K,
              Xi = Xi, W = flsl_fit$W, gp = flsl_fit$gp, reduction = flsl_fit$reduction,
              elapsed_seconds = elapsed, flsl_fit = flsl_fit)
  class(out) <- "fLSL_benchmark"
  out
}
