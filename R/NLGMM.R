# ==============================================================
# Direct nonlinear GMM (benchmark of Section 4.1): minimization of the
# projected GMM criterion (3.3) by L-BFGS-B with the same probability, Z and
# Q as IFLSL-LS. rho_starts = 0 gives the single-start version (NLGMM-S);
# the default five starting values with a polishing run give NLGMM-M.
# The criterion and gradient are cached so that L-BFGS-B, which requests fn
# and gr at the same point, triggers a single Jacobian evaluation.
# ==============================================================

#' Direct nonlinear GMM for the functional spatial-index logit
#'
#' @param y Binary response vector of length n; defaults to the response
#'   stored in \code{flsl_fit}.
#' @param flsl_fit A fitted \code{"FLSL"} object.
#' @param rho_bound Bound of the spatial parameter (default 0.99).
#' @param rho_starts Starting values of rho; \code{0} gives NLGMM-S, the
#'   default \code{c(-0.9, -0.5, 0, 0.5, 0.9)} gives NLGMM-M.
#' @param maxit Maximum number of L-BFGS-B iterations per start.
#' @param kkt_tol Projected-gradient tolerance used to qualify a solution.
#' @param scale_step Step scaling of the projected-gradient mapping.
#' @param refine Logical; polish the best solutions with a second run.
#' @return An object of class \code{"NLGMM"}.
#' @export
NLGMM <- function(y = flsl_fit$y, flsl_fit, rho_bound = 0.99,
                  rho_starts = c(-0.9, -0.5, 0, 0.5, 0.9),
                  maxit = 2000L, kkt_tol = 1e-5, scale_step = 1,
                  refine = TRUE) {
  if (!inherits(flsl_fit, "FLSL")) stop("flsl_fit must be an FLSL object")
  Xi <- as.matrix(flsl_fit$Xi)
  W <- as.matrix(flsl_fit$W)
  Ws <- .as_sparse_W(flsl_fit$W)
  Z <- as.matrix(flsl_fit$Z); y <- as.numeric(y)
  n <- nrow(Xi); K <- ncol(Xi); d <- K + 2L
  tr <- if (is.null(flsl_fit$train)) seq_len(n) else flsl_fit$train
  if (length(y) != n) stop("Invalid binary y")
  ytr <- y[tr]
  if (any(!is.finite(ytr)) || !all(ytr %in% c(0, 1))) stop("Invalid binary y")
  if (length(rho_bound) != 1L || !is.finite(rho_bound) || rho_bound <= 0 || rho_bound >= 1)
    stop("Invalid rho_bound")
  if (length(kkt_tol) != 1L || !is.finite(kkt_tol) || kkt_tol <= 0 ||
      length(scale_step) != 1L || !is.finite(scale_step) || scale_step <= 0)
    stop("Invalid tolerances")
  if (length(maxit) != 1L || !is.finite(maxit) || maxit < 1 || maxit != floor(maxit))
    stop("Invalid maxit")
  if (length(refine) != 1L || is.na(refine) || !is.logical(refine))
    stop("refine must be TRUE or FALSE")
  if (!is.numeric(rho_starts) || !length(rho_starts) || any(!is.finite(rho_starts)))
    stop("Invalid rho_starts")
  if (nrow(W) != n || ncol(W) != n || nrow(Z) != n)
    stop("Incompatible spatial/instrument dimensions")
  Ztr <- Z[tr, , drop = FALSE]; ntr <- length(tr)
  qz <- qr(Ztr)
  if (qz$rank < ncol(Ztr)) stop("Instrument matrix rank deficient")
  lower <- c(rep(-Inf, d - 1L), -rho_bound + 1e-7)
  upper <- c(rep(Inf, d - 1L), rho_bound - 1e-7)
  initial <- c(flsl_fit$alpha0, as.numeric(flsl_fit$b0))
  if (length(initial) != d - 1L || any(!is.finite(initial))) stop("Invalid FLSL start")

  evaluate <- function(par) {
    e <- .flsl_general_gradient(par[1L], par[seq_len(K) + 1L], par[d], Xi, W, Ws)
    u <- ytr - as.numeric(e$P[tr])
    pu <- as.numeric(qr.fitted(qz, u))
    pg <- as.matrix(qr.fitted(qz, e$G[tr, , drop = FALSE]))
    list(value = sum(pu^2) / ntr,
         gradient = as.numeric(-2 * crossprod(pg, pu) / ntr),
         pg = pg, probability = as.numeric(e$P), residuals = u)
  }
  cache <- new.env(parent = emptyenv())
  evaluate_cached <- function(par) {
    if (!is.null(cache$par) && length(cache$par) == length(par) &&
        all(cache$par == par)) return(cache$val)
    val <- evaluate(par)
    cache$par <- par; cache$val <- val
    val
  }
  mapping <- function(par, grad) {
    projected <- pmax(lower, pmin(upper, par - scale_step * grad))
    (par - projected) / scale_step
  }
  starts <- lapply(unique(pmax(lower[d], pmin(upper[d], rho_starts))),
                   function(r) c(initial, r))
  fn <- function(p) evaluate_cached(p)$value
  gr <- function(p) evaluate_cached(p)$gradient
  run_optim <- function(s, phase) {
    tryCatch({
      opt <- stats::optim(s, fn, gr, method = "L-BFGS-B", lower = lower, upper = upper,
                          control = list(maxit = as.integer(maxit), factr = 1e3, pgtol = 1e-11))
      ev <- evaluate_cached(opt$par)
      pg <- sqrt(sum(mapping(opt$par, ev$gradient)^2))
      list(par = opt$par, ev = ev, pg = pg, rank = qr(ev$pg)$rank,
           code = opt$convergence, message = as.character(opt$message),
           counts = opt$counts[["function"]], phase = phase)
    }, error = function(e) list(par = NULL, ev = NULL, pg = Inf, rank = NA_integer_,
                                code = 999L, message = conditionMessage(e),
                                counts = NA_real_, phase = phase))
  }
  t0 <- proc.time()[["elapsed"]]
  fits <- lapply(seq_along(starts), function(i) run_optim(starts[[i]], paste0("start_", i)))
  if (refine) {
    finite <- which(vapply(fits, function(f) !is.null(f$ev) && is.finite(f$ev$value), logical(1)))
    if (length(finite)) {
      vals <- vapply(fits[finite], function(f) f$ev$value, numeric(1))
      good <- finite[vapply(fits[finite], function(f)
        f$code == 0L && f$pg <= kkt_tol && isTRUE(f$rank == d), logical(1))]
      indices <- finite[which.min(vals)]
      if (length(good))
        indices <- unique(c(indices, good[which.min(vapply(fits[good], function(f) f$ev$value, numeric(1)))]))
      for (j in indices) fits[[length(fits) + 1L]] <- run_optim(fits[[j]]$par, paste0("polish_", j))
    }
  }
  attempts <- do.call(rbind, lapply(seq_along(fits), function(i) {
    f <- fits[[i]]
    data.frame(attempt = i, phase = f$phase,
      objective = if (is.null(f$ev)) NA_real_ else f$ev$value,
      rho = if (is.null(f$par)) NA_real_ else f$par[d],
      projected_gradient_norm = f$pg, rank = f$rank,
      convergence_code = f$code, function_evaluations = f$counts,
      kkt_ok = is.finite(f$pg) && f$pg <= kkt_tol,
      message = f$message, stringsAsFactors = FALSE)
  }))
  valid <- which(is.finite(attempts$objective))
  if (!length(valid)) {
    out <- list(method = "NLGMM", converged = FALSE, failed = TRUE, status = "all_failed",
                attempts = attempts, elapsed_seconds = proc.time()[["elapsed"]] - t0,
                flsl_fit = flsl_fit)
    class(out) <- "NLGMM"
    return(out)
  }
  qualified <- valid[attempts$convergence_code[valid] == 0L &
                       attempts$kkt_ok[valid] & attempts$rank[valid] %in% d]
  pool <- if (length(qualified)) qualified else valid
  chosen <- pool[which.min(attempts$objective[pool])]
  raw_best <- valid[which.min(attempts$objective[valid])]
  f <- fits[[chosen]]
  theta <- stats::setNames(as.numeric(f$par), c("alpha", paste0("b", seq_len(K)), "rho"))
  b <- stats::setNames(unname(theta[seq_len(K) + 1L]), paste0("b", seq_len(K)))
  beta <- if (!is.null(flsl_fit$reduction) && !is.null(flsl_fit$gp))
    as.numeric(fda::eval.fd(flsl_fit$gp, flsl_fit$reduction$PCAcoef) %*% b) else NULL
  success <- chosen %in% qualified
  out <- list(method = "NLGMM", instrument_order = flsl_fit$instrument_order,
       theta_hat = theta, alpha_hat = unname(theta[1L]), b_hat = b,
       rho_hat = unname(theta[d]), beta_hat = beta,
       probability = f$ev$probability, fitted.values = f$ev$probability,
       residuals = f$ev$residuals, objective = f$ev$value,
       gradient_norm = sqrt(sum(f$ev$gradient^2)),
       projected_gradient_norm = f$pg, rank = f$rank, required_rank = d,
       converged = success, failed = !success,
       status = if (success) "converged" else "not_kkt_qualified",
       convergence_code = f$code, boundary = abs(theta[d]) >= rho_bound - 1e-5,
       iterations = sum(attempts$function_evaluations, na.rm = TRUE),
       selected_attempt = chosen, raw_best_attempt = raw_best,
       raw_best_objective = attempts$objective[raw_best],
       objective_gap_from_raw_best = f$ev$value - attempts$objective[raw_best],
       attempts = attempts, elapsed_seconds = proc.time()[["elapsed"]] - t0,
       y = y, train = tr, n = n, n_train = ntr,
       Xi = Xi, W = flsl_fit$W, Z = Z, gp = flsl_fit$gp, reduction = flsl_fit$reduction,
       ncomp = K, flsl_fit = flsl_fit)
  class(out) <- "NLGMM"
  out
}
