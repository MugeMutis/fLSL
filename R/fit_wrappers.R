#' Fit IFLSL-LS directly from functional observations
#'
#' Computes the FLSL fit (FPCA, nonspatial start, one-step estimator and its
#' test) and then the IFLSL-LS iteration from it.
#'
#' @param y Binary response vector.
#' @param x Matrix of functional observations.
#' @param W Spatial weight matrix.
#' @param instrument_order Spatial lag order of the instruments (default 1).
#' @param rho_bound,maxit,max_halving,armijo,grad_tol,boundary_tol,ridge
#'   Settings of Algorithm 1; see \code{IFLSL_LS()}.
#' @param ... Additional arguments passed to \code{FLSL()} (for example
#'   \code{gp}, \code{nbasis}, \code{threshold} or \code{train}).
#' @return An object of class \code{"IFLSL_LS"}.
#' @export
fit_IFLSL_LS <- function(y, x, W, instrument_order = 1L,
                         rho_bound = .99, maxit = 1500L,
                         max_halving = 40L, armijo = 1e-4,
                         grad_tol = 1e-5, boundary_tol = 1e-7, ridge = 1e-10, ...) {
  base <- FLSL(y = y, x = x, W = W, instrument_order = instrument_order, ...)
  IFLSL_LS(y = y, flsl_fit = base, rho_bound = rho_bound,
           maxit = maxit, max_halving = max_halving,
           armijo = armijo, grad_tol = grad_tol,
           boundary_tol = boundary_tol, ridge = ridge)
}

#' Fit direct nonlinear GMM directly from functional observations
#'
#' @inheritParams fit_IFLSL_LS
#' @param rho_bound Bound of the spatial parameter (default 0.99).
#' @param rho_starts,maxit,kkt_tol,scale_step,refine Settings of
#'   \code{NLGMM()}.
#' @return An object of class \code{"NLGMM"}.
#' @export
fit_NLGMM <- function(y, x, W, instrument_order = 1L,
                      rho_bound = .99, rho_starts = c(-.9, -.5, 0, .5, .9),
                      maxit = 2000L, kkt_tol = 1e-5, scale_step = 1,
                      refine = TRUE, ...) {
  base <- FLSL(y = y, x = x, W = W, instrument_order = instrument_order, ...)
  NLGMM(y = y, flsl_fit = base, rho_bound = rho_bound,
        rho_starts = rho_starts, maxit = maxit,
        kkt_tol = kkt_tol, scale_step = scale_step, refine = refine)
}
