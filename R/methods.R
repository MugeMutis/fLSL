# Print and coefficient methods.

.print_fit_header <- function(x, title) {
  cat(title, "\n", sep = "")
  ntr <- if (is.null(x$n_train)) NA else x$n_train
  cat(sprintf("  units: %d (training: %d), FPCA components: K = %d\n",
              as.integer(x$n), as.integer(ntr), as.integer(x$ncomp)))
}

#' @export
print.FLSL <- function(x, digits = 4, ...) {
  .print_fit_header(x, "Functional linearized spatial logit (FLSL, one-step)")
  cat("  rho_hat =", format(x$rho_hat, digits = digits),
      if (!isTRUE(x$admissible)) "(inadmissible, |rho| >= 1)" else "", "\n")
  if (!is.null(x$wald))
    cat(sprintf("  Inversion-free Wald test of rho = 0 (conditional independence): %s, p = %s\n",
                format(x$wald$statistic, digits = digits), format(x$wald$p_value, digits = digits)))
  invisible(x)
}

#' @export
print.IFLSL_LS <- function(x, digits = 4, ...) {
  .print_fit_header(x, "Iterated functional linearized spatial logit (IFLSL-LS)")
  cat("  rho_hat =", format(x$rho_hat, digits = digits), "\n")
  cat(sprintf("  status: %s after %d iterations; criterion = %s\n",
              x$status, as.integer(x$iterations), format(x$objective, digits = digits)))
  invisible(x)
}

#' @export
print.NLGMM <- function(x, digits = 4, ...) {
  cat("Direct nonlinear GMM (NLGMM)\n")
  cat("  rho_hat =", format(x$rho_hat, digits = digits), "; status:", x$status, "\n")
  invisible(x)
}

#' @export
print.fLSL_benchmark <- function(x, digits = 4, ...) {
  cat("Benchmark estimator:", x$method, "\n")
  cat("  rho_hat =", format(x$rho_hat, digits = digits), "; status:", x$status, "\n")
  invisible(x)
}

#' @export
coef.FLSL <- function(object, ...) object$theta_hat

#' @export
coef.IFLSL_LS <- function(object, ...) object$theta_hat

#' @export
coef.NLGMM <- function(object, ...) object$theta_hat

#' @export
coef.fLSL_benchmark <- function(object, ...) object$theta_hat

#' @export
coef.FLR <- function(object, ...) object$coefficients
