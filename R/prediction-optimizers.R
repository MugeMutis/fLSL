# Prediction methods for optimizer-based fits (IFLSL-LS, NLGMM, QMLE/NLS).
# The FLSL object stored in the fit supplies the training FPCA
# transformation; the probabilities are computed by the same engine as for
# FLSL. FPCA is never refitted on test curves.

.predict_optimized_flsl <- function(object, newdata = NULL, scores = NULL,
                                    Wnew = NULL, type = c("probability", "class", "link", "all"),
                                    cutoff = 0.5) {
  type <- match.arg(type)
  if (!isTRUE(object$converged) || isTRUE(object$failed)) {
    stop("The optimization did not converge successfully; prediction is unavailable.",
         call. = FALSE)
  }
  base <- object$flsl_fit
  if (is.null(base) || !inherits(base, "FLSL")) {
    stop("The training FLSL fit is unavailable.", call. = FALSE)
  }
  if (is.null(object$alpha_hat) || is.null(object$b_hat) ||
      is.null(object$rho_hat)) {
    stop("Optimized parameter estimates are missing.", call. = FALSE)
  }
  if (length(object$alpha_hat) != 1L || length(object$rho_hat) != 1L ||
      length(object$b_hat) != base$ncomp ||
      any(!is.finite(c(object$alpha_hat, object$b_hat, object$rho_hat)))) {
    stop("Optimized parameters are invalid or inconsistent with training FPCA.",
         call. = FALSE)
  }
  base$alpha_hat <- as.numeric(object$alpha_hat)
  base$b_hat <- as.numeric(object$b_hat)
  base$rho_hat <- as.numeric(object$rho_hat)
  predict.FLSL(base, newdata = newdata, scores = scores, Wnew = Wnew,
               type = type, cutoff = cutoff)
}

#' Predict from a converged IFLSL-LS fit
#' @param object A fitted object from \code{IFLSL_LS()} or \code{fit_IFLSL_LS()}.
#' @param newdata New functional observations on the training grid.
#' @param scores Alternatively, scores in the training FPCA basis.
#' @param Wnew Weight matrix of the prediction sample; not used when
#'   \code{newdata} and \code{scores} are \code{NULL}, in which case the
#'   fitted probabilities of all units are returned.
#' @param type One of probability, class, link, all.
#' @param cutoff Classification threshold.
#' @param ... Unused.
#' @export
predict.IFLSL_LS <- function(object, newdata = NULL, scores = NULL, Wnew = NULL,
                             type = c("probability", "class", "link", "all"),
                             cutoff = 0.5, ...) {
  .predict_optimized_flsl(object, newdata = newdata, scores = scores,
                          Wnew = Wnew, type = type, cutoff = cutoff)
}

#' Predict from a converged direct nonlinear GMM fit
#' @inheritParams predict.IFLSL_LS
#' @param object A fitted object from \code{NLGMM()} or \code{fit_NLGMM()}.
#' @export
predict.NLGMM <- function(object, newdata = NULL, scores = NULL, Wnew = NULL,
                          type = c("probability", "class", "link", "all"),
                          cutoff = 0.5, ...) {
  .predict_optimized_flsl(object, newdata = newdata, scores = scores,
                          Wnew = Wnew, type = type, cutoff = cutoff)
}

#' Predict from a converged QMLE or NLS benchmark fit
#' @inheritParams predict.IFLSL_LS
#' @param object A fitted object from \code{spatial_benchmark()}.
#' @export
predict.fLSL_benchmark <- function(object, newdata = NULL, scores = NULL, Wnew = NULL,
                                   type = c("probability", "class", "link", "all"),
                                   cutoff = 0.5, ...) {
  .predict_optimized_flsl(object, newdata = newdata, scores = scores,
                          Wnew = Wnew, type = type, cutoff = cutoff)
}
