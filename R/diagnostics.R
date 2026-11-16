# ==============================================================
# Diagnostics used in the application (Section 5)
#
#   identification_r2()    : R^2 of each column of W Xi on [1, Xi]; the
#                            spatial parameter is separated from the
#                            coefficient function only if W Xi carries
#                            information beyond Xi (Appendix B).
#   fitted_probability()   : probabilities of a fit at all units.
#   pearson_residuals()    : Pearson residuals at the training units.
#   moran_test()           : Moran's I with a permutation p-value.
#   residual_correlogram() : binned correlation of standardized residuals by
#                            distance, with an optional permutation band;
#                            used to choose the HAC bandwidth.
# ==============================================================

#' Identification diagnostic for the spatial parameter
#'
#' Regresses each column of the spatially lagged scores \eqn{W\Xi} on
#' \eqn{[1_n, \Xi]} and returns the \eqn{R^2} values. Values close to one
#' indicate that the lagged scores carry little information beyond the scores
#' themselves, so that the spatial parameter is weakly identified.
#'
#' @param object A fitted \code{"FLSL"} object, or any fit that stores one in
#'   \code{flsl_fit}.
#' @return Named vector of \eqn{R^2} values, one per FPCA component.
#' @export
identification_r2 <- function(object) {
  base <- if (inherits(object, "FLSL")) object else object$flsl_fit
  if (is.null(base)) stop("An FLSL fit is required.", call. = FALSE)
  Xi <- as.matrix(base$Xi); WXi <- as.matrix(base$WXi); X <- cbind(1, Xi)
  r2 <- vapply(seq_len(ncol(Xi)), function(k) {
    f <- stats::lm.fit(X, WXi[, k])
    1 - sum(f$residuals^2) / sum((WXi[, k] - mean(WXi[, k]))^2)
  }, numeric(1))
  stats::setNames(r2, paste0("PC", seq_along(r2)))
}

#' Fitted probabilities at all units
#'
#' Returns the probabilities (2.3) at the estimate for all units of the
#' design, including units whose responses were withheld through
#' \code{train}; these are the predictions of the fixed-location
#' training/test design.
#'
#' @param object A fitted \code{"FLR"}, \code{"FLSL"}, \code{"IFLSL_LS"},
#'   \code{"NLGMM"} or \code{"fLSL_benchmark"} object.
#' @return Numeric vector of length n.
#' @export
fitted_probability <- function(object) {
  if (inherits(object, "FLR")) return(object$probability)
  if (inherits(object, "FLSL")) {
    if (abs(object$rho_hat) >= 1)
      stop("The one-step estimate is inadmissible (|rho| >= 1).", call. = FALSE)
    return(.flsl_probability(object$alpha_hat, object$b_hat, object$rho_hat,
                             object$Xi, as.matrix(object$W))$P)
  }
  if (is.null(object$probability)) stop("The fit has no probabilities.", call. = FALSE)
  object$probability
}

#' Pearson residuals of a fit
#'
#' @inheritParams fitted_probability
#' @return Numeric vector of length n with Pearson residuals at the training
#'   units and \code{NA} elsewhere.
#' @export
pearson_residuals <- function(object) {
  p <- fitted_probability(object)
  y <- if (inherits(object, "FLR")) object$initial_fit$y else object$y
  tr <- if (is.null(object$train)) seq_along(p) else object$train
  r <- rep(NA_real_, length(p))
  r[tr] <- (y[tr] - p[tr]) / sqrt(p[tr] * (1 - p[tr]))
  r
}

#' Moran's I with a permutation test
#'
#' @param x Numeric vector (for example Pearson residuals).
#' @param W n x n weight matrix (dense or sparse).
#' @param nperm Number of random permutations (default 999).
#' @return A list with Moran's I, its expectation under no autocorrelation,
#'   the one-sided permutation p-value and the number of permutations.
#' @export
moran_test <- function(x, W, nperm = 999L) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) stop("'x' must be finite.", call. = FALSE)
  if (nrow(W) != length(x) || ncol(W) != length(x)) stop("'W' must be n x n.", call. = FALSE)
  s0 <- sum(W)
  moran <- function(z) {
    z <- z - mean(z)
    length(z) / s0 * as.numeric(crossprod(z, as.numeric(W %*% z))) / sum(z^2)
  }
  I0 <- moran(x)
  Ip <- vapply(seq_len(nperm), function(i) moran(sample(x)), numeric(1))
  list(I = I0, expected = -1 / (length(x) - 1),
       p_value = (1 + sum(Ip >= I0)) / (nperm + 1), nperm = nperm)
}

#' Residual correlogram
#'
#' Average product of standardized residuals over pairs of units in distance
#' classes, with an optional permutation band under no spatial correlation.
#' It can be used to choose the HAC bandwidth as the distance beyond which
#' the residual correlation stays within the band.
#'
#' @param x Numeric vector of residuals.
#' @param coords n x 2 matrix of coordinates.
#' @param breaks Distance class boundaries (in the units of \code{coords}).
#' @param nperm Number of permutations for the band (0 for none).
#' @param level Coverage of the permutation band.
#' @return Data frame with the distance classes, number of pairs, the
#'   correlation and, when \code{nperm > 0}, the band.
#' @details Uses all pairs of units and is therefore intended for samples of
#'   up to a few thousand units.
#' @export
residual_correlogram <- function(x, coords, breaks, nperm = 0L, level = 0.95) {
  x <- as.numeric(x); coords <- as.matrix(coords)
  if (any(!is.finite(x)) || nrow(coords) != length(x))
    stop("'x' must be finite with one value per row of 'coords'.", call. = FALSE)
  D <- as.matrix(stats::dist(coords)); up <- upper.tri(D)
  cls <- cut(D[up], breaks, include.lowest = TRUE)
  ii <- row(D)[up]; jj <- col(D)[up]; rm(D)
  correl <- function(r) {
    z <- (r - mean(r)) / sqrt(mean((r - mean(r))^2))
    as.numeric(tapply(z[ii] * z[jj], cls, mean))
  }
  out <- data.frame(lower = utils::head(breaks, -1L), upper = breaks[-1L],
                    n_pairs = as.numeric(table(cls)), correlation = correl(x))
  if (nperm > 0) {
    perm <- replicate(nperm, correl(sample(x)))
    a <- (1 - level) / 2
    out$band_lower <- apply(perm, 1, stats::quantile, a, na.rm = TRUE)
    out$band_upper <- apply(perm, 1, stats::quantile, 1 - a, na.rm = TRUE)
  }
  out
}
