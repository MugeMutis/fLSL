# ==============================================================
# Functional effects (Sections 2.2 and 3.6)
#
#   ADE(t) = c_D(theta) beta(t),  AIE(t) = c_I(theta) beta(t),
#   ATE(t) = {c_D(theta) + c_I(theta)} beta(t),
#   c_D = n^{-1} sum_i p_i(1-p_i)(S)_ii / s_i,
#   c_D + c_I = n^{-1} sum_i p_i(1-p_i) (S 1)_i / s_i,
# averaged over all units of the observed design, with pointwise
# delta-method bands (3.12) based on vcov_theta(). The gradient of the
# multipliers is obtained by central differences.
# ==============================================================

#' Effect multipliers of the functional spatial-index logit
#'
#' @param theta Parameter vector (alpha, b_1, ..., b_K, rho).
#' @param Xi n x K matrix of FPCA scores.
#' @param W Row-standardized n x n weight matrix.
#' @return Named vector \code{c(direct, indirect, total)}.
#' @export
effect_multipliers <- function(theta, Xi, W) {
  Xi <- as.matrix(Xi); W <- as.matrix(W)
  K <- ncol(Xi); rho <- theta[K + 2L]
  if (length(theta) != K + 2L || !is.finite(rho) || abs(rho) >= 1)
    return(c(direct = NA_real_, indirect = NA_real_, total = NA_real_))
  m <- as.numeric(theta[1] + Xi %*% theta[1 + seq_len(K)])
  n <- length(m)
  S <- if (rho == 0) diag(n) else .spatial_inverse(rho, W)
  sigma <- sqrt(rowSums(S * S))
  P <- stats::plogis(as.numeric(S %*% m) / sigma)
  weight <- P * (1 - P) / sigma
  cd <- mean(weight * diag(S))
  ct <- mean(weight * rowSums(S))
  c(direct = cd, indirect = ct - cd, total = ct)
}

#' Functional effect curves with pointwise bands
#'
#' Computes the estimated coefficient function and the average direct,
#' indirect and total functional effect curves (2.8)-(2.10) with pointwise
#' delta-method bands (3.12). The bands treat the observed design and the
#' FPCA basis as fixed and are pointwise, not simultaneous.
#'
#' @param object A fitted \code{"FLSL"}, \code{"IFLSL_LS"}, \code{"NLGMM"},
#'   \code{"fLSL_benchmark"} or \code{"FLR"} object fitted from curves (so
#'   that the FPCA basis is available).
#' @param vcov Optional result of \code{vcov_theta()}; computed from
#'   \code{type}, \code{coords} and \code{bandwidth} when \code{NULL}.
#' @inheritParams vcov_theta
#' @param level Nominal level of the bands (default 0.95).
#' @param h Step of the central differences.
#' @return An object of class \code{"fLSL_effects"}: a list with the data
#'   frame \code{curves} (columns curve, t, estimate, se, lower, upper), the
#'   multipliers and the indirect share \eqn{c_I/(c_D + c_I)}.
#' @details The multipliers require the dense spatial multiplier at the
#'   estimate and are therefore limited to sample sizes for which an n x n
#'   inverse is feasible.
#' @export
functional_effects <- function(object, vcov = NULL, type = c("independent", "hac"),
                               coords = NULL, bandwidth = NULL, level = 0.95, h = 1e-5) {
  type <- match.arg(type)
  base <- if (inherits(object, c("FLSL", "FLR"))) object else object$flsl_fit
  red <- object$reduction; gp <- object$gp
  if (is.null(red) || is.null(gp))
    stop("The FPCA basis is required; fit the model from curves.", call. = FALSE)
  if (is.null(vcov)) vcov <- vcov_theta(object, type = type, coords = coords, bandwidth = bandwidth)
  theta <- vcov$theta; V <- vcov$V; n <- vcov$n; d <- length(theta); K <- d - 2L
  if (abs(theta[d]) >= 1) stop("The estimate is inadmissible (|rho| >= 1).", call. = FALSE)
  Xi <- as.matrix(object$Xi)
  Wd <- if (inherits(object, "FLR")) diag(nrow(Xi)) else as.matrix(base$W)
  Phi <- fda::eval.fd(gp, red$PCAcoef)
  c0 <- effect_multipliers(theta, Xi, Wd)
  Jc <- matrix(0, 3L, d)
  for (j in seq_len(d)) {
    if (inherits(object, "FLR") && j == d) next
    e <- numeric(d); e[j] <- h
    Jc[, j] <- (effect_multipliers(theta + e, Xi, Wd) - effect_multipliers(theta - e, Xi, Wd)) / (2 * h)
  }
  z <- stats::qnorm(1 - (1 - level) / 2)
  b <- as.numeric(Phi %*% theta[1 + seq_len(K)])
  Gb <- cbind(0, Phi, 0)
  one <- function(lab, est, Gr) {
    se <- sqrt(pmax(rowSums((Gr %*% V) * Gr), 0) / n)
    data.frame(curve = lab, t = gp, estimate = est, se = se,
               lower = est - z * se, upper = est + z * se)
  }
  curves <- rbind(one("beta", b, Gb),
                  one("ADE", c0[[1]] * b, outer(b, Jc[1, ]) + c0[[1]] * Gb),
                  one("AIE", c0[[2]] * b, outer(b, Jc[2, ]) + c0[[2]] * Gb),
                  one("ATE", c0[[3]] * b, outer(b, Jc[3, ]) + c0[[3]] * Gb))
  out <- list(curves = curves, multipliers = c0,
              indirect_share = unname(c0[["indirect"]] / c0[["total"]]),
              level = level, type = vcov$type, bandwidth = vcov$bandwidth)
  class(out) <- "fLSL_effects"
  out
}

#' @export
print.fLSL_effects <- function(x, digits = 4, ...) {
  cat("Functional effect curves (", 100 * x$level, "% pointwise bands, ",
      x$type, " covariance)\n", sep = "")
  m <- x$multipliers
  cat(sprintf("  c_D = %s, c_I = %s, c_D + c_I = %s, indirect share = %s\n",
              format(m[["direct"]], digits = digits), format(m[["indirect"]], digits = digits),
              format(m[["total"]], digits = digits), format(x$indirect_share, digits = digits)))
  invisible(x)
}

#' @export
plot.fLSL_effects <- function(x, ...) {
  op <- graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2.2, 1), font.main = 1)
  on.exit(graphics::par(op))
  labs <- c(beta = "beta(t)", ADE = "ADE(t)", AIE = "AIE(t)", ATE = "ATE(t)")
  for (cv in names(labs)) {
    s <- x$curves[x$curves$curve == cv, ]
    graphics::plot(s$t, s$estimate, type = "n", ylim = range(c(s$lower, s$upper)),
                   xlab = "t", ylab = labs[[cv]], main = labs[[cv]], ...)
    graphics::polygon(c(s$t, rev(s$t)), c(s$lower, rev(s$upper)), col = "grey85", border = NA)
    graphics::lines(s$t, s$estimate, lwd = 2)
    graphics::abline(h = 0, lty = 3)
  }
  invisible(x)
}
