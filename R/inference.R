# ==============================================================
# Inference for theta (Sections 3.4-3.5)
#
#   hac_omega()  : spatial HAC estimate (3.11) of the moment covariance,
#                  Bartlett kernel, sparse sum over pairs within the
#                  bandwidth (no n x n matrix is formed), negative
#                  eigenvalues set to zero.
#   vcov_theta() : covariance V of sqrt(n)(theta_hat - theta), under
#                  conditional independence (3.9) or with the HAC Omega:
#                    FLSL       V^(0) at theta^(0) (Corollary 3.1);
#                    IFLSL-LS,  projected-GMM sandwich at theta_hat;
#                    NLGMM
#                    QMLE       Bernoulli sandwich (score contributions);
#                    NLS        NLS sandwich;
#                    FLR        logit information (rho fixed at 0).
#   wald_rho()   : Wald test of H0: rho = 0 from vcov_theta().
# ==============================================================

#' Spatial HAC estimate of a moment covariance
#'
#' Computes \eqn{n^{-1}\sum_i\sum_j k(d_{ij}/\bar d)\, U_i U_j^\top} with the
#' Bartlett kernel \eqn{k(x) = (1 - x)_+}, equation (3.11) of the paper, where
#' the rows \eqn{U_i} are moment contributions such as \eqn{\hat u_i z_i}.
#' Only pairs within the bandwidth enter, so the estimate is a sparse sum.
#' Negative eigenvalues of the finite-sample estimate are set to zero.
#'
#' @param U n x q matrix of moment contributions.
#' @param coords n x 2 matrix of coordinates in the units of \code{bandwidth}
#'   (for geographic data, project to kilometres or metres first).
#' @param bandwidth Bandwidth \eqn{\bar d > 0}.
#' @param chunk Number of query points processed at a time.
#' @return The q x q matrix, with attribute \code{"mean_neighbours"}.
#' @details The neighbour search uses \code{dbscan::frNN} when the
#'   \pkg{dbscan} package is installed and blockwise distances otherwise.
#' @export
hac_omega <- function(U, coords, bandwidth, chunk = 1000L) {
  U <- as.matrix(U); coords <- as.matrix(coords)
  n <- nrow(U)
  if (!is.numeric(U) || any(!is.finite(U))) stop("'U' must be a finite numeric matrix.", call. = FALSE)
  if (nrow(coords) != n || !is.numeric(coords) || any(!is.finite(coords)))
    stop("'coords' must be a finite numeric matrix with nrow(U) rows.", call. = FALSE)
  if (length(bandwidth) != 1L || !is.finite(bandwidth) || bandwidth <= 0)
    stop("'bandwidth' must be a positive number.", call. = FALSE)
  use_dbscan <- requireNamespace("dbscan", quietly = TRUE)
  Om <- crossprod(U); n_pairs <- 0
  for (s in seq(1L, n, by = chunk)) {
    rows <- s:min(n, s + chunk - 1L)
    if (use_dbscan) {
      nn <- dbscan::frNN(coords, eps = bandwidth, query = coords[rows, , drop = FALSE],
                         sort = FALSE)
      len <- lengths(nn$id)
      i <- rep(seq_along(rows), len)
      j <- unlist(nn$id, use.names = FALSE)
      dist <- unlist(nn$dist, use.names = FALSE)
    } else {
      D2 <- matrix(0, length(rows), n)
      for (k in seq_len(ncol(coords)))
        D2 <- D2 + outer(coords[rows, k], coords[, k], "-")^2
      idx <- which(D2 < bandwidth^2, arr.ind = TRUE)
      i <- idx[, 1L]; j <- idx[, 2L]; dist <- sqrt(D2[idx])
    }
    keep <- rows[i] != j & dist < bandwidth
    if (any(keep)) {
      n_pairs <- n_pairs + sum(keep)
      Kc <- Matrix::sparseMatrix(i = i[keep], j = j[keep], x = 1 - dist[keep] / bandwidth,
                                 dims = c(length(rows), n))
      Om <- Om + crossprod(U[rows, , drop = FALSE], as.matrix(Kc %*% U))
    }
  }
  Om <- (Om + t(Om)) / 2 / n
  e <- eigen(Om, symmetric = TRUE)
  Om <- e$vectors %*% diag(pmax(e$values, 0), length(e$values)) %*% t(e$vectors)
  dimnames(Om) <- list(colnames(U), colnames(U))
  attr(Om, "mean_neighbours") <- n_pairs / n
  Om
}

#' Parameter vector of a fit
#' @keywords internal
.theta_of <- function(object) {
  th <- c(object$alpha_hat, as.numeric(object$b_hat), object$rho_hat)
  names(th) <- c("alpha", paste0("b", seq_along(object$b_hat)), "rho")
  th
}

#' Training coordinates for HAC
#' @keywords internal
.train_coords <- function(coords, n, tr) {
  if (is.null(coords)) stop("'coords' must be supplied for HAC inference.", call. = FALSE)
  coords <- as.matrix(coords)
  if (nrow(coords) == n) return(coords[tr, , drop = FALSE])
  if (nrow(coords) == length(tr)) return(coords)
  stop("'coords' must have one row per unit (or per training unit).", call. = FALSE)
}

#' Covariance of the estimated parameters
#'
#' Returns the asymptotic covariance V of \eqn{\sqrt{n}(\hat\theta - \theta)}
#' for a fitted model, estimated under conditional independence (Theorem
#' 3.1(iii)) or with the spatial HAC moment covariance (3.11). For an FLSL
#' fit, V is the inversion-free \eqn{\hat V^{(0)}} evaluated at the nonspatial
#' logit fit (Corollary 3.1), so that it is available without forming the
#' spatial multiplier; with \code{type = "hac"} the residuals of the
#' nonspatial fit are used.
#'
#' @param object A fitted \code{"FLSL"}, \code{"IFLSL_LS"}, \code{"NLGMM"},
#'   \code{"fLSL_benchmark"} or \code{"FLR"} object.
#' @param type \code{"independent"} (default) or \code{"hac"}.
#' @param coords Coordinates of the units, required for \code{type = "hac"}.
#' @param bandwidth HAC bandwidth, in the units of \code{coords}.
#' @return A list of class \code{"fLSL_vcov"} with \code{V}, the standard
#'   errors \code{se = sqrt(diag(V)/n)}, \code{theta}, \code{n} (number of
#'   training units), \code{type} and \code{bandwidth}.
#' @export
vcov_theta <- function(object, type = c("independent", "hac"),
                       coords = NULL, bandwidth = NULL) {
  type <- match.arg(type)
  hac <- identical(type, "hac")
  if (inherits(object, "FLR")) {
    Xi <- as.matrix(object$Xi); K <- ncol(Xi); n <- nrow(Xi); tr <- seq_len(n)
    B <- cbind(1, Xi); p <- object$probability; w <- pmax(p * (1 - p), 1e-10)
    u <- object$residuals
    Hi <- solve(crossprod(B * sqrt(w)) / n)
    Vd <- if (hac) Hi %*% hac_omega(B * u, .train_coords(coords, n, tr), bandwidth) %*% Hi else Hi
    V <- matrix(0, K + 2L, K + 2L); V[seq_len(K + 1L), seq_len(K + 1L)] <- Vd
    theta <- c(object$coefficients, rho = 0)
  } else if (inherits(object, "FLSL")) {
    tr <- object$train; n <- length(tr); K <- object$ncomp
    Z <- object$Z[tr, , drop = FALSE]; G <- object$G[tr, , drop = FALSE]
    Om <- if (hac) hac_omega(Z * object$u0, .train_coords(coords, object$n, tr), bandwidth) else NULL
    V <- .flsl_wald(Z = Z, G0 = G, q0 = object$q0[tr], rho_hat = object$rho_hat, Omega = Om)$V
    theta <- object$theta_hat
  } else if (inherits(object, c("IFLSL_LS", "NLGMM", "fLSL_benchmark"))) {
    if (is.null(object$theta_hat)) stop("The fit has no parameter estimate.", call. = FALSE)
    base <- object$flsl_fit; tr <- object$train; n <- length(tr); K <- object$ncomp
    theta <- object$theta_hat
    if (abs(theta[["rho"]]) >= 1) stop("The estimate is inadmissible (|rho| >= 1).", call. = FALSE)
    o <- .flsl_general_gradient(theta[1], theta[1 + seq_len(K)], theta[K + 2L],
                                as.matrix(base$Xi), as.matrix(base$W), .as_sparse_W(base$W))
    G <- o$G[tr, , drop = FALSE]; P <- o$P[tr]; u <- object$y[tr] - P
    w <- pmax(P * (1 - P), 1e-10)
    xy <- if (hac) .train_coords(coords, object$n, tr) else NULL
    est <- if (inherits(object, "fLSL_benchmark")) object$method else "GMM"
    V <- switch(est,
      GMM = {
        Z <- base$Z[tr, , drop = FALSE]
        A <- solve(crossprod(Z) / n); J <- -crossprod(Z, G) / n
        Om <- if (hac) hac_omega(Z * u, xy, bandwidth) else crossprod(Z * w, Z) / n
        Mi <- solve(crossprod(J, A %*% J))
        Mi %*% t(J) %*% A %*% Om %*% A %*% J %*% Mi
      },
      QMLE = {
        H <- crossprod(G / sqrt(w)) / n; Hi <- solve(H)
        if (hac) Hi %*% hac_omega(G * (u / w), xy, bandwidth) %*% Hi else Hi
      },
      NLS = {
        Ai <- solve(crossprod(G) / n)
        Om <- if (hac) hac_omega(G * u, xy, bandwidth) else crossprod(G * sqrt(w)) / n
        Ai %*% Om %*% Ai
      })
  } else stop("Unsupported object class.", call. = FALSE)
  V <- (V + t(V)) / 2
  nm <- c("alpha", paste0("b", seq_len(nrow(V) - 2L)), "rho")
  dimnames(V) <- list(nm, nm)
  out <- list(V = V, se = stats::setNames(sqrt(pmax(diag(V), 0) / n), nm),
              theta = stats::setNames(as.numeric(theta), nm), n = n, type = type,
              bandwidth = if (hac) bandwidth else NULL)
  class(out) <- "fLSL_vcov"
  out
}

#' Wald test of no spatial propagation
#'
#' Tests \eqn{H_0: \rho = 0} by \eqn{W_n = n\hat\rho^2/\hat V_{\rho\rho}}.
#' For an FLSL fit this is the inversion-free test of Corollary 3.1; with
#' \code{type = "hac"} it uses the spatial HAC covariance (Section 3.5).
#'
#' @inheritParams vcov_theta
#' @param level Nominal level of the test (default 0.05).
#' @return A list with the statistic, p-value, decision, \code{rho_hat},
#'   its standard error and the covariance object.
#' @export
wald_rho <- function(object, type = c("independent", "hac"), coords = NULL,
                     bandwidth = NULL, level = 0.05) {
  type <- match.arg(type)
  if (inherits(object, "FLR")) stop("FLR has no spatial parameter.", call. = FALSE)
  vc <- vcov_theta(object, type = type, coords = coords, bandwidth = bandwidth)
  rho <- unname(vc$theta["rho"]); se <- unname(vc$se["rho"])
  stat <- if (is.finite(se) && se > 0) (rho / se)^2 else NA_real_
  list(statistic = stat,
       p_value = if (is.finite(stat)) stats::pchisq(stat, 1, lower.tail = FALSE) else NA_real_,
       reject = if (is.finite(stat)) stat > stats::qchisq(1 - level, 1) else NA,
       level = level, rho_hat = rho, se_rho = se, type = type,
       bandwidth = vc$bandwidth, vcov = vc)
}
