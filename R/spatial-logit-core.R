# ==============================================================
# Core probability and Jacobian of the functional spatial-index logit
#
#   p_i(theta) = Lambda( (S_rho m)_i / s_i(rho) ),
#   S_rho = (I - rho W)^{-1},  s_i(rho)^2 = (S_rho S_rho')_ii,
#   m = alpha 1 + Xi b.
#
# Two entry points:
#   .flsl_probability()      probabilities only (one n x n inverse)
#   .flsl_general_gradient() probabilities and the n x (K+2) Jacobian
#                            (one inverse and one dense n x n product)
#
# At rho = 0 both functions avoid any n x n operation.
# W must be the row-standardized dense matrix; Ws is an optional sparse
# copy of the same matrix, used only for the product W %*% S.
# ==============================================================

.elapsed <- function() proc.time()[["elapsed"]]

.check_core_inputs <- function(alpha, b, rho, Xi, W) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      !is.numeric(rho) || length(rho) != 1L || !is.finite(rho))
    stop("'alpha' and 'rho' must be finite numeric scalars.", call. = FALSE)
  if (abs(rho) >= 1)
    stop("The spatial parameter must satisfy |rho| < 1.", call. = FALSE)
  if (length(b) != ncol(Xi) || any(!is.finite(b)))
    stop("Length of 'b' must equal ncol(Xi) and 'b' must be finite.", call. = FALSE)
  if (nrow(W) != nrow(Xi) || ncol(W) != nrow(Xi))
    stop("'W' must be an n x n matrix.", call. = FALSE)
  invisible(TRUE)
}

.spatial_inverse <- function(rho, W) {
  A <- -rho * W
  diag(A) <- diag(A) + 1
  solve(A)
}

#' Probabilities of the spatial-index logit (no Jacobian)
#' @keywords internal
.flsl_probability <- function(alpha, b, rho, Xi, W) {
  Xi <- as.matrix(Xi); b <- as.numeric(b)
  alpha <- as.numeric(alpha); rho <- as.numeric(rho)
  .check_core_inputs(alpha, b, rho, Xi, W)
  m <- as.numeric(alpha + Xi %*% b)
  if (rho == 0) {
    sigma <- rep(1, length(m))
    eta <- m
  } else {
    S <- .spatial_inverse(rho, W)
    sigma <- sqrt(rowSums(S * S))
    eta <- as.numeric(S %*% m) / sigma
  }
  P <- stats::plogis(eta)
  list(P = P, eta = eta, sigma = sigma, m = m)
}

#' Probabilities and analytic Jacobian of the spatial-index logit
#'
#' Column order of G: alpha, b_1, ..., b_K, rho.
#' d s_i / d rho = (T S')_ii / s_i with T = S W S, computed as
#' rowSums(T * S) / s_i, so that S S' is never formed.
#' @keywords internal
.flsl_general_gradient <- function(alpha, b, rho, Xi, W, Ws = NULL) {
  Xi <- as.matrix(Xi); b <- as.numeric(b)
  alpha <- as.numeric(alpha); rho <- as.numeric(rho)
  .check_core_inputs(alpha, b, rho, Xi, W)
  n <- nrow(Xi); K <- ncol(Xi)
  B <- cbind(1, Xi)
  m <- as.numeric(B %*% c(alpha, b))
  Wmul <- function(X) {
    if (is.null(Ws)) W %*% X else as.matrix(Ws %*% X)
  }

  if (rho == 0) {
    sigma <- rep(1, n)
    Sm <- m
    SB <- B
    deta_rho <- as.numeric(Wmul(m))          # T = W, s'(0) = 0
  } else {
    S <- .spatial_inverse(rho, W)
    sigma <- sqrt(rowSums(S * S))
    Sm <- as.numeric(S %*% m)
    SB <- S %*% B
    Tmat <- S %*% Wmul(S)                    # dS/drho = S W S
    dsigma <- rowSums(Tmat * S) / sigma
    dSm <- as.numeric(Tmat %*% m)
    deta_rho <- dSm / sigma - Sm * dsigma / sigma^2
  }

  eta <- Sm / sigma
  P <- stats::plogis(eta)
  q <- P * (1 - P)
  G <- cbind(SB * (q / sigma), q * deta_rho)
  colnames(G) <- c("alpha", paste0("b", seq_len(K)), "rho")

  list(P = P, q = q, eta = eta, m = m, sigma = sigma, Sm = Sm, G = G)
}

#' Sparse copy of a weight matrix (NULL when Matrix is unavailable)
#' @keywords internal
.as_sparse_W <- function(W) {
  if (inherits(W, "sparseMatrix")) return(W)
  if (!requireNamespace("Matrix", quietly = TRUE)) return(NULL)
  Matrix::Matrix(as.matrix(W), sparse = TRUE)
}
