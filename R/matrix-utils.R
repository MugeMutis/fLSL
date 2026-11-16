#' Symmetric square root of a positive-definite matrix
#'
#' Computes the symmetric square root and inverse square root of a
#' real symmetric positive-definite matrix using an eigendecomposition.
#'
#' This is an internal utility used in the functional dimension-reduction
#' routines to construct the square root of the B-spline Gram matrix.
#'
#' @param A A numeric square matrix.
#' @param tol Numeric relative tolerance used to assess numerical
#'   positive definiteness. Default is
#'   \code{sqrt(.Machine$double.eps)}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{B}}{The symmetric square root of \code{A}.}
#'   \item{\code{Binv}}{The symmetric inverse square root of \code{A}.}
#'   \item{\code{eigenvalues}}{The eigenvalues of the symmetrized matrix.}
#' }
#'
#' @keywords internal
.matrix_sqrt_spd <- function(
    A,
    tol = sqrt(.Machine$double.eps)
) {
  
  ## ============================================================
  ## 1. Input checks
  ## ============================================================
  
  if (!is.matrix(A)) {
    A <- as.matrix(A)
  }
  
  if (!is.numeric(A)) {
    stop(
      "'A' must be a numeric matrix.",
      call. = FALSE
    )
  }
  
  if (nrow(A) != ncol(A)) {
    stop(
      "'A' must be a square matrix.",
      call. = FALSE
    )
  }
  
  if (nrow(A) < 1L) {
    stop(
      "'A' must have positive dimension.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(A))) {
    stop(
      "'A' must contain only finite values.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      !is.finite(tol) ||
      tol <= 0) {
    
    stop(
      "'tol' must be a single positive finite number.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 2. Enforce numerical symmetry
  ## ============================================================
  
  A_sym <- (A + t(A)) / 2
  
  
  ## ============================================================
  ## 3. Eigendecomposition
  ## ============================================================
  
  ee <- eigen(
    A_sym,
    symmetric = TRUE
  )
  
  eigenvalues <- as.numeric(
    ee$values
  )
  
  V <- ee$vectors
  
  
  ## ============================================================
  ## 4. Numerical positive-definiteness check
  ## ============================================================
  
  max_eigenvalue <- max(
    eigenvalues
  )
  
  if (!is.finite(max_eigenvalue) ||
      max_eigenvalue <= 0) {
    
    stop(
      "'A' must be positive definite.",
      call. = FALSE
    )
  }
  
  cutoff <- tol * max_eigenvalue
  
  if (any(eigenvalues <= cutoff)) {
    
    stop(
      paste0(
        "'A' must be numerically positive definite. ",
        "At least one eigenvalue is smaller than or equal to ",
        "the numerical cutoff."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 5. Symmetric square root
  ##
  ## A^(1/2) = V diag(sqrt(lambda)) V'
  ## ============================================================
  
  B <- sweep(
    V,
    MARGIN = 2L,
    STATS = sqrt(eigenvalues),
    FUN = "*"
  ) %*% t(V)
  
  
  ## ============================================================
  ## 6. Symmetric inverse square root
  ##
  ## A^(-1/2) = V diag(1 / sqrt(lambda)) V'
  ## ============================================================
  
  Binv <- sweep(
    V,
    MARGIN = 2L,
    STATS = 1 / sqrt(eigenvalues),
    FUN = "*"
  ) %*% t(V)
  
  
  ## ============================================================
  ## 7. Remove negligible numerical asymmetry
  ## ============================================================
  
  B <- (B + t(B)) / 2
  
  Binv <- (Binv + t(Binv)) / 2
  
  
  ## ============================================================
  ## 8. Return
  ## ============================================================
  
  list(
    B = B,
    Binv = Binv,
    eigenvalues = eigenvalues
  )
}