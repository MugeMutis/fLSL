# ==============================================================
# Functional dimension reduction
#
# This file contains:
#   getFPCA()
#
#
# The function returns objects inheriting from
# "flsl_reduction".
# ==============================================================


#' Functional principal component analysis
#'
#' Performs functional principal component analysis (FPCA) after
#' representing the observed curves in a B-spline basis.
#'
#' The B-spline coefficient matrix is centered using the training
#' sample mean function. Let C denote the centered coefficient matrix
#' and
#'
#' \deqn{
#'   J = \int B(t) B(t)^\top dt
#' }
#'
#' denote the B-spline Gram matrix. The transformed coefficient
#' matrix is
#'
#' \deqn{
#'   C J^{1/2}.
#' }
#'
#' Ordinary FPCA is then obtained from the empirical covariance
#' matrix of the transformed coefficients.
#'
#' If \code{ncomp = NULL}, the number of retained components is
#' selected as the smallest number whose cumulative fraction of
#' variance explained is at least \code{threshold}. If
#' \code{ncomp} is supplied, exactly that number of components is
#' retained.
#'
#' @param data Numeric matrix containing functional observations in
#'   rows and grid points in columns.
#' @param nbasis Number of B-spline basis functions.
#' @param gp Numeric vector containing the functional observation grid.
#' @param ncomp Optional fixed number of retained components.
#' @param threshold Fraction-of-variance-explained threshold used when
#'   \code{ncomp = NULL}.
#' @param eig_tol Relative numerical tolerance used when determining
#'   positive eigenvalues.
#'
#' @return An object of classes \code{"flsl_fpca"} and
#'   \code{"flsl_reduction"}.
#'
#' @export
getFPCA <- function(
    data,
    nbasis = 15,
    gp = NULL,
    ncomp = NULL,
    threshold = 0.95,
    eig_tol = 1e-10
) {
  
  ## ============================================================
  ## 1. Input checks
  ## ============================================================
  
  if (!is.matrix(data)) {
    data <- as.matrix(data)
  }
  
  if (!is.numeric(data) ||
      any(!is.finite(data))) {
    
    stop(
      "'data' must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  if (n < 2L) {
    stop(
      "'data' must contain at least two functional observations.",
      call. = FALSE
    )
  }
  
  if (p < 2L) {
    stop(
      "'data' must contain at least two grid points.",
      call. = FALSE
    )
  }
  
  
  if (is.null(gp)) {
    gp <- seq(
      0,
      1,
      length.out = p
    )
  }
  
  gp <- as.numeric(gp)
  
  if (length(gp) != p ||
      any(!is.finite(gp))) {
    
    stop(
      "'gp' must be a finite numeric vector with length equal to ncol(data).",
      call. = FALSE
    )
  }
  
  if (is.unsorted(gp, strictly = TRUE)) {
    stop(
      "'gp' must be strictly increasing.",
      call. = FALSE
    )
  }
  
  
  if (!is.numeric(nbasis) ||
      length(nbasis) != 1L ||
      !is.finite(nbasis) ||
      nbasis != as.integer(nbasis) ||
      nbasis < 4L) {
    
    stop(
      "'nbasis' must be an integer greater than or equal to 4.",
      call. = FALSE
    )
  }
  
  nbasis <- as.integer(nbasis)
  
  
  if (nbasis > p) {
    stop(
      "'nbasis' cannot exceed the number of observed grid points.",
      call. = FALSE
    )
  }
  
  
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold <= 0 ||
      threshold > 1) {
    
    stop(
      "'threshold' must lie in (0, 1].",
      call. = FALSE
    )
  }
  
  
  if (!is.numeric(eig_tol) ||
      length(eig_tol) != 1L ||
      !is.finite(eig_tol) ||
      eig_tol <= 0) {
    
    stop(
      "'eig_tol' must be a single positive finite number.",
      call. = FALSE
    )
  }
  
  
  if (!is.null(ncomp)) {
    
    if (!is.numeric(ncomp) ||
        length(ncomp) != 1L ||
        !is.finite(ncomp) ||
        ncomp != as.integer(ncomp) ||
        ncomp < 1L) {
      
      stop(
        "'ncomp' must be NULL or a positive integer.",
        call. = FALSE
      )
    }
    
    ncomp <- as.integer(ncomp)
  }
  
  
  ## ============================================================
  ## 2. B-spline representation
  ## ============================================================
  
  bs_basis <- fda::create.bspline.basis(
    rangeval = range(gp),
    nbasis = nbasis
  )
  
  
  fd_fit <- fda::smooth.basisPar(
    argvals = gp,
    y = t(data),
    fdobj = bs_basis,
    Lfdobj = NULL,
    lambda = 0
  )
  
  
  fd_object <- fd_fit$fd
  
  
  ## Coefficient matrix:
  ## observations x B-spline basis functions
  
  coef_matrix <- t(
    fd_object$coefs
  )
  
  
  ## ============================================================
  ## 3. Functional mean and centered coefficient matrix
  ## ============================================================
  
  mean_coef <- colMeans(
    coef_matrix
  )
  
  
  centered_coef <- sweep(
    coef_matrix,
    MARGIN = 2L,
    STATS = mean_coef,
    FUN = "-"
  )
  
  
  mean_function <- fda::fd(
    coef = matrix(
      mean_coef,
      ncol = 1L
    ),
    basisobj = bs_basis
  )
  
  
  ## ============================================================
  ## 4. B-spline Gram matrix
  ##
  ## J = integral B(t) B(t)' dt
  ## ============================================================
  
  J <- fda::inprod(
    bs_basis,
    bs_basis
  )
  
  J <- as.matrix(J)
  
  sqrt_J <- .matrix_sqrt_spd(
    J
  )
  
  J_half <- sqrt_J$B
  J_inv_half <- sqrt_J$Binv
  
  
  ## ============================================================
  ## 5. Transform to Euclidean functional geometry
  ##
  ## Z_X = C J^(1/2)
  ## ============================================================
  
  transformed_data <- centered_coef %*%
    J_half
  
  
  ## ============================================================
  ## 6. Ordinary covariance matrix
  ## ============================================================
  
  covariance_matrix <- crossprod(
    transformed_data
  ) / (n - 1)
  
  
  covariance_matrix <- (
    covariance_matrix +
      t(covariance_matrix)
  ) / 2
  
  
  ## ============================================================
  ## 7. Eigendecomposition
  ## ============================================================
  
  ee <- eigen(
    covariance_matrix,
    symmetric = TRUE
  )
  
  
  eigenvalues <- as.numeric(
    ee$values
  )
  
  eigenvectors <- ee$vectors
  
  
  ## Numerical cleanup of tiny negative covariance eigenvalues.
  
  max_abs_eigenvalue <- max(
    abs(eigenvalues)
  )
  
  eigen_cutoff <- eig_tol *
    max(1, max_abs_eigenvalue)
  
  
  eigenvalues[
    eigenvalues < 0 &
      abs(eigenvalues) <= eigen_cutoff
  ] <- 0
  
  
  if (any(eigenvalues < -eigen_cutoff)) {
    stop(
      "The FPCA covariance matrix has a materially negative eigenvalue.",
      call. = FALSE
    )
  }
  
  
  positive_indices <- which(
    eigenvalues > eigen_cutoff
  )
  
  
  if (length(positive_indices) < 1L) {
    stop(
      "No positive FPCA eigenvalues were found.",
      call. = FALSE
    )
  }
  
  
  positive_values <- eigenvalues[
    positive_indices
  ]
  
  
  ## ============================================================
  ## 8. Determine number of components
  ## ============================================================
  
  cumulative_prop <- cumsum(
    positive_values
  ) / sum(
    positive_values
  )
  
  
  if (is.null(ncomp)) {
    
    K <- which(
      cumulative_prop >= threshold
    )[1L]
    
    selection_rule <-
      "fraction_of_variance_explained"
    
  } else {
    
    if (ncomp > length(positive_indices)) {
      
      stop(
        paste0(
          "'ncomp' exceeds the number of numerically positive ",
          "FPCA eigenvalues (",
          length(positive_indices),
          ")."
        ),
        call. = FALSE
      )
    }
    
    K <- ncomp
    
    selection_rule <-
      "fixed_number_of_components"
  }
  
  
  K <- as.integer(
    K
  )
  
  
  selected_indices <- positive_indices[
    seq_len(K)
  ]
  
  
  selected_values <- eigenvalues[
    selected_indices
  ]
  
  
  V <- eigenvectors[
    ,
    selected_indices,
    drop = FALSE
  ]
  
  
  ## ============================================================
  ## 9. Recover L2-orthonormal functional components
  ##
  ## If v_k is Euclidean-normalized, then
  ##
  ## a_k = J^(-1/2) v_k
  ##
  ## and phi_k(t) = B(t)' a_k.
  ## ============================================================
  
  PCAcoef_matrix <- J_inv_half %*%
    V
  
  
  PCAcoef <- fda::fd(
    coef = PCAcoef_matrix,
    basisobj = bs_basis
  )
  
  
  ## ============================================================
  ## 10. Canonical component scores
  ##
  ## Xi = C J A
  ##
  ## which is algebraically identical to
  ##
  ## Xi = C J^(1/2) V.
  ## ============================================================
  
  PCAscore <- centered_coef %*%
    J %*%
    PCAcoef_matrix
  
  
  colnames(PCAscore) <- paste0(
    "PC",
    seq_len(K)
  )
  
  
  ## ============================================================
  ## 11. Internal geometry checks
  ## ============================================================
  
  gram_components <- crossprod(
    PCAcoef_matrix,
    J %*% PCAcoef_matrix
  )
  
  
  orthonormality_error <- max(
    abs(
      gram_components -
        diag(K)
    )
  )
  
  
  score_matrix_check <- transformed_data %*%
    V
  
  
  score_reconstruction_error <- max(
    abs(
      PCAscore -
        score_matrix_check
    )
  )
  
  
  if (!is.finite(orthonormality_error) ||
      orthonormality_error > 1e-7) {
    
    stop(
      paste0(
        "FPCA functional orthonormality check failed. ",
        "Maximum error = ",
        signif(orthonormality_error, 6),
        "."
      ),
      call. = FALSE
    )
  }
  
  
  if (!is.finite(score_reconstruction_error) ||
      score_reconstruction_error > 1e-7) {
    
    stop(
      paste0(
        "FPCA score reconstruction check failed. ",
        "Maximum error = ",
        signif(score_reconstruction_error, 6),
        "."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 12. Return reduction object
  ## ============================================================
  
  out <- list(
    
    method =
      "FPCA",
    
    call =
      match.call(),
    
    ## Data geometry
    
    gp =
      gp,
    
    n =
      n,
    
    p =
      p,
    
    nbasis =
      nbasis,
    
    ## B-spline representation
    
    bs_basis =
      bs_basis,
    
    fd =
      fd_object,
    
    coefficients =
      coef_matrix,
    
    centered_coef =
      centered_coef,
    
    mean_function =
      mean_function,
    
    mean_coef =
      mean_coef,
    
    ## Gram geometry
    
    gram_matrix =
      J,
    
    gram_sqrt =
      J_half,
    
    gram_inv_sqrt =
      J_inv_half,
    
    ## FPCA quantities
    
    covariance_matrix =
      covariance_matrix,
    
    eigenvalues =
      eigenvalues,
    
    positive_eigenvalues =
      positive_values,
    
    selected_eigenvalues =
      selected_values,
    
    selected_indices =
      selected_indices,
    
    eigenvectors =
      V,
    
    PCAcoef =
      PCAcoef,
    
    PCAcoef_matrix =
      PCAcoef_matrix,
    
    PCAscore =
      PCAscore,
    
    ## Dimension selection
    
    n_positive =
      length(positive_indices),
    
    ncomp =
      K,
    
    threshold =
      threshold,
    
    cumulative_prop =
      cumulative_prop,
    
    selection_rule =
      selection_rule,
    
    eig_tol =
      eig_tol,
    
    ## Numerical diagnostics
    
    orthonormality_error =
      orthonormality_error,
    
    score_reconstruction_error =
      score_reconstruction_error
  )
  
  
  class(out) <- c(
    "flsl_fpca",
    "flsl_reduction"
  )
  
  
  out
}



