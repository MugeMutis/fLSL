# ==============================================================
# Functional-data utilities
#
# This file contains utilities used to project new functional
# observations onto an already fitted FPCA representation.
# ==============================================================


#' Project functional observations onto a fitted reduction basis
#'
#' Projects new functional observations onto the component functions
#' obtained from a fitted \code{getFPCA()}
#' object.
#'
#' The projection uses the B-spline basis, functional mean, Gram
#' matrix, and component functions estimated from the training data.
#' No dimension-reduction model is re-estimated for the new data.
#'
#' If C_new denotes the B-spline coefficient matrix of the new
#' curves, mu_c denotes the training mean coefficient vector, J is
#' the training B-spline Gram matrix, and A contains the B-spline
#' coefficients of the retained component functions, then the
#' projected scores are
#'
#' \deqn{
#'   Xi_new = (C_new - mu_c) J A.
#' }
#'
#' @param object A fitted functional reduction object returned by
#'   \code{getFPCA()}.
#' @param newdata Numeric matrix containing new functional
#'   observations in rows and grid points in columns.
#' @param gp Optional numeric vector containing the observation grid
#'   for \code{newdata}. If omitted, the training grid stored in
#'   \code{object} is used.
#'
#' @return A numeric matrix containing projected component scores.
#'
#' @keywords internal
project_scores <- function(
    object,
    newdata,
    gp = NULL
) {
  
  ## ============================================================
  ## 1. Reduction-object checks
  ## ============================================================
  
  if (!inherits(
    object,
    "flsl_reduction"
  )) {
    
    stop(
      paste0(
        "'object' must be a functional reduction object returned ",
        "by getFPCA()."
      ),
      call. = FALSE
    )
  }
  
  
  required_components <- c(
    "gp",
    "bs_basis",
    "mean_coef",
    "gram_matrix",
    "PCAcoef_matrix",
    "ncomp"
  )
  
  
  missing_components <- setdiff(
    required_components,
    names(object)
  )
  
  
  if (length(missing_components) > 0L) {
    
    stop(
      paste0(
        "The reduction object is missing required component(s): ",
        paste(
          missing_components,
          collapse = ", "
        ),
        "."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 2. New functional-data checks
  ## ============================================================
  
  if (!is.matrix(newdata)) {
    newdata <- as.matrix(
      newdata
    )
  }
  
  
  if (!is.numeric(newdata) ||
      any(!is.finite(newdata))) {
    
    stop(
      "'newdata' must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  
  
  if (nrow(newdata) < 1L) {
    
    stop(
      "'newdata' must contain at least one functional observation.",
      call. = FALSE
    )
  }
  
  
  if (ncol(newdata) < 2L) {
    
    stop(
      "'newdata' must contain at least two grid points.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 3. Observation grid
  ## ============================================================
  
  if (is.null(gp)) {
    
    gp_new <- object$gp
    
  } else {
    
    gp_new <- as.numeric(
      gp
    )
  }
  
  
  if (length(gp_new) != ncol(newdata) ||
      any(!is.finite(gp_new))) {
    
    stop(
      paste0(
        "'gp' must be a finite numeric vector with length equal ",
        "to ncol(newdata)."
      ),
      call. = FALSE
    )
  }
  
  
  if (is.unsorted(
    gp_new,
    strictly = TRUE
  )) {
    
    stop(
      "'gp' must be strictly increasing.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 4. Check compatibility with the training B-spline basis
  ##
  ## New curves may in principle be observed on a different grid,
  ## provided that the grid lies inside the domain of the training
  ## basis. The training basis itself is never re-estimated.
  ## ============================================================
  
  training_range <- object$bs_basis$rangeval
  
  
  grid_tol <- sqrt(
    .Machine$double.eps
  )
  
  
  if (min(gp_new) < training_range[1L] - grid_tol ||
      max(gp_new) > training_range[2L] + grid_tol) {
    
    stop(
      paste0(
        "The grid supplied for 'newdata' lies outside the domain ",
        "of the training functional basis."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 5. Smooth new curves using the TRAINING B-spline basis
  ##
  ## This is important:
  ##
  ## - no new basis is constructed;
  ## - no new FPCA is estimated;
  ## - the training representation is preserved.
  ## ============================================================
  
  fd_new_fit <- fda::smooth.basisPar(
    argvals = gp_new,
    y = t(newdata),
    fdobj = object$bs_basis,
    Lfdobj = NULL,
    lambda = 0
  )
  
  
  fd_new <- fd_new_fit$fd
  
  
  ## observations x B-spline coefficients
  
  coef_new <- t(
    fd_new$coefs
  )
  
  
  ## ============================================================
  ## 6. Compatibility checks
  ## ============================================================
  
  if (ncol(coef_new) != length(object$mean_coef)) {
    
    stop(
      paste0(
        "The coefficient dimension of 'newdata' is incompatible ",
        "with the fitted functional reduction object."
      ),
      call. = FALSE
    )
  }
  
  
  J <- as.matrix(
    object$gram_matrix
  )
  
  
  A <- as.matrix(
    object$PCAcoef_matrix
  )
  
  
  if (nrow(J) != ncol(J) ||
      nrow(J) != ncol(coef_new)) {
    
    stop(
      "The stored Gram matrix has incompatible dimensions.",
      call. = FALSE
    )
  }
  
  
  if (nrow(A) != ncol(coef_new)) {
    
    stop(
      "The stored component coefficient matrix has incompatible dimensions.",
      call. = FALSE
    )
  }
  
  
  if (ncol(A) != object$ncomp) {
    
    stop(
      "The stored component dimension is internally inconsistent.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 7. Center using the TRAINING functional mean
  ##
  ## C_new,c = C_new - mean_training
  ## ============================================================
  
  centered_coef_new <- sweep(
    coef_new,
    MARGIN = 2L,
    STATS = object$mean_coef,
    FUN = "-"
  )
  
  
  ## ============================================================
  ## 8. Canonical functional projection
  ##
  ## Xi_new = C_new,c J A
  ##
  ## This is the coefficient/Gram representation of
  ##
  ## < X_new - mu_train, phi_k >_{L2}.
  ##
  ## It is used for both ordinary FPCA .
  ## ============================================================
  
  scores <- centered_coef_new %*%
    J %*%
    A
  
  
  ## ============================================================
  ## 9. Column names
  ## ============================================================
  
  ## Preserve exactly the score names used during training.
  ## In particular, ordinary FPCA uses PC1, PC2, ... (not SPC1, ...).
  training_names <- colnames(object$PCAscore)
  if (!is.null(training_names) &&
      length(training_names) == object$ncomp) {
    colnames(scores) <- training_names
  } else {
    colnames(scores) <- paste0("PC", seq_len(object$ncomp))
  }


  ## ============================================================
  ## 10. Numerical check
  ## ============================================================
  
  if (any(!is.finite(scores))) {
    
    stop(
      "Projection of the new functional observations produced non-finite scores.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 11. Return
  ## ============================================================
  
  scores
}