# ==============================================================
# Non-spatial functional logistic regression
# ==============================================================


#' Functional logistic regression
#'
#' Fits a non-spatial scalar-on-function logistic regression model
#' using functional principal component scores.
#'
#' The fitted model is
#'
#' \deqn{
#'   logit(P_i) =
#'   alpha + Xi_i^\top b,
#' }
#'
#' where Xi contains the retained FPCA scores. The estimated slope
#' function is reconstructed as
#'
#' \deqn{
#'   beta_hat(t) =
#'   sum_k b_hat_k phi_hat_k(t).
#' }
#'
#' This model serves as the non-spatial benchmark for FLSL and
#' IFLSL.
#'
#' @param y Binary response vector.
#' @param x Numeric matrix containing functional observations in rows
#'   and grid points in columns.
#' @param gp Optional numeric vector containing the functional grid.
#' @param nbasis Number of B-spline basis functions.
#' @param ncomp Optional fixed number of retained FPCA components.
#' @param threshold Fraction-of-variance-explained threshold used when
#'   \code{ncomp = NULL}.
#' @param eig_tol Numerical tolerance used by \code{getFPCA()}.
#'
#' @return An object of class \code{"FLR"}.
#'
#' @export
FLR <- function(
    y,
    x,
    gp = NULL,
    nbasis = 15,
    ncomp = NULL,
    threshold = 0.95,
    eig_tol = 1e-10
) {
  
  ## ============================================================
  ## 1. Input checks
  ## ============================================================
  
  if (!is.numeric(y) && !is.logical(y)) {
    stop("'y' must be numeric or logical, coded 0/1.", call. = FALSE)
  }

  y <- as.numeric(
    y
  )
  
  
  if (any(!is.finite(y))) {
    stop(
      "'y' must contain only finite values.",
      call. = FALSE
    )
  }
  
  
  if (!all(y %in% c(0, 1))) {
    stop(
      "'y' must be a binary vector containing only 0 and 1.",
      call. = FALSE
    )
  }
  
  
  if (!is.matrix(x)) {
    x <- as.matrix(
      x
    )
  }
  
  
  if (!is.numeric(x) ||
      any(!is.finite(x))) {
    
    stop(
      "'x' must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  
  
  n <- nrow(
    x
  )
  
  
  if (length(y) != n) {
    
    stop(
      "length(y) must equal nrow(x).",
      call. = FALSE
    )
  }
  
  
  if (length(unique(y)) < 2L) {
    
    stop(
      "'y' must contain both response classes.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 2. FPCA reduction
  ##
  ## FLR is intentionally the ordinary-FPCA non-spatial
  ## benchmark in the proposed simulation design.
  ## ============================================================
  
  red <- getFPCA(
    data = x,
    nbasis = nbasis,
    gp = gp,
    ncomp = ncomp,
    threshold = threshold,
    eig_tol = eig_tol
  )
  
  
  Xi <- as.matrix(
    red$PCAscore
  )
  
  
  K <- ncol(
    Xi
  )
  
  
  ## ============================================================
  ## 3. Logistic regression design
  ## ============================================================
  
  design <- cbind(
    intercept = 1,
    Xi
  )
  
  
  colnames(design) <- c(
    "intercept",
    paste0(
      "b",
      seq_len(K)
    )
  )
  
  
  ## ============================================================
  ## 4. Fit logistic regression
  ##
  ## glm.fit() avoids constructing a temporary model formula and
  ## gives us direct control over the design matrix.
  ## ============================================================
  
  fit <- stats::glm.fit(
    x = design,
    y = y,
    family = stats::binomial(
      link = "logit"
    )
  )
  
  
  if (!isTRUE(fit$converged)) {
    stop("Non-spatial functional logit failed to converge.", call. = FALSE)
  }

  if (fit$rank < ncol(design)) {
    stop("The non-spatial functional logit design is rank deficient.",
         call. = FALSE)
  }

  if (isTRUE(fit$boundary)) {
    warning("The initial logit fit reached a boundary; check for separation.",
            call. = FALSE)
  }

  coefficients <- as.numeric(
    fit$coefficients
  )
  
  
  names(coefficients) <- colnames(
    design
  )
  
  
  if (any(!is.finite(coefficients))) {
    
    stop(
      paste0(
        "The initial functional logistic regression produced ",
        "non-finite coefficient estimates."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 5. Extract estimates
  ## ============================================================
  
  alpha_hat <- unname(
    coefficients["intercept"]
  )
  
  
  b_hat <- unname(
    coefficients[
      paste0(
        "b",
        seq_len(K)
      )
    ]
  )
  
  
  names(b_hat) <- paste0(
    "b",
    seq_len(K)
  )
  
  
  ## ============================================================
  ## 6. Reconstruct beta_hat(t)
  ## ============================================================
  
  gp_fit <- red$gp
  
  
  Phi <- fda::eval.fd(
    gp_fit,
    red$PCAcoef
  )
  
  
  beta_hat <- as.vector(
    Phi %*%
      b_hat
  )
  
  
  ## ============================================================
  ## 7. Training linear predictor and probabilities
  ## ============================================================
  
  eta_hat <- as.vector(
    alpha_hat +
      Xi %*% b_hat
  )
  
  
  probability <- stats::plogis(
    eta_hat
  )
  
  
  residuals <- y -
    probability
  
  
  ## ============================================================
  ## 8. Return fitted object
  ## ============================================================
  
  out <- list(
    
    call =
      match.call(),
    
    method =
      "FLR-FPCA",
    
    ## Parameter estimates
    
    alpha_hat =
      alpha_hat,
    
    b_hat =
      b_hat,
    
    coefficients =
      coefficients,
    
    beta_hat =
      beta_hat,
    
    ## Training fitted quantities
    
    eta =
      eta_hat,
    
    probability =
      probability,
    
    fitted.values =
      probability,
    
    residuals =
      residuals,
    
    ## Functional reduction
    
    reduction_method =
      "FPCA",
    
    reduction =
      red,
    
    Xi =
      Xi,
    
    ncomp =
      K,
    
    gp =
      gp_fit,
    
    ## Logistic fit
    
    initial_fit =
      fit,
    
    ## Settings
    
    nbasis =
      nbasis,
    
    threshold =
      threshold,
    
    eig_tol =
      eig_tol
  )
  
  
  class(out) <- "FLR"
  
  
  out
}