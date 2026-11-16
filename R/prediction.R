# ==============================================================
# Prediction methods
# With newdata (or scores) and Wnew, predictions are for a test system with
# its own weight matrix Wnew; in the simulation design of Section 4.1 the test
# sample has the same locations, so Wnew = W. Without newdata and scores, the
# fitted probabilities of all units of the design are returned, including
# units whose responses were withheld through 'train' (Section 5).
# ==============================================================


# ==============================================================
# Internal spatial-logit prediction engine
# ==============================================================

#' Spatial-logit prediction engine
#'
#' Internal helper used by the FLSL and IFLSL prediction methods.
#'
#' @param alpha_hat Estimated intercept.
#' @param b_hat Estimated functional score coefficients.
#' @param rho_hat Estimated spatial autoregressive parameter.
#' @param Xi_new Functional score matrix for the prediction sample.
#' @param Wnew Spatial weight matrix for the prediction sample.
#' @param cutoff Probability cutoff for binary classification.
#'
#' @return A list containing probabilities, classes, spatial
#'   indices, and intermediate prediction quantities.
#'
#' @keywords internal
.predict_spatial_logit <- function(
    alpha_hat,
    b_hat,
    rho_hat,
    Xi_new,
    Wnew,
    cutoff = 0.5
) {
  
  ## ------------------------------------------------------------
  ## 1. Functional scores
  ## ------------------------------------------------------------
  
  Xi_new <- as.matrix(Xi_new)
  
  if (!is.numeric(Xi_new) ||
      any(!is.finite(Xi_new))) {
    stop(
      "'Xi_new' must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  
  if (nrow(Xi_new) < 1L) {
    stop(
      "The prediction sample must contain at least one observation.",
      call. = FALSE
    )
  }
  
  b_hat <- as.numeric(b_hat)
  
  if (length(b_hat) != ncol(Xi_new) ||
      any(!is.finite(b_hat))) {
    stop(
      paste0(
        "Length of 'b_hat' must equal the number of ",
        "columns of 'Xi_new'."
      ),
      call. = FALSE
    )
  }
  
  if (!is.numeric(alpha_hat) ||
      length(alpha_hat) != 1L ||
      !is.finite(alpha_hat)) {
    stop(
      "'alpha_hat' must be a finite scalar.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(rho_hat) ||
      length(rho_hat) != 1L ||
      !is.finite(rho_hat)) {
    stop(
      "'rho_hat' must be a finite scalar.",
      call. = FALSE
    )
  }
  
  n_new <- nrow(Xi_new)
  
  
  ## ------------------------------------------------------------
  ## 2. Spatial weight matrix
  ## ------------------------------------------------------------
  
  if (!is.matrix(Wnew)) {
    Wnew <- as.matrix(Wnew)
  }
  
  if (!is.numeric(Wnew) ||
      any(!is.finite(Wnew))) {
    stop(
      "'Wnew' must be a finite numeric matrix.",
      call. = FALSE
    )
  }
  
  if (nrow(Wnew) != n_new ||
      ncol(Wnew) != n_new) {
    stop(
      paste0(
        "'Wnew' must be a ",
        n_new,
        " x ",
        n_new,
        " matrix."
      ),
      call. = FALSE
    )
  }
  
  Wnew <- normalize_weights(Wnew)
  
  
  ## ------------------------------------------------------------
  ## 3. Classification cutoff
  ## ------------------------------------------------------------
  
  if (!is.numeric(cutoff) ||
      length(cutoff) != 1L ||
      !is.finite(cutoff) ||
      cutoff <= 0 ||
      cutoff >= 1) {
    stop(
      "'cutoff' must be a single finite number in (0, 1).",
      call. = FALSE
    )
  }
  
  
  ## ------------------------------------------------------------
  ## 4. Non-spatial systematic component
  ##
  ## m_hat = alpha_hat 1 + Xi_new b_hat
  ## ------------------------------------------------------------
  
  m_hat <- as.vector(
    alpha_hat +
      Xi_new %*% b_hat
  )
  
  
  ## ------------------------------------------------------------
  ## 5. Spatial reduced-form operator
  ##
  ## S_rho = (I - rho_hat Wnew)^(-1)
  ## ------------------------------------------------------------
  
  if (abs(rho_hat) >= 1) {
    stop("Spatial prediction requires |rho_hat| < 1; the one-step estimate is outside the admissible parameter range.",
         call. = FALSE)
  }

  A_rho <- diag(n_new) -
    rho_hat * Wnew
  
  if (rcond(A_rho) <
      sqrt(.Machine$double.eps)) {
    stop(
      "'I - rho_hat * Wnew' is numerically singular.",
      call. = FALSE
    )
  }
  
  S_rho <- solve(A_rho)
  
  
  ## ------------------------------------------------------------
  ## 6. Reduced-form covariance
  ##
  ## Sigma_rho = S_rho S_rho'
  ##
  ## sigma_i =
  ## sqrt{ diag(Sigma_rho)_i }
  ## ------------------------------------------------------------
  
  Sigma_rho <- tcrossprod(
    S_rho
  )
  
  sigma2 <- diag(
    Sigma_rho
  )
  
  if (any(!is.finite(sigma2)) ||
      any(sigma2 <= 0)) {
    stop(
      "Invalid reduced-form spatial variances were obtained.",
      call. = FALSE
    )
  }
  
  sigma <- sqrt(
    sigma2
  )
  
  
  ## ------------------------------------------------------------
  ## 7. Reduced-form systematic component
  ##
  ## S_rho m_hat
  ## ------------------------------------------------------------
  
  spatial_systematic <- as.vector(
    S_rho %*% m_hat
  )
  
  
  ## ------------------------------------------------------------
  ## 8. Standardized spatial index
  ##
  ## eta_hat_i =
  ## (S_rho m_hat)_i / sigma_i
  ## ------------------------------------------------------------
  
  eta_hat <- spatial_systematic /
    sigma
  
  
  ## ------------------------------------------------------------
  ## 9. Predicted probabilities
  ## ------------------------------------------------------------
  
  probability <- stats::plogis(
    eta_hat
  )
  
  
  ## ------------------------------------------------------------
  ## 10. Binary predictions
  ## ------------------------------------------------------------
  
  class_prediction <- as.integer(
    probability >= cutoff
  )
  
  
  ## ------------------------------------------------------------
  ## 11. Return
  ## ------------------------------------------------------------
  
  list(
    probability = probability,
    class = class_prediction,
    eta = eta_hat,
    spatial_systematic = spatial_systematic,
    m = m_hat,
    sigma = sigma,
    scores = Xi_new,
    rho_hat = rho_hat,
    cutoff = cutoff
  )
}



# ==============================================================
# FLR prediction
# ==============================================================

#' Predict from a fitted FLR model
#'
#' Generates predictions from a fitted functional logistic
#' regression model using the functional representation estimated
#' from the training sample.
#'
#' @param object A fitted object of class \code{"FLR"}.
#' @param newdata Optional numeric matrix of new functional
#'   observations, with observations in rows and grid points in
#'   columns.
#' @param scores Optional matrix of precomputed functional scores.
#' @param type Type of prediction. One of \code{"probability"},
#'   \code{"class"}, \code{"link"}, or \code{"all"}.
#' @param cutoff Probability cutoff used for binary classification.
#' @param ... Additional arguments currently unused.
#'
#' @return Predicted probabilities, binary classes, linear
#'   predictors, or all prediction quantities.
#'
#' @export
predict.FLR <- function(
    object,
    newdata = NULL,
    scores = NULL,
    type = c(
      "probability",
      "class",
      "link",
      "all"
    ),
    cutoff = 0.5,
    ...
) {
  
  ## ------------------------------------------------------------
  ## 1. Check fitted object
  ## ------------------------------------------------------------
  
  if (!inherits(object, "FLR")) {
    stop(
      "'object' must be a fitted FLR model.",
      call. = FALSE
    )
  }
  
  type <- match.arg(type)
  
  
  ## ------------------------------------------------------------
  ## 2. Obtain functional scores
  ## ------------------------------------------------------------
  
  if (is.null(newdata) &&
      is.null(scores)) {
    stop(
      "Either 'newdata' or 'scores' must be supplied.",
      call. = FALSE
    )
  }
  
  if (!is.null(newdata) &&
      !is.null(scores)) {
    stop(
      "Supply either 'newdata' or 'scores', but not both.",
      call. = FALSE
    )
  }
  
  if (!is.null(newdata)) {
    
    Xi_new <- project_scores(
      object = object$reduction,
      newdata = newdata
    )
    
  } else {
    
    if (!is.matrix(scores)) {
      scores <- as.matrix(scores)
    }
    
    if (!is.numeric(scores) ||
        any(!is.finite(scores))) {
      stop(
        "'scores' must be a finite numeric matrix.",
        call. = FALSE
      )
    }
    
    Xi_new <- scores
  }
  
  Xi_new <- as.matrix(Xi_new)
  
  
  ## ------------------------------------------------------------
  ## 3. Score dimensions
  ## ------------------------------------------------------------
  
  if (nrow(Xi_new) < 1L) {
    stop(
      "The prediction sample must contain at least one observation.",
      call. = FALSE
    )
  }
  
  if (ncol(Xi_new) != object$ncomp) {
    stop(
      paste0(
        "The prediction score matrix must contain ",
        object$ncomp,
        " columns."
      ),
      call. = FALSE
    )
  }
  
  
  ## ------------------------------------------------------------
  ## 4. Non-spatial linear predictor
  ##
  ## eta_hat =
  ## alpha_hat + Xi_new b_hat
  ## ------------------------------------------------------------
  
  eta_hat <- as.vector(
    object$alpha_hat +
      Xi_new %*% object$b_hat
  )
  
  
  ## ------------------------------------------------------------
  ## 5. Probability
  ## ------------------------------------------------------------
  
  probability <- stats::plogis(
    eta_hat
  )
  
  
  ## ------------------------------------------------------------
  ## 6. Classification cutoff
  ## ------------------------------------------------------------
  
  if (!is.numeric(cutoff) ||
      length(cutoff) != 1L ||
      !is.finite(cutoff) ||
      cutoff <= 0 ||
      cutoff >= 1) {
    stop(
      "'cutoff' must be a single finite number in (0, 1).",
      call. = FALSE
    )
  }
  
  class_hat <- as.integer(
    probability >= cutoff
  )
  
  
  ## ------------------------------------------------------------
  ## 7. Return requested type
  ## ------------------------------------------------------------
  
  if (type == "probability") {
    return(probability)
  }
  
  if (type == "class") {
    return(class_hat)
  }
  
  if (type == "link") {
    return(eta_hat)
  }
  
  list(
    probability = probability,
    class = class_hat,
    link = eta_hat,
    scores = Xi_new,
    cutoff = cutoff
  )
}



# ==============================================================
# FLSL prediction
# ==============================================================

#' Predict from a fitted FLSL model
#'
#' Generates spatial binary-response predictions from a fitted
#' functional linearized spatial logit model.
#'
#' Prediction uses the exact nonlinear reduced-form probability
#' evaluated at the estimated FLSL parameters.
#'
#' @param object A fitted object of class \code{"FLSL"}.
#' @param newdata Optional numeric matrix of new functional
#'   observations.
#' @param scores Optional matrix of precomputed functional scores.
#' @param Wnew Spatial weight matrix for the prediction sample. When
#'   \code{newdata} and \code{scores} are both \code{NULL}, \code{Wnew} is
#'   not used and the fitted probabilities of all units are returned.
#' @param type Type of prediction. One of \code{"probability"},
#'   \code{"class"}, \code{"link"}, or \code{"all"}.
#' @param cutoff Probability cutoff used for binary classification.
#' @param ... Additional arguments currently unused.
#'
#' @return Predicted probabilities, binary classes, standardized
#'   spatial indices, or all prediction quantities.
#'
#' @export
predict.FLSL <- function(
    object,
    newdata = NULL,
    scores = NULL,
    Wnew = NULL,
    type = c(
      "probability",
      "class",
      "link",
      "all"
    ),
    cutoff = 0.5,
    ...
) {
  
  ## ------------------------------------------------------------
  ## 1. Check fitted object
  ## ------------------------------------------------------------
  
  if (!inherits(object, "FLSL")) {
    stop(
      "'object' must be a fitted FLSL model.",
      call. = FALSE
    )
  }
  
  type <- match.arg(type)
  
  
  ## ------------------------------------------------------------
  ## 2. Obtain functional scores
  ## ------------------------------------------------------------
  
  in_sample <- is.null(newdata) && is.null(scores)
  if (in_sample) {
    scores <- object$Xi
    Wnew <- object$W
  }
  
  if (!is.null(newdata) &&
      !is.null(scores)) {
    stop(
      "Supply either 'newdata' or 'scores', but not both.",
      call. = FALSE
    )
  }
  
  if (!is.null(newdata)) {
    
    if (is.null(object$reduction)) {
      stop(
        "The training dimension-reduction object is missing.",
        call. = FALSE
      )
    }
    
    Xi_new <- project_scores(
      object = object$reduction,
      newdata = newdata
    )
    
  } else {
    
    if (!is.matrix(scores)) {
      scores <- as.matrix(scores)
    }
    
    if (!is.numeric(scores) ||
        any(!is.finite(scores))) {
      stop(
        "'scores' must be a finite numeric matrix.",
        call. = FALSE
      )
    }
    
    Xi_new <- scores
  }
  
  Xi_new <- as.matrix(Xi_new)
  
  
  ## ------------------------------------------------------------
  ## 3. Check score dimensions
  ## ------------------------------------------------------------
  
  if (nrow(Xi_new) < 1L) {
    stop(
      "The prediction sample must contain at least one observation.",
      call. = FALSE
    )
  }
  
  if (ncol(Xi_new) != object$ncomp) {
    stop(
      paste0(
        "The prediction score matrix must contain ",
        object$ncomp,
        " columns."
      ),
      call. = FALSE
    )
  }
  
  
  ## ------------------------------------------------------------
  ## 4. Wnew is required
  ## ------------------------------------------------------------
  
  if (is.null(Wnew)) {
    stop(
      "'Wnew' must be supplied with 'newdata' or 'scores'.",
      call. = FALSE
    )
  }
  
  
  ## ------------------------------------------------------------
  ## 5. Common nonlinear spatial prediction
  ## ------------------------------------------------------------
  
  pred <- .predict_spatial_logit(
    alpha_hat = object$alpha_hat,
    b_hat = object$b_hat,
    rho_hat = object$rho_hat,
    Xi_new = Xi_new,
    Wnew = Wnew,
    cutoff = cutoff
  )
  
  
  ## ------------------------------------------------------------
  ## 6. Return requested type
  ## ------------------------------------------------------------
  
  if (type == "probability") {
    return(pred$probability)
  }
  
  if (type == "class") {
    return(pred$class)
  }
  
  if (type == "link") {
    return(pred$eta)
  }
  
  pred
}
