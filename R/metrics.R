# ==============================================================
# PERFORMANCE METRICS
# For probabilistic binary predictions, use predicted probabilities (not classes)
# for Brier score, log loss, RMSE, and MAE.
# Report convergence/invalid-prediction rates separately.
# ==============================================================


# --------------------------------------------------------------
# Internal validation helper
# --------------------------------------------------------------

.check_metric_vectors <- function(
    observed,
    predicted
) {
  
  if (!is.numeric(observed) ||
      !is.numeric(predicted) ||
      !is.null(dim(observed)) ||
      !is.null(dim(predicted))) {
    stop(
      "'observed' and 'predicted' must be numeric.",
      call. = FALSE
    )
  }
  
  if (length(observed) != length(predicted)) {
    stop(
      paste0(
        "'observed' and 'predicted' must have ",
        "the same length."
      ),
      call. = FALSE
    )
  }
  
  if (length(observed) == 0L) {
    stop(
      "Metric vectors must not be empty.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(observed)) ||
      any(!is.finite(predicted))) {
    stop(
      "Metric vectors must contain only finite values.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


# ==============================================================
# BETA ISE
# ==============================================================

#' Integrated squared error for the functional coefficient
#'
#' Computes
#'
#' \deqn{
#'   \int [\hat{\beta}(t)-\beta(t)]^2 dt
#' }
#'
#' using trapezoidal numerical integration.
#'
#' This quantity is an ISE for a single Monte Carlo replication.
#' Averaging it over replications gives the MISE.
#'
#' @param beta_hat Estimated functional coefficient evaluated on gp.
#' @param beta_true True functional coefficient evaluated on gp.
#' @param gp Functional grid.
#'
#' @return Integrated squared error.
#'
#' @export

beta_ise <- function(
    beta_hat,
    beta_true,
    gp
) {
  
  beta_hat <- as.numeric(beta_hat)
  beta_true <- as.numeric(beta_true)
  gp <- as.numeric(gp)
  
  if (length(beta_hat) != length(beta_true) ||
      length(beta_hat) != length(gp)) {
    stop(
      paste0(
        "'beta_hat', 'beta_true', and 'gp' ",
        "must have the same length."
      ),
      call. = FALSE
    )
  }
  
  if (length(gp) < 2L) {
    stop(
      "'gp' must contain at least two grid points.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(beta_hat)) ||
      any(!is.finite(beta_true)) ||
      any(!is.finite(gp))) {
    stop(
      "Inputs must contain only finite values.",
      call. = FALSE
    )
  }
  
  if (is.unsorted(gp, strictly = TRUE)) {
    stop(
      "'gp' must be strictly increasing.",
      call. = FALSE
    )
  }
  
  squared_error <- (
    beta_hat -
      beta_true
  )^2
  
  h <- diff(gp)
  
  sum(
    h *
      (
        squared_error[-length(squared_error)] +
          squared_error[-1L]
      ) / 2
  )
}


# ==============================================================
# RMSE
# ==============================================================

#' Root mean squared error
#'
#' @param observed Observed values.
#' @param predicted Predicted values.
#'
#' @return Root mean squared error.
#'
#' @export

rmse <- function(
    observed,
    predicted
) {
  
  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  
  .check_metric_vectors(
    observed,
    predicted
  )
  
  sqrt(
    mean(
      (
        observed -
          predicted
      )^2
    )
  )
}


# ==============================================================
# MAE
# ==============================================================

#' Mean absolute error
#'
#' @param observed Observed values.
#' @param predicted Predicted values.
#'
#' @return Mean absolute error.
#'
#' @export

mae <- function(
    observed,
    predicted
) {
  
  observed <- as.numeric(observed)
  predicted <- as.numeric(predicted)
  
  .check_metric_vectors(
    observed,
    predicted
  )
  
  mean(
    abs(
      observed -
        predicted
    )
  )
}


# ==============================================================
# BRIER SCORE
# ==============================================================

#' Brier score
#'
#' Computes the mean squared difference between observed binary
#' outcomes and predicted probabilities.
#'
#' Lower values indicate better probabilistic prediction.
#'
#' @param observed Binary response containing 0 and 1.
#' @param probability Predicted probabilities.
#'
#' @return Brier score.
#'
#' @export

brier_score <- function(
    observed,
    probability
) {
  
  observed <- as.numeric(observed)
  probability <- as.numeric(probability)
  
  .check_metric_vectors(
    observed,
    probability
  )
  
  if (!all(observed %in% c(0, 1))) {
    stop(
      "'observed' must contain only 0 and 1.",
      call. = FALSE
    )
  }
  
  if (any(probability < 0 |
          probability > 1)) {
    stop(
      "'probability' must contain values in [0, 1].",
      call. = FALSE
    )
  }
  
  mean(
    (
      observed -
        probability
    )^2
  )
}


# ==============================================================
# LOG LOSS
# ==============================================================

#' Binary log loss
#'
#' Computes the negative Bernoulli log-likelihood per observation.
#'
#' Lower values indicate better probabilistic prediction.
#'
#' @param observed Binary response containing 0 and 1.
#' @param probability Predicted probabilities.
#' @param eps Numerical clipping constant.
#'
#' @return Binary log loss.
#'
#' @export

log_loss <- function(
    observed,
    probability,
    eps = 1e-15
) {
  
  observed <- as.numeric(observed)
  probability <- as.numeric(probability)
  
  .check_metric_vectors(
    observed,
    probability
  )
  
  if (!all(observed %in% c(0, 1))) {
    stop(
      "'observed' must contain only 0 and 1.",
      call. = FALSE
    )
  }
  
  if (any(probability < 0 |
          probability > 1)) {
    stop(
      "'probability' must contain values in [0, 1].",
      call. = FALSE
    )
  }
  
  if (!is.numeric(eps) ||
      length(eps) != 1L ||
      !is.finite(eps) ||
      eps <= 0 ||
      eps >= 0.5) {
    stop(
      "'eps' must be a finite scalar in (0, 0.5).",
      call. = FALSE
    )
  }
  
  probability <- pmin(
    pmax(
      probability,
      eps
    ),
    1 - eps
  )
  
  -mean(
    observed *
      log(probability) +
      (1 - observed) *
      log(1 - probability)
  )
}


# ==============================================================
# CCR
# ==============================================================

#' Correct classification rate
#'
#' Converts predicted probabilities to binary predictions using
#' a fixed classification cutoff and computes the proportion of
#' correctly classified observations.
#'
#' Higher values indicate better classification performance.
#'
#' @param observed Binary response containing 0 and 1.
#' @param probability Predicted probabilities.
#' @param cutoff Classification threshold. Default is 0.5.
#'
#' @return Correct classification rate in [0, 1].
#'
#' @export

ccr <- function(
    observed,
    probability,
    cutoff = 0.5
) {
  
  observed <- as.numeric(observed)
  probability <- as.numeric(probability)
  
  .check_metric_vectors(
    observed,
    probability
  )
  
  if (!all(observed %in% c(0, 1))) {
    stop(
      "'observed' must contain only 0 and 1.",
      call. = FALSE
    )
  }
  
  if (any(probability < 0 |
          probability > 1)) {
    stop(
      "'probability' must contain values in [0, 1].",
      call. = FALSE
    )
  }
  
  if (!is.numeric(cutoff) ||
      length(cutoff) != 1L ||
      !is.finite(cutoff) ||
      cutoff <= 0 ||
      cutoff >= 1) {
    stop(
      "'cutoff' must be a finite scalar in (0, 1).",
      call. = FALSE
    )
  }
  
  predicted_class <- as.integer(
    probability >= cutoff
  )
  
  mean(
    predicted_class ==
      observed
  )
}


# ==============================================================
# ROC AUC
# ==============================================================

#' Area under the ROC curve
#'
#' Computes the binary ROC AUC from predicted probabilities using
#' the rank-sum representation of the Mann-Whitney statistic.
#'
#' Ties receive half credit through average ranks.
#'
#' Higher values indicate better discrimination.
#'
#' @param observed Binary response containing 0 and 1.
#' @param probability Predicted probabilities.
#'
#' @return ROC AUC in [0, 1]. Returns NA_real_ when only one
#'   response class is present.
#'
#' @export

roc_auc <- function(
    observed,
    probability
) {
  
  observed <- as.numeric(observed)
  probability <- as.numeric(probability)
  
  .check_metric_vectors(
    observed,
    probability
  )
  
  if (!all(observed %in% c(0, 1))) {
    stop(
      "'observed' must contain only 0 and 1.",
      call. = FALSE
    )
  }
  
  if (any(probability < 0 |
          probability > 1)) {
    stop(
      "'probability' must contain values in [0, 1].",
      call. = FALSE
    )
  }
  
  n_positive <- sum(
    observed == 1
  )
  
  n_negative <- sum(
    observed == 0
  )
  
  # ROC AUC is undefined when the test sample contains
  # observations from only one class.
  if (n_positive == 0L ||
      n_negative == 0L) {
    return(NA_real_)
  }
  
  probability_rank <- rank(
    probability,
    ties.method = "average"
  )
  
  rank_sum_positive <- sum(
    probability_rank[
      observed == 1
    ]
  )
  
  auc <- (
    rank_sum_positive -
      n_positive *
      (n_positive + 1) / 2
  ) / (
    n_positive *
      n_negative
  )
  
  as.numeric(auc)
}