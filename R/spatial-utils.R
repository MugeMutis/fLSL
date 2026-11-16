# ==============================================================
# Spatial weight-matrix utilities
# ==============================================================


#' Normalize a spatial weight matrix
#'
#' Row-standardizes a spatial weight matrix so that each row sums
#' to one. The diagonal is set to zero before normalization.
#'
#' This is an internal utility used by the spatial regression,
#' prediction, dimension-reduction, and simulation routines.
#'
#' @param W A numeric square spatial weight matrix.
#' @param tol Numerical tolerance used to determine whether row sums
#'   are effectively zero.
#'
#' @return A row-standardized spatial weight matrix with zero diagonal.
#'
#' @keywords internal
normalize_weights <- function(
    W,
    tol = sqrt(.Machine$double.eps)
) {
  
  ## ============================================================
  ## 0. Sparse weight matrices (Matrix package) keep their format
  ## ============================================================

  if (inherits(W, "sparseMatrix")) {
    return(.normalize_weights_sparse(W, tol = tol))
  }

  ## ============================================================
  ## 1. Input checks
  ## ============================================================

  if (!is.matrix(W)) {
    W <- as.matrix(W)
  }

  if (!is.numeric(W)) {
    stop(
      "'W' must be a numeric matrix.",
      call. = FALSE
    )
  }

  if (nrow(W) != ncol(W)) {
    stop(
      "'W' must be a square matrix.",
      call. = FALSE
    )
  }
  
  if (nrow(W) < 2L) {
    stop(
      "'W' must have dimension at least 2 x 2.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(W))) {
    stop(
      "'W' must contain only finite values.",
      call. = FALSE
    )
  }
  
  if (any(W < 0)) {
    stop(
      "'W' must contain non-negative weights.",
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
  ## 2. Enforce zero diagonal
  ## ============================================================
  
  diag(W) <- 0
  
  
  ## ============================================================
  ## 3. Check row sums
  ## ============================================================
  
  rs <- rowSums(W)
  
  if (any(!is.finite(rs))) {
    stop(
      "Non-finite row sums encountered in 'W'.",
      call. = FALSE
    )
  }
  
  if (any(rs <= tol)) {
    
    bad_rows <- which(
      rs <= tol
    )
    
    stop(
      paste0(
        "'W' contains row(s) with zero or numerically zero ",
        "off-diagonal weight: ",
        paste(bad_rows, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 4. Row standardization
  ## ============================================================
  
  W_normalized <- W / rs
  
  
  ## ============================================================
  ## 5. Numerical cleanup
  ## ============================================================
  
  diag(W_normalized) <- 0
  
  
  ## ============================================================
  ## 6. Integrity check
  ## ============================================================
  
  normalized_row_sums <- rowSums(
    W_normalized
  )
  
  if (max(abs(normalized_row_sums - 1)) > 1e-8) {
    
    stop(
      "Row normalization of 'W' failed numerically.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 7. Return
  ## ============================================================
  
  W_normalized
}



#' Row-standardize a sparse weight matrix
#'
#' Sparse counterpart of normalize_weights(): sets the diagonal to zero,
#' checks nonnegativity and positive row sums, and row-standardizes
#' without densifying the matrix.
#'
#' @keywords internal
.normalize_weights_sparse <- function(W, tol = sqrt(.Machine$double.eps)) {
  if (!requireNamespace("Matrix", quietly = TRUE))
    stop("Package 'Matrix' is required for sparse weight matrices.", call. = FALSE)
  n <- nrow(W)
  if (n != ncol(W) || n < 2L)
    stop("'W' must be a square matrix of dimension at least 2.", call. = FALSE)
  if (inherits(W, "symmetricMatrix") || inherits(W, "triangularMatrix"))
    W <- methods::as(W, "generalMatrix")
  W <- methods::as(W, "CsparseMatrix")
  vals <- W@x
  if (!is.numeric(vals) || any(!is.finite(vals)) || any(vals < 0))
    stop("'W' must contain finite, nonnegative weights.", call. = FALSE)
  W <- W - Matrix::Diagonal(n = n, x = Matrix::diag(W))
  W <- Matrix::drop0(W)
  rs <- Matrix::rowSums(W)
  if (any(!is.finite(rs)) || any(rs <= tol))
    stop("'W' contains row(s) with zero or numerically zero off-diagonal weight: ",
         paste(utils::head(which(rs <= tol), 10L), collapse = ", "), ".",
         call. = FALSE)
  W <- Matrix::Diagonal(n = n, x = 1 / rs) %*% W
  if (max(abs(Matrix::rowSums(W) - 1)) > 1e-8)
    stop("Row normalization of sparse 'W' failed numerically.", call. = FALSE)
  W
}



#' Construct a spatial weight matrix
#'
#' Constructs a spatial weight matrix from two-dimensional coordinates
#' using one of three weighting schemes:
#'
#' \itemize{
#'   \item \code{"knn"}: k-nearest-neighbor weights;
#'   \item \code{"idw"}: inverse-distance weights;
#'   \item \code{"exp"}: exponentially decaying distance weights.
#' }
#'
#' The diagonal is set to zero and the resulting matrix is
#' row-standardized.
#'
#' For \code{"knn"}, each observation is connected to its \code{k}
#' nearest neighbors with equal raw weights. The resulting matrix is
#' generally asymmetric.
#'
#' For \code{"idw"}, the raw off-diagonal weights are
#'
#' \deqn{
#'   w_{ij}^{(0)} = d_{ij}^{-\mathrm{power}}.
#' }
#'
#' For \code{"exp"}, the raw off-diagonal weights are
#'
#' \deqn{
#'   w_{ij}^{(0)} = \exp(-\mathrm{decay}\,d_{ij}).
#' }
#'
#' @param coords A numeric matrix or data frame with observations in
#'   rows and two coordinate columns.
#' @param type Spatial weighting scheme. One of \code{"knn"},
#'   \code{"idw"}, or \code{"exp"}.
#' @param k Number of nearest neighbors when \code{type = "knn"}.
#' @param power Positive power parameter for inverse-distance weights.
#' @param decay Positive decay parameter for exponential weights.
#'
#' @return A row-standardized spatial weight matrix.
#'
#' @examples
#' set.seed(123)
#' coords <- cbind(
#'   x = runif(50),
#'   y = runif(50)
#' )
#'
#' W <- generate_W(
#'   coords = coords,
#'   type = "knn",
#'   k = 4
#' )
#'
#' @export
generate_W <- function(
    coords,
    type = c("knn", "idw", "exp"),
    k = 4,
    power = 2,
    decay = 20
) {
  
  ## ============================================================
  ## 1. Match weighting scheme
  ## ============================================================
  
  type <- match.arg(
    type
  )
  
  
  ## ============================================================
  ## 2. Coordinate checks
  ## ============================================================
  
  if (!is.matrix(coords)) {
    coords <- as.matrix(coords)
  }
  
  if (!is.numeric(coords)) {
    stop(
      "'coords' must be numeric.",
      call. = FALSE
    )
  }
  
  if (ncol(coords) != 2L) {
    stop(
      "'coords' must contain exactly two columns.",
      call. = FALSE
    )
  }
  
  if (nrow(coords) < 2L) {
    stop(
      "'coords' must contain at least two observations.",
      call. = FALSE
    )
  }
  
  if (any(!is.finite(coords))) {
    stop(
      "'coords' must contain only finite values.",
      call. = FALSE
    )
  }
  
  
  n <- nrow(
    coords
  )
  
  
  ## ============================================================
  ## 3. Reject duplicate spatial locations
  ##
  ## Duplicate locations create zero off-diagonal distances and
  ## therefore make inverse-distance weighting undefined.
  ## Keeping this check for every weighting scheme also ensures
  ## that the spatial network is constructed from unique locations.
  ## ============================================================
  
  if (anyDuplicated(
    data.frame(
      x = coords[, 1L],
      y = coords[, 2L]
    )
  ) > 0L) {
    
    stop(
      "'coords' contains duplicate spatial locations.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 4. Pairwise Euclidean distance matrix
  ##
  ## The simulation coordinates are generated on [0,1]^2, so the
  ## Euclidean metric is the intended metric in the simulation DGP.
  ## ============================================================
  
  D <- as.matrix(
    stats::dist(coords)
  )
  
  diag(D) <- Inf
  
  
  ## ============================================================
  ## 5. Construct raw weights
  ## ============================================================
  
  if (type == "knn") {
    
    ## ----------------------------------------------------------
    ## 5.1 K-nearest-neighbor weights
    ## ----------------------------------------------------------
    
    if (!is.numeric(k) ||
        length(k) != 1L ||
        !is.finite(k) ||
        k != as.integer(k)) {
      
      stop(
        "'k' must be a single finite integer.",
        call. = FALSE
      )
    }
    
    k <- as.integer(
      k
    )
    
    if (k < 1L || k >= n) {
      
      stop(
        "'k' must satisfy 1 <= k < n.",
        call. = FALSE
      )
    }
    
    
    W_raw <- matrix(
      0,
      nrow = n,
      ncol = n
    )
    
    
    for (i in seq_len(n)) {
      
      neighbors_i <- order(
        D[i, ],
        decreasing = FALSE
      )[seq_len(k)]
      
      W_raw[
        i,
        neighbors_i
      ] <- 1
    }
    
    
  } else if (type == "idw") {
    
    ## ----------------------------------------------------------
    ## 5.2 Inverse-distance weights
    ## ----------------------------------------------------------
    
    if (!is.numeric(power) ||
        length(power) != 1L ||
        !is.finite(power) ||
        power <= 0) {
      
      stop(
        "'power' must be a single positive finite number.",
        call. = FALSE
      )
    }
    
    
    W_raw <- D^(-power)
    
    diag(W_raw) <- 0
    
    
  } else {
    
    ## ----------------------------------------------------------
    ## 5.3 Exponential-distance weights
    ## ----------------------------------------------------------
    
    if (!is.numeric(decay) ||
        length(decay) != 1L ||
        !is.finite(decay) ||
        decay <= 0) {
      
      stop(
        "'decay' must be a single positive finite number.",
        call. = FALSE
      )
    }
    
    
    W_raw <- exp(
      -decay * D
    )
    
    diag(W_raw) <- 0
  }
  
  
  ## ============================================================
  ## 6. Raw-weight checks
  ## ============================================================
  
  if (any(!is.finite(W_raw))) {
    stop(
      "Construction of the spatial weights produced non-finite values.",
      call. = FALSE
    )
  }
  
  if (any(W_raw < 0)) {
    stop(
      "Construction of the spatial weights produced negative values.",
      call. = FALSE
    )
  }
  
  
  ## ============================================================
  ## 7. Row standardization
  ## ============================================================
  
  W <- normalize_weights(
    W_raw
  )
  
  
  ## ============================================================
  ## 8. Return
  ## ============================================================
  
  W
}