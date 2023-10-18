#' simMVData
#'
#' \code{simMVData} simulates multi-trait data using given row and column covariance matrices.
#'
#' @param rowCov Row covariance matrix (e.g. kinship or identity).
#' @param colCov Column covariance matrix (e.g. genetic covariance matrix).
#'
#' @return A data matrix.
#' @keywords internal
#'
simMVData <- function(rowCov, colCov) {
  n <- nrow(rowCov)
  p <- ncol(colCov)
  U <- MatrixRoot(rowCov)
  V <- MatrixRoot(colCov)
  W <- matrix(stats::rnorm(n = n * p), ncol = p)
  D  <- U %*% W %*% V
  return(D)
}
