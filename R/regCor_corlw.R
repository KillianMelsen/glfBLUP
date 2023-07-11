#' regCor_corlw
#'
#' \code{.regCor_corlw} calculates the penalized correlation matrix using an identity
#' matrix target.
#'
#' @param R A correlation matrix.
#' @param penalty The regularization penalty.
#'
#' @return A penalized correlation matrix.
#' @keywords internal
#'
regCor_corlw <- function(R, penalty) {
      return((1 - penalty) * R + penalty * diag(dim(R)[1]))
}
