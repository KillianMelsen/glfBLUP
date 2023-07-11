#' MatrixRoot
#'
#' Calculates the root of a matrix.
#'
#' @param x A matrix
#'
#' @return The root of the matrix \code{x}.
#' @keywords internal
#'
MatrixRoot <- function(x) {

  #x <- (x + t(x)) / 2
  x.eig <- eigen(as.matrix(x), symmetric=TRUE)
  if (length(x.eig$values) > 1) {
    x.sqrt <- x.eig$vectors %*% diag(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  } else {
    x.sqrt <- x.eig$vectors %*% matrix(sqrt(x.eig$values)) %*% solve(x.eig$vectors)
  }
  return(x.sqrt)
}
