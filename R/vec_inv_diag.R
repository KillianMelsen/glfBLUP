#' vec.inv.diag
#'
#' @param x Something
#' @param y Something
#'
#' @return Something
#'
#' @keywords internal
vec.inv.diag <- function(x,y) {

  p <- length(x)
  n <- length(y)
  z <- sapply(x,function(x_i){return(1 / (rep(1,n) + x_i * y))})
  return(z)
}
