#' RF
#'
#' \code{RF} is a wrapper around \code{FMradio::RF} that redundancy filters a correlation matrix
#' given a specified correlation threshold.
#'
#' @param cormat A correlation matrix.
#' @param tau Threshold value. The filtered matrix will have no pair of features with an absolute
#' correlation larger than this threshold. Defaults to \code{0.95}.
#' @param verbose Boolean.
#'
#' @return A redundancy filtered correlation matrix.
#' @export
#'
RF <- function(cormat, tau = 0.95, verbose = TRUE) {

  # # Checks:
  # stopifnot("cormat contains non-numeric values!" = all(is.numeric(data[, -1])),
  #           "cormat is not symmetric!" = isSymmetric(cormat),
  #           "cormat has diagonal entries other than 1!" = all(diag(cormat) == 1))

  cormat.filtered <- FMradio::RF(cormat, t = tau)

  if (verbose) {
    cat(sprintf("%d out of %d features remaining...\n", nrow(cormat.filtered), nrow(cormat)))
  }

  return(cormat.filtered)
}
