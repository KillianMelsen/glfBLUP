#' RF
#'
#' \code{RF} is a wrapper around \code{FMradio::RF} that redundancy filters a correlation matrix
#' given a specified correlation threshold.
#'
#' @param cormat A correlation matrix.
#' @param sec An optional argument for a secondary feature data matrix that allows immediate filtering of the data matrix.
#' @param tau Threshold value. The filtered matrix will have no pair of features with an absolute
#' correlation larger than this threshold. Defaults to \code{0.95}.
#' @param verbose Boolean.
#'
#' @return A redundancy filtered correlation matrix or a list containing the RF correlation matrix and RG secondary feature matrix.
#' @export
#'
RF <- function(cormat, sec = NULL, tau = 0.95, verbose = TRUE) {

  # # Checks:
  # stopifnot("cormat contains non-numeric values!" = all(is.numeric(data[, -1])),
  #           "cormat is not symmetric!" = isSymmetric(cormat),
  #           "cormat has diagonal entries other than 1!" = all(diag(cormat) == 1))

  cormat.filtered <- FMradio::RF(cormat, t = tau)
  if (verbose) {cat(sprintf("%d out of %d features remaining...\n", nrow(cormat.filtered), nrow(cormat)))}

  if (is.null(sec)) {
    return(cormat.filtered)
  } else {
    sec.subset <- subsetRF(sec, cormat.filtered)
    return(list(cormat.RF = cormat.filtered,
                sec.RF = sec.subset))
  }




}
