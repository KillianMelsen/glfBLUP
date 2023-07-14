#' redundancyFilter
#'
#' \code{redundancyFilter} is a wrapper around \code{FMradio::RF} that redundancy filters a correlation matrix
#' given a specified correlation threshold.
#'
#' @param data Data matrix to be filtered based on genetic correlations. First column contains genotype identifiers.
#' @param tau Threshold value. The filtered matrix will have no pair of features with an absolute
#' correlation larger than this threshold. Defaults to \code{0.95}.
#' @param verbose Logical enabling verbose mode.
#'
#' @return A list containing the filtered correlation matrix and filtered secondary feature matrix.
#' @export
#'
redundancyFilter <- function(data, tau = 0.95, verbose = TRUE) {

  # # Checks:
  # stopifnot("cormat contains non-numeric values!" = all(is.numeric(data[, -1])),
  #           "cormat is not symmetric!" = isSymmetric(cormat),
  #           "cormat has diagonal entries other than 1!" = all(diag(cormat) == 1))

  Rg <- stats::cov2cor(covSS(data)$Sg)
  Rg.RF <- FMradio::RF(Rg, t = tau)
  if (verbose) {cat(sprintf("%d out of %d features remaining...\n", ncol(Rg.RF), ncol(Rg)))}
  data.subset <- subsetRF(data, Rg.RF)

  return(list(Rg.RF = Rg.RF,
              data.RF = data.subset))
}
