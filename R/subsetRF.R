#' subsetRF
#'
#' \code{subsetRF} is a wrapper around \code{FMradio::subSet} that subsets a secondary feature dataframe using
#' a redundancy filtered correlation matrix.
#'
#' @param sec A dataframe containing the secondary features to be filtered. The first column contains the genotypes.
#' @param filter The redundancy filtered correlation matrix containing the features that should be kept.
#'
#' @return A redundancy filtered secondary feature dataframe.
#' @export
#'
subsetRF <- function(sec, filter) {
  # stopifnot("filter contains non-numeric values!" = all(is.numeric(filter)),
  #           "filter is not symmetric!" = isSymmetric(filter),
  #           "filter has diagonal entries other than 1!" = all(diag(filter) == 1))

  sec.subset <- as.data.frame(FMradio::subSet(as.matrix(sec[, 2:ncol(sec)]), filter))
  sec.subset <- cbind(sec$G, sec.subset)
  names(sec.subset)[1] <- "G"

  return(sec.subset)
}
