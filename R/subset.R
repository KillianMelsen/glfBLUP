#' subset
#'
#' \code{subset} is a wrapper around \code{FMradio::subSet} that subsets a dataframe using
#' a redundancy filtered correlation matrix.
#'
#' @param data A dataframe containing the features to be filtered. The first column contains the genotypes, the last
#' column contains the focal trait.
#' @param filter The redundancy filtered correlation matrix containing the features that should be kept.
#'
#' @return A redundancy filtered dataframe.
#' @export
#'
subset <- function(data, filter) {
  # stopifnot("filter contains non-numeric values!" = all(is.numeric(filter)),
  #           "filter is not symmetric!" = isSymmetric(filter),
  #           "filter has diagonal entries other than 1!" = all(diag(filter) == 1))

  data.subset <- as.data.frame(FMradio::subSet(as.matrix(data[, 2:(ncol(data) - 1)]), filter))
  data.subset <- cbind(data$G, data.subset, data$Y)
  names(data.subset)[1] <- "G"
  names(data.subset)[ncol(data.subset)] <- "Y"

  return(data.subset)
}
