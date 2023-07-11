#' genoMeans
#'
#' \code{.genoMeans} calculates the per-genotype means for secondary or focal trait data.
#'
#' @param data Dataframe containing the genotypes in the first column.
#'
#' @return A dataframe with genotypic means.
#' @keywords internal
#'
genoMeans <- function(data) {

  # Calculate genotypic means:
  data.means <- stats::aggregate(data[, 2:ncol(data)], list(data[, 1]), mean)

  # Specify correct column names and make sure the genotype column is character:
  names(data.means)[1] <- names(data)[1]
  data.means$G <- as.character(data.means$G)

  return(data.means)
}
