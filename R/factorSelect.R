#' factorSelect
#'
#' \code{factorSelect} uses the specified procedure to determine which factors should be used in genomic prediction.
#'
#' @param data A dataframe containing the genotypes in the first column, followed by the factors, and finally the focal trait.
#' @param procedure Only \code{"leaps"} for now.
#' @param verbose Boolean
#'
#' @return A vector with the names of the selected factors.
#' @export
#'
factorSelect <- function(data, procedure = "leaps", verbose = TRUE) {

  if (procedure == "leaps") {
    data.train <- droplevels(stats::na.omit(data))
    data.train.means <- genotypeMeans(data.train)

    factors <- names(data.train.means)[2:(ncol(data.train.means) - 1)]

    leaps <- leaps::leaps(x = as.matrix(data.train.means[, 2:(ncol(data.train.means) - 1)]),
                          y = as.matrix(data.train.means[, ncol(data.train.means)]),
                          method = "adjr2",
                          names = factors,
                          int = TRUE)

    select.logical <- leaps$which[which(leaps$adjr2 == max(leaps$adjr2)),]
    select <- factors[select.logical]
  }

  if (verbose) {
    cat(sprintf("Selected factors using %s procedure are %s...\n", procedure, paste(select, collapse = " ")))
  }

  return(select)
}
