#' createFolds
#'
#' \code{createFolds} creates folds out of a list of genotypes.
#'
#' @param genos Character vector containing the genotypes.
#' @param folds Number of folds to be created
#'
#' @return A list of vectors with the genotypes for each fold.
#' @export
#'
createFolds <- function(genos, folds = 5) {

  # Checks:
  stopifnot("Number of folds must be at least 2!" = folds > 1,
            "Number of folds cannot exceed the number of genotypes!" = folds <= length(genos))

  # Creating folds:
  n.genos <- length(genos)
  fold <- max(min(ceiling(folds), n.genos), 2)
  fold <- rep(1:fold, ceiling(n.genos/fold))[1:n.genos]
  shuffle <- sample(unique(genos), n.genos)
  folds <- split(shuffle, as.factor(fold))

  return(folds)
}
