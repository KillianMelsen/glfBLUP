#' train2test
#'
#' A function to convert training set BLUPs to test set BLUPs using the kinship.
#'
#' @param BLUPs A vector containing the training set focal trait BLUPs.
#' @param K The kinship matrix.
#' @param train.set A character vector containing the names of the training set genotypes.
#' @param test.set A character vector containing the names of the test set genotypes.
#'
#' @return A vector containing the test set BLUPs.
#'
#' @keywords internal
#'
train2test <- function(BLUPs, K, train.set, test.set) {

  Koo <- K[train.set, train.set]
  Kto <- K[test.set, train.set]
  focalBLUP_test <- Kto %*% solve(Koo) %*% BLUPs
  if (is.null(dim(BLUPs))) {
    focalBLUP_test <- as.numeric(focalBLUP_test)
    names(focalBLUP_test) <- test.set
  } else {
    rownames(focalBLUP_test) <- test.set
  }

  return(focalBLUP_test)
}
