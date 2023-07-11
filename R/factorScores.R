#' factorScores
#'
#' \code{factorScores} calculates the factor scores for each individual based on the loadings and uniquenesses, as well as the latent dimension.
#'
#' @param data Dataframe containing the secondary features used to estimate the factor scores. Genotypes should be in the first column.
#' @param loadings Loadings matrix.
#' @param uniquenesses Uniquenesses matrix.
#' @param m Latent dimension.
#' @param type Type of factor scores to be calculated.
#'
#' @return A dataframe with the factor scores.
#' @export
#'
factorScores <- function(data, loadings, uniquenesses, m, type = "thomson") {

  # Checks:
  stopifnot("Mismatch between loadings matrix and latent dimension!" = m == dim(loadings)[2])

  datamat <- as.matrix(data[, -1])

  # Filling the dataframe:
  psi.i <- diag(1/diag(uniquenesses))

  if (type == "thomson"){
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% solve(diag(m) + t(loadings) %*% psi.i %*% loadings))
  }
  if (type == "bartlett"){
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% solve(t(loadings) %*% psi.i %*% loadings))
  }
  if (type == "anderson"){
    G <- t(loadings) %*% psi.i %*% loadings
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% expm::sqrtm(solve((G %*% (diag(m) + G)))))
  }

  names(scores) <- paste0("F", 1:m)
  scores <- cbind(data$G, scores)

  return(scores)
}
