#' factorScores
#'
#' \code{factorScores} calculates the factor scores for each individual based on the loadings and
#' uniquenesses, as well as the latent dimension. IMPORTANT: for \code{"genetic-thomson"} factor scores, both the loadings, uniquenesses,
#' as well as the (regularized) residual covariance matrix must be on the covariance scale.
#'
#' @param data Dataframe containing the secondary features used to estimate the factor scores. Genotypes should be in the first column.
#' @param loadings Loadings matrix.
#' @param uniquenesses Uniquenesses matrix.
#' @param m Latent dimension.
#' @param type Type of factor scores to be calculated.
#' @param Se \code{NULL} by default. Must be specified if \code{"genetic-thomson"} is used as the type.
#'
#' @return A dataframe with the factor scores.
#' @export
#'
factorScores <- function(data, loadings, uniquenesses, m, type = "thomson", Se = NULL) {

  # Checks:
  stopifnot("Mismatch between loadings matrix and latent dimension!" = m == dim(loadings)[2])
  if (type == "genetic-thomson") {
    stopifnot("Se must be specified for 'genetic-thomson' factor scores" = !is.null(Se))
  }

  datamat <- as.matrix(data[, -1])

  # Filling the dataframe:
  if (type == "thomson"){
    psi.i <- diag(1/diag(uniquenesses))
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% solve(diag(m) + t(loadings) %*% psi.i %*% loadings))
  } else if (type == "bartlett"){
    psi.i <- diag(1/diag(uniquenesses))
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% solve(t(loadings) %*% psi.i %*% loadings))
  } else if (type == "anderson"){
    psi.i <- diag(1/diag(uniquenesses))
    G <- t(loadings) %*% psi.i %*% loadings
    scores <- as.data.frame(datamat %*% psi.i %*% loadings %*% expm::sqrtm(solve((G %*% (diag(m) + G)))))
  } else if (type == "genetic-thomson") {
    A.i <- solve(uniquenesses + Se)
    scores <- as.data.frame(datamat %*% A.i %*% loadings %*% solve(diag(m) + t(loadings) %*% A.i %*% loadings))
  }

  names(scores) <- paste0("F", 1:m)
  scores <- cbind(data$G, scores)

  return(scores)
}
