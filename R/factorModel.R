#' factorModel
#'
#' \code{factorModel} fits the ML factor model (FAD) using a correlation matrix and a specified latent dimension.
#'
#' @param data The datamatrix used to estimate the (regularized) correlation matrix used to fit the factor model. Only used to determine the number of
#' individuals (in case of a phenotypic or residual correlation matrix) or number of genotypes (in case of a genetic correlation matrix).
#' @param cormat The correlation matrix.
#' @param m Latent dimension.
#' @param lower Lower bound for the uniquenesses.
#' @param what What type of correlation matrix will be used to fit the factor model.
#' @param verbose Boolean
#'
#' @return A list with the varimax rotated loadings, uniquenesses, rotation matrix, and latent dimension.
#' @export
#'
factorModel <- function(data, cormat, m = NULL, lower = 0.1, what = "genetic", verbose = TRUE) {

  # Checks:
  stopifnot("what must be either 'genetic' or 'residual'!" = what %in% c("genetic", "residual"))

  # MP-bound:
  if (what == "residual") {
    n.obs <- nrow(data)
    ev.thr <- (1 + sqrt(ncol(cormat) / n.obs))^2
  } else if (what == "genetic") {
    n.obs <- length(unique(data$G))
    ev.thr <- (1 + sqrt(ncol(cormat) / n.obs))^2
  }

  if (verbose) {cat(sprintf("M-P bound is %.3f...\n", round(ev.thr, 3)))}

  # Determining latent dimension:
  if (is.null(m)) {
    m <- sum(eigen(cormat)$values > ev.thr)
  }
  if (verbose) {cat(sprintf("Fitting factor model with %d factors...\n", m))}

  fit <- fad::fad(factors = m, covmat = cormat, rotation = "varimax", lower = lower, n.obs = n.obs)

  colnames(fit$loadings) <- paste0("F", 1:m)
  rownames(fit$loadings) <- names(fit$uniquenesses) <- rownames(cormat)
  rownames(fit$rotmat) <- colnames(fit$rotmat) <- colnames(fit$loadings)

  return(list(loadings = fit$loadings,
              uniquenesses = diag(fit$uniquenesses),
              rotmat = fit$rotmat,
              m = m))
}


