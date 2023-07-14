#' covSS
#'
#' \code{covSS} calculates the phenotypic, genetic, and residual covariance matrix. If bot \code{sec} and \code{foc}
#' are specified, individual replicates must be in the same rows in each dataframe!
#'
#' @param data Dataframe containing the genotypes in the first column and features/traits for which (co)variances should
#' be estimated in the other columns.
#' @param genoMeans Genotypic means. \code{NULL} by default and calculated internally. Can be specified to avoid
#' repeated calculations if covSS is called for the same genotypes multiple times.
#' @param reps Number of replicates. Also \code{NULL} by default but can be specified for the same reasons
#' as \code{genoMeans}.
#' @param sepExp Boolean specifying whether the focal trait was measured on different plants than the secondary
#' features. Only use if there is a focal trait. Sets the residual covariances between secondary features or factors
#' and the focal to \code{0}.
#' @param verbose Boolean
#'
#' @return A list containing the phenotypic, genetic, and residual covariance matrices, as well as a vector
#' containing the plot level heritabilities.
#' @export
#'
covSS <- function(data, genoMeans = NULL, reps = NULL, sepExp = FALSE, verbose = TRUE) {

  # # Checks:
  # stopifnot("NA values in the data!" = all(!is.na(data)),
  #           "Non-numeric values in the data!" = all(is.numeric(data[, -1])))

  # stopifnot("No data supplied!" = is.null(sec) & is.null(foc))
  # if (!is.null(sec) & is.null(foc)) {
  #   stopifnot("genoMeans does not match sec!" = nrow(genoMeans) == length(unique(sec$G)))
  #   data <- sec
  # } else if (is.null(sec) & !is.null(foc)) {
  #   stopifnot("genoMeans does not match foc!" = nrow(genoMeans) == length(unique(foc$G)))
  #   data <- foc
  # } else if (!is.null(sec) & !is.null(foc)) {
  #   stopifnot("genoMeans$sec does not match sec!" = nrow(genoMeans$sec) == length(unique(sec$G)),
  #             "genoMeans$foc does not match foc!" = nrow(genoMeans$foc) == length(unique(foc$G)))
  #   data <- cbind(sec, foc[, -1])
  # }

  # Determine number of traits (number of columns - 1 because of genotype column):
  n.features <- ncol(data) - 1

  # Calculate genotype means if genoMeans is NULL:
  if (is.null(genoMeans)) {genoMeans <- genotypeMeans(data)}

  # Determine number of replicates if n.rep is NULL:
  if (is.null(reps)) {
    n.rep.vector <- as.integer(table(data[, 1]))
    reps <- (sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector)) / (length(n.rep.vector) - 1)
  }

  # Calculate total sums of squares:
  overallMeans <- rep(NA, n.features)
  for (i in 2:(n.features + 1)) {
    overallMeans[i - 1] <- mean(data[, i])
  }
  diffs.T <- as.matrix(data[, -1]) - kronecker(matrix(rep(1, nrow(data))), t(matrix(overallMeans)))
  SS.T <- t(diffs.T) %*% (diffs.T)

  # Calculate residual sums of squares:
  # Create matrix with genotype means that is conformable to data:
  genoMeans.matched <- genoMeans[match(data[, 1], genoMeans[, 1]),]
  rownames(genoMeans.matched) <- rownames(data)

  diffs.E <- as.matrix(data[, -1]) - as.matrix(genoMeans.matched[, -1])
  SS.E <- t(diffs.E) %*% (diffs.E)

  # Calculate genotypic sums of squares:
  SS.G <- SS.T - SS.E

  # Calculate genotypic and residual mean squares:
  MS.G <- SS.G / (length(unique(data[, 1])) - 1)
  MS.E <- SS.E / (nrow(data) - length(unique(data[, 1])))

  # Set residual covariances between secondary features or factors and focal
  # trait to 0 if measured on different plants:
  if (sepExp) {MS.E[nrow(MS.E), 1:(ncol(MS.E) - 1)] <- MS.E[1:(nrow(MS.E) - 1), ncol(MS.E)] <- 0}

  # Calculate genetic covariance matrix:
  Sg <- (MS.G - MS.E) / reps

  # Ensure symmetry:
  Sg <- (Sg + t(Sg)) / 2

  # Enforce PDness if specified:
  if (min(eigen(Sg)$values) < 0) {
    Sg <- as.matrix(Matrix::nearPD(Sg, keepDiag = FALSE)$mat)
    if (verbose) {cat("covSS() using nearPD because min(d(MS.G)) < max(d(MS.E))...\n")}
  }

  # Something else...:
  temp <- which(diag(Sg) < 0)
  if (length(temp) > 0) {
    if (verbose) {cat("Resolving negative variances in Sg...\n")}
    Sg[temp, ] <- 0
    Sg[, temp] <- 0
    if (min(eigen(Sg)$values) < 0) {
      Sg <- as.matrix(Matrix::nearPD(Sg, keepDiag = TRUE)$mat)
    }
  }

  h2s <- diag(Sg)/(diag(Sg + MS.E))
  names(h2s) <- colnames(Sg)

  # Calculate phenotypic covariance matrix and return all matrices including the heritabilities:
  return(list(Sp = Sg + MS.E,
              Sg = Sg,
              Se = MS.E,
              h2s = h2s))
}
