#' regCor
#'
#' \code{regCor} regularizes a phenotypic, genetic, or residual correlation matrix
#' using k-fold cross-validation.
#'
#' @param data A dataframe containing the genotypes in the first column and secondary features
#' in the other columns.
#' @param folds A list of vectors containing the genotypes for each fold, or a number specifying how
#' many folds must be created.
#' @param what Which correlation matrix must be regularized? Can be \code{"phenotypic"}, \code{"genetic"}, or \code{"residual"}.
#' @param dopar Boolean specifying whether to use parallelization or not.
#' @param use_nearPD Boolean whether to use \code{Matrix::nearPD} to enforce PDness or not.
#' @param penalty An optional argument containing a user-specified penalty between \code{0} and \code{1}. Specifying one overrides the data-driven regularization.
#' @param verbose Boolean.
#'
#' @return A list containing the optimal penalty and the regularized correlation matrix.
#' @export
#'
regCor <- function(data, folds = 5, what = "genetic", dopar = FALSE, use_nearPD = TRUE, penalty = NULL, verbose = TRUE) {

  # # Checks:
  # stopifnot("NA values in the data!" = all(!is.na(data)),
  #           "Non-numeric values in the data!" = all(is.numeric(data[, -1])),
  #           "folds is neither a list or a number!" = class(folds) %in% c("numeric", "integer", "list"),
  #           "what must be either 'phenotypic', 'genetic', or 'residual'!" = what %in% c("phenotypic", "genetic", "residual"))

  if (inherits(folds, "list")) {
    stopifnot("Number of folds must be at least 2!" = length(folds) > 1,
              "Number of folds cannot exceed the number of genotypes!" = length(folds) <= length(unique(data$G)))
  } else if (inherits(folds, c("numeric", "integer"))) {
    stopifnot("Number of folds must be at least 2!" = folds > 1,
              "Number of folds cannot exceed the number of genotypes!" = folds <= length(unique(data$G)))
  }

  # Calculate genotype means:
  genoMeans <- genoMeans(data)

  # Determine number of reps:
  n.rep.vector <- as.integer(table(data$G))
  reps <- (sum(n.rep.vector) - sum(n.rep.vector^2)/sum(n.rep.vector))/(length(n.rep.vector) - 1)

  # Data-driven regularization:
  if (is.null(penalty)) {

    # Determining folds if they were not specified in a list:
    if (!inherits(folds, "list")) {
      if (verbose) {cat("Creating folds...\n")}
      folds <- createFolds(genos = unique(data$G), folds = folds)
      }
    if (verbose) {cat("Determining optimal penalty value...\n")}

    # Setting up parallelization if required:
    if (dopar) {
      cores <- min(5, parallel::detectCores() - 2)
      doParallel::registerDoParallel(cores = cores)
      if (verbose) {cat(sprintf("Parallelizing the optimization process using %d cores...\n", cores))}
    }

    # Optimization:
    optPenalty <- stats::optim(0.5, regCor_kcvl, method = "Brent", lower = 0, upper = 1,
                               data = data, folds = folds, genoMeans = genoMeans, reps = reps,
                               dopar = dopar, use_nearPD = use_nearPD, what = what)$par

    if (dopar) {doParallel::stopImplicitCluster()}
  } else {
    # Use specified penalty if required:
    optPenalty <- penalty
  }

  covmats <- covSS(data = data, genoMeans = genoMeans, reps = reps, use_nearPD = use_nearPD)

  if (verbose) {cat(sprintf("Optimal penalty value set at %.3f...\n", round(optPenalty, 3)))}

  if (what == "phenotypic") {
    return(list(optPen = optPenalty,
                optCor = regCor_corlw(stats::cov2cor(covmats$Sp), optPenalty),
                Sp = covmats$Sp))
  } else if (what == "genetic") {
    return(list(optPen = optPenalty,
                optCor = regCor_corlw(stats::cov2cor(covmats$Sg), optPenalty),
                Sg = covmats$Sg))
  } else if (what == "residual") {
    return(list(optPen = optPenalty,
                optCor = regCor_corlw(stats::cov2cor(covmats$Se), optPenalty),
                Se = covmats$Se))
  }
}
