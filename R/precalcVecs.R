#' precalcVecs
#'
#' precalcVecs is an internal function that pre-calculates the vectors (possibly in parallel) used
#' in the efficient regularization of the phenotypic/genetic/residual correlation matrix.
#'
#' @param data A dataframe containing the genotypes in the first column and secondary features
#' in the other columns.
#' @param genoMeans A dataframe containing the genotypic means, with the genotypes in the first column.
#' @param reps Number of replicates
#' @param what Pre-calculation of vectors for the genetic or residual correlation matrix?
#' @param folds List of vectors with the genotypes in the folds
#' @param n.folds Number of folds
#' @param verbose Verbosity
#' @param dopar Parallel calculation or not
#'
#' @return A list of lists containing the required vectors for the regularization
#' @keywords internal
#'
precalcVecs <- function(data, genoMeans, reps, what, folds, n.folds, verbose, dopar) {

  if (verbose) {cat("Pre-calculating vectors for each fold...\n")}

  # Parallel pre-calculation of vectors
  if (dopar) {

    cores <- min(n.folds, parallel::detectCores() - 2)
    doParallel::registerDoParallel(cores = cores)
    cat(sprintf("Parallelizing the pre-calculation process using %d cores...\n", cores))

    vecs <- foreach::foreach(i = 1:n.folds, .packages = c("gfBLUP")) %dopar% {

      # Making a dataframe for 4/5 folds:
      R_dataframe <- droplevels(data[!(data$G %in% folds[[i]]), , drop = FALSE])

      # Calculating eigen decomp of R:
      R.ED <- ifelse(what == "genetic",
                     eigen(stats::cov2cor(covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps)$Sg)),
                     eigen(stats::cov2cor(covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps)$Se)))

      # Making a dataframe for the left out fold:
      S_dataframe <- droplevels(data[data$G %in% folds[[i]], , drop = FALSE])

      # Calculating A:
      S <- ifelse(what == "genetic",
                  stats::cov2cor(covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps)$Sg),
                  stats::cov2cor(covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps)$Se))

      A.diag <-  diag(t(R.ED$vectors) %*% S %*% R.ED$vectors)

      return(list(R.EV = R.ED$values, A.diag = A.diag, nf = dim(S_dataframe)[1]))

    }
    doParallel::stopImplicitCluster()
    return(vecs)

  } else if (!dopar) {

    vecs <- vector("list", n.folds)
    for (i in 1:n.folds) {

      # Making a dataframe for 4/5 folds:
      R_dataframe <- droplevels(data[!(data$G %in% folds[[i]]), , drop = FALSE])

      # Calculating eigen decomp of R:
      R.ED <- ifelse(what == "genetic",
                     eigen(stats::cov2cor(covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps)$Sg)),
                     eigen(stats::cov2cor(covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps)$Se)))

      # Making a dataframe for the left out fold:
      S_dataframe <- droplevels(data[data$G %in% folds[[i]], , drop = FALSE])

      # Calculating A:
      S <- ifelse(what == "genetic",
                  stats::cov2cor(covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps)$Sg),
                  stats::cov2cor(covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps)$Se))

      A.diag <-  diag(t(R.ED$vectors) %*% S %*% R.ED$vectors)

      # Adding both to the list of pre-calculated matrices:
      vecs[[i]] <- list(R.EV = R.ED$values, A.diag = A.diag, nf = dim(S_dataframe)[1])
    }
    return(vecs)
  }
}


