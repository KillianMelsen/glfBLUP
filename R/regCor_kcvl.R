#' regCor_kcvl
#'
#' \code{.regCor_kcvl} provides the function that should be optimized to find the optimal regularization penalty.
#'
#' @param penalty Defaults to \code{NULL} but can be specified to override the data-driven regularization.
#' @param data A dataframe containing the genotypes in the first column and the secondary features in the other columns.
#' @param folds Either a numeric value specifying the number of folds that should be used or a list of vectors each containing the genotypes for each fold.
#' @param genoMeans Dataframe with genotypic means.
#' @param reps Number of replicates.
#' @param dopar Boolean specifying whether to use parallelization.
#' @param what Type of correlation matrix that will be regularized.
#'
#' @return The cross-validated negative log-likelihood to be minimized.
#' @keywords internal
#' @importFrom foreach %dopar%
#'
regCor_kcvl <- function(penalty, data, folds, genoMeans, reps, dopar, what) {
  if (!(dopar)) {
    cvLL <- 0
    for (f in 1:length(folds)){
      # Make dataframe for 4/5 folds:
      R_dataframe <- data[!(data[, 1] %in% folds[[f]]), , drop = FALSE]
      R_dataframe$G <- factor(as.character(R_dataframe$G))
      secondaries <- names(R_dataframe)[2:length(names(R_dataframe))]
      R_dataframe[secondaries] <- lapply(R_dataframe[secondaries], as.numeric)

      # Calculating R:
      covmats <- covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps, verbose = FALSE)
      if (what == "phenotypic") {
        R <- stats::cov2cor(covmats$Sp)
      } else if (what == "genetic") {
        R <- stats::cov2cor(covmats$Sg)
      } else if (what == "residual") {
        R <- stats::cov2cor(covmats$Se)
      }

      # Making dataframe for the left out fold:
      S_dataframe <- data[data[, 1] %in% folds[[f]], , drop = FALSE]
      S_dataframe$G <- factor(as.character(S_dataframe$G))
      secondaries <- names(S_dataframe)[2:length(names(S_dataframe))]
      S_dataframe[secondaries] <- lapply(S_dataframe[secondaries], as.numeric)

      # Calculating S:
      covmats <- covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps, verbose = FALSE)
      if (what == "phenotypic") {
        S <- stats::cov2cor(covmats$Sp)
      } else if (what == "genetic") {
        S <- stats::cov2cor(covmats$Sg)
      } else if (what == "residual") {
        S <- stats::cov2cor(covmats$Se)
      }

      if (what == "genetic") {
        nf <- length(folds[[f]])
      } else {
        nf <- dim(data[data[,1] %in% folds[[f]], , drop = FALSE])[1]
      }

      cvLL <- cvLL + nf * regCor_ll(S, regCor_corlw(R, penalty))
    }
    return(cvLL/length(folds))

  } else if (dopar) {
    ### PARALLELIZATION
    cvLL <- foreach::foreach(f = 1:length(folds), .combine = "+", .packages = c("glfBLUP")) %dopar% {
      # Make dataframe for 4/5 folds:
      R_dataframe <- data[!(data[,1] %in% folds[[f]]), , drop = FALSE]
      R_dataframe$G <- factor(as.character(R_dataframe$G))
      secondaries <- names(R_dataframe)[2:length(names(R_dataframe))]
      R_dataframe[secondaries] <- lapply(R_dataframe[secondaries], as.numeric)

      # Calculating R:
      covmats <- glfBLUP::covSS(data = R_dataframe, genoMeans = genoMeans, reps = reps, verbose = FALSE)
      if (what == "phenotypic") {
        R <- stats::cov2cor(covmats$Sp)
      } else if (what == "genetic") {
        R <- stats::cov2cor(covmats$Sg)
      } else if (what == "residual") {
        R <- stats::cov2cor(covmats$Se)
      }

      # Making dataframe for the left out fold:
      S_dataframe <- data[data[,1] %in% folds[[f]], , drop = FALSE]
      S_dataframe$G <- factor(as.character(S_dataframe$G))
      secondaries <- names(S_dataframe)[2:length(names(S_dataframe))]
      S_dataframe[secondaries] <- lapply(S_dataframe[secondaries], as.numeric)

      # Calculating S:
      covmats <- glfBLUP::covSS(data = S_dataframe, genoMeans = genoMeans, reps = reps, verbose = FALSE)
      if (what == "phenotypic") {
        S <- stats::cov2cor(covmats$Sp)
      } else if (what == "genetic") {
        S <- stats::cov2cor(covmats$Sg)
      } else if (what == "residual") {
        S <- stats::cov2cor(covmats$Se)
      }

      if (what == "genetic") {
        nf <- length(folds[[f]])
      } else {
        nf <- dim(data[data[,1] %in% folds[[f]], , drop = FALSE])[1]
      }

      LL <- nf * regCor_ll(S, regCor_corlw(R, penalty))
      return(LL)
    }
    return(cvLL/length(folds))
  }
}
