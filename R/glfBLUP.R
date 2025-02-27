#' glfBLUP
#'
#' \code{glfBLUP} uses genetic factors to produce multi-trait genomic predictions.
#'
#' @param data The dataframe containing the secondary features (factors) and the focal trait. The column
#' order should be G-F1-...-Fm-Y.
#' @param selection Optional argument. Output from \code{glfBLUP::factorSelect}. All factors are used if left unspecified.
#' @param K The kinship matrix with the genotype identifiers as row and column names.
#' @param sepExp Boolean indicating whether the focal trait and secondary features were measured on different plants/plots.
#' @param returnBLUPs Boolean indicating whether a matrix with BLUPs should be returned.
#' @param verbose Boolean
#'
#' @return A list with the test set predictions, factor selection, and the genetic and residual covariance matrices of factors
#' and focal trait, as well as their plot-level heritabilities.
#' @export
#'
glfBLUP <- function(data, selection = NULL, K, sepExp = FALSE, returnBLUPs = FALSE, verbose = TRUE) {

  # Determining training and test set as well as the scenario (CV1 or CV2):
  train.set <- as.character(unique(data$G[!is.na(data$Y)]))
  test.set <- as.character(setdiff(unique(data$G), train.set))

  if (all(is.na(data[which(data$G %in% test.set), 2:(ncol(data) - 1)]))) {
    CV <- "CV1"
  } else if (!any(is.na(data[which(data$G %in% test.set), 2:(ncol(data) - 1)]))) {
    CV <- "CV2"
  } else {
    stop("Either all secondary traits must be present for the test set, or missing...")
  }

  # Number of replicates:
  n.rep.vector <- as.integer(table(data$G))
  reps <- (sum(n.rep.vector) - sum(n.rep.vector^2) / sum(n.rep.vector)) / (length(n.rep.vector) - 1)

  # Subsetting data to the selected factors if specified:
  if (!is.null(selection)) {
    data.subset <- data[, c("G", selection, "Y")]
  } else {
    selection <- names(data)[2:(length(names(data)) - 1)]
    data.subset <- data[, c("G", selection, "Y")]
  }

  # Determining covariances between factors and focal trait:
  if (verbose) {cat("Determining covariances between factors and focal trait...\n")}
  covmats <- covSS(data.subset[which(data.subset$G %in% train.set),], sepExp = sepExp, verbose = TRUE)

  # Reducing pseudoCRD data to genotypic means:
  data.subset.means <- genotypeMeans(data.subset)
  data.subset.means <- data.subset.means[match(c(test.set, train.set), data.subset.means$G),]
  rownames(data.subset.means) <- c(test.set, train.set)

  # Calculating training set BLUPs:
  BLUPs_train <- gBLUPo(Y = as.matrix(data.subset.means[train.set, -1]),
                        K = K,
                        Sg = covmats$Sg,
                        Se = covmats$Se / reps,
                        train.set = train.set)

  if (CV == "CV1") {
    # Calculating focal trait test set BLUPs:
    focalBLUP_test <- train2test(BLUPs_train$focal_gBLUP, K = K, train.set = train.set, test.set = test.set)
  } else if (CV == "CV2") {
    # Calculating all test set BLUPs:
    allBLUPs_test <- train2test(BLUPs_train$all_gBLUPs, K = K, train.set = train.set, test.set = test.set)

    BLUP_Y_test <- allBLUPs_test[, "Y"]
    BLUP_S_test <- allBLUPs_test[, setdiff(colnames(allBLUPs_test), "Y")]
    Ktt <- K[test.set, test.set]
    Kto <- K[test.set, train.set]
    Koo <- K[train.set, train.set]
    Vc <- kronecker(covmats$Sg[selection, selection], solve(Ktt)) + kronecker(covmats$Se[selection, selection] / reps, diag(dim(Ktt)[1]))

    # Final CV2 test set focal trait BLUP calculation following two-step approach (Runcie & Cheng, 2019):
    focalBLUP_test <- as.numeric(matrix(BLUP_Y_test) + kronecker(t(as.matrix(covmats$Sg[selection, "Y"])), solve(Ktt)) %*%
                                   solve(Vc) %*% (matrix(as.matrix(data.subset.means[test.set, selection])) - matrix(BLUP_S_test)))

    names(focalBLUP_test) <- test.set
  }

  if (returnBLUPs & CV == "CV2") {
    BLUPs <- rbind(BLUPs_train$all_gBLUPs, allBLUPs_test)
  } else if (returnBLUPs & CV == "CV1") {
    allBLUPs_test <- train2test(BLUPs_train$all_gBLUPs, K = K, train.set = train.set, test.set = test.set)
    BLUPs <- rbind(BLUPs_train$all_gBLUPs, allBLUPs_test)
  }

  h2s <- diag(covmats$Sg) / diag(covmats$Sg + covmats$Se)
  names(h2s) <- colnames(covmats$Sg)

  if (returnBLUPs) {
    ret <- list(preds = focalBLUP_test,
                factorSelection = selection,
                Sg = covmats$Sg,
                Se = covmats$Se,
                h2s = h2s,
                BLUPs = BLUPs)
  } else {
    ret <- list(preds = focalBLUP_test,
                factorSelection = selection,
                Sg = covmats$Sg,
                Se = covmats$Se,
                h2s = h2s)
  }
  return(ret)
}
