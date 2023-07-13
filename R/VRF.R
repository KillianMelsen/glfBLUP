#' VRF
#'
#' \code{VRF} variance and redundancy filters a secondary feature data matrix. First, features with
#' negative estimated genetic and/or residual variances are filtered out. Then, the remaining features
#' are redundancy filtered based on absolute genetic correlations. Note that only training data is used for the filtering.
#'
#' @param sec A secondary feature dataframe where the first column contains the genotypes.
#' @param train.set A character vector with the names of the training set genotypes.
#' @param tau Threshold value. The filtered matrix will have no pair of features with an absolute
#' correlation larger than this threshold. Defaults to \code{0.95}.
#' @param verbose Boolean.
#'
#' @return A list with the VR filtered genetic correlation and covariance matrices, as well as the VR filtered secondary data matrix.
#' @export
#'
VRF <- function(sec, train.set, tau = 0.95, verbose = TRUE) {

  # # Checks:
  # stopifnot("cormat contains non-numeric values!" = all(is.numeric(data[, -1])),
  #           "cormat is not symmetric!" = isSymmetric(cormat),
  #           "cormat has diagonal entries other than 1!" = all(diag(cormat) == 1))

  # Variance filtering (i.e. discard features with negative genetic or residual variances):
  covmats <- covSS(sec[which(sec$G %in% train.set),], use_nearPD = TRUE)
  keep <- diag(covmats$Sg) > 0 & diag(covmats$Se) > 0
  if (verbose) {cat(sprintf("%d out of %d features discarded due to negative variances...\n", sum(!keep), ncol(covmats$Sg)))}
  covmats$Sp <- covmats$Sp[keep, keep]
  covmats$Sg <- covmats$Sg[keep, keep]
  covmats$Se <- covmats$Se[keep, keep]

  # Redundancy filtering:
  cormat.filtered <- FMradio::RF(stats::cov2cor(covmats$Sg), t = tau)
  if (verbose) {cat(sprintf("%d out of %d features remaining...\n", nrow(cormat.filtered), nrow(covmats$Sg)))}

  D <- covmats$Sg[rownames(cormat.filtered), colnames(cormat.filtered)]
  D <- sqrt(diag(D))

  sec.subset <- subsetRF(sec, cormat.filtered)
  return(list(cormat.RF = cormat.filtered,
              covmat.RF = outer(D, D) * cormat.filtered,
              sec.RF = sec.subset))
}
