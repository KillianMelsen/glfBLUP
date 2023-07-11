#' gBLUPo
#'
#' A fast implementation to calculate the training set gBLUPs.
#'
#' @param Y A matrix containing the mean factor and focal trait values. MUST NOT CONTAIN THE GENOTYPE COLUMN.
#' @param K The kinship matrix.
#' @param X Matrix of covariates.
#' @param B Matrix of coefficients.
#' @param Sg Possibly subsetted genetic covariance matrix of factors and focal trait.
#' @param Se Possibly subsetted residual covariance matrix of factors and focal trait.
#' @param train.set Character vector that contains the names of the training set genotypes.
#'
#' @return A list containing the training set gBLUPs for all factors/traits, or just the focal trait.
#' @keywords internal
#'
gBLUPo <- function(Y, K, X = matrix(rep(1, nrow(K[train.set, train.set]))),
                   B = t(matrix(rep(0, ncol(Y)))), Sg, Se, train.set) {

  p <- ncol(Y)
  Koo <- K[train.set, train.set]
  C <- solve(Sg)
  D <- solve(Se)
  nc <- ncol(X)

  R <- solve(Koo)
  w <- eigen(R)
  Lambda.R <- diag(w$values)
  U <- w$vectors

  D.sqrt.inv <- MatrixRoot(solve(D))
  w <- eigen(D.sqrt.inv %*% C %*% D.sqrt.inv)
  Q.1 <- w$vectors
  Lambda.1 <- w$values

  if (nc > 0) {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% (Y - X %*% B) %*% MatrixRoot(D) %*% Q.1)
  } else {
    S.1     <- vec.inv.diag(x = Lambda.1, y = diag(Lambda.R)) * (t(U) %*% Y %*% MatrixRoot(D) %*% Q.1)
  }

  D.sqrt.inv <- MatrixRoot(solve(D))
  w <- eigen(D.sqrt.inv %*% C %*% D.sqrt.inv)
  Q.1 <- w$vectors

  # All training gBLUPs:
  MU <- matrix(kronecker(D.sqrt.inv %*% Q.1, U) %*% matrix(S.1), ncol = p)
  colnames(MU) <- colnames(Y)
  rownames(MU) <- rownames(Y)
  focal <- MU[, "Y"]
  focal <- as.numeric(focal)
  names(focal) <- rownames(MU)
  return(list(all_gBLUPs = MU,
              focal_gBLUP = focal))
}
