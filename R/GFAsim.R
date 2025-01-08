#' GFAsim
#'
#' \code{GFAsim} simulates genetic latent factor driven data, including the actual latent factor scores for benchmarking purposes.
#'
#' @param K Kinship matrix.
#' @param r Number of replicates per genotype.
#' @param n.LSF Number of latent signal factors (i.e. factors that the focal trait loads on).
#' @param n.LNF Number of latent noise factors.
#' @param LSF.rho Between latent signal factor correlation.
#' @param LNF.rho Between latent noise factor correlation.
#' @param LSNF.rho Correlation between latent signal and noise factors.
#' @param S.per.LF Number of observed features per latent factor.
#' @param Y.psi Focal trait unique variance (i.e. 1 minus the communality).
#' @param L.min Minimum absolute loading for the observed features.
#' @param L.max Maximum absolute loading for the observed features.
#' @param S.sg2 Observed feature genetic variance.
#' @param S.se2 Observed feature residual variance.
#' @param Y.sg2 Focal trait genetic variance.
#' @param Y.se2 Focal trait residual variance.
#' @param resCors Should there be random residual correlations among secondary features and between sec. features and the focal trait?
#'
#' @return A list with the "real" simulated dataframe, benchmark dataframe, and the corresponding genetic and residual covariance matrices.
#' The list also contains the prediction target and simulated BLUEs for all features.
#' @export
#'
GFAsim <- function(K, r, n.LSF, n.LNF, LSF.rho = 0, LNF.rho = 0, LSNF.rho = 0, S.per.LF = 10, Y.psi,
                   L.min = 0.3, L.max = 0.8, S.sg2, S.se2, Y.sg2, Y.se2, resCors = TRUE) {

  n.LF <- n.LSF + n.LNF

  # Factor covariance phi:
  phi <- diag(n.LF)

  phi[1:n.LSF, 1:n.LSF][upper.tri(phi[1:n.LSF, 1:n.LSF])] <-
    phi[1:n.LSF, 1:n.LSF][lower.tri(phi[1:n.LSF, 1:n.LSF])] <- LSF.rho

  phi[(n.LSF + 1):nrow(phi), (n.LSF + 1):ncol(phi)][upper.tri(phi[(n.LSF + 1):nrow(phi), (n.LSF + 1):ncol(phi)])] <-
    phi[(n.LSF + 1):nrow(phi), (n.LSF + 1):ncol(phi)][lower.tri(phi[(n.LSF + 1):nrow(phi), (n.LSF + 1):ncol(phi)])] <- LNF.rho

  phi[1:n.LSF, (n.LSF + 1):ncol(phi)] <- phi[(n.LSF + 1):nrow(phi), 1:n.LSF] <- LSNF.rho

  # Unique variances, assuming loadings between L.min and L.max for the observed secondary features:
  psi <- diag(x = c(stats::runif(n.LF * S.per.LF, 1 - L.max^2, 1 - L.min^2), Y.psi))

  # Communalities:
  comms <- 1 - diag(psi)

  # Loadings matrix:
  L <- matrix(0, n.LF * S.per.LF + 1, n.LF)
  colnames(L) <- paste0("F", 1:n.LF)
  rownames(L) <- c(paste0("S", 1:(n.LF * S.per.LF)), "Y")

  # Secondary features:
  for (i in 1:n.LF) {
    L[((i - 1) * S.per.LF + 1):(i * S.per.LF), i] <-
      sign(stats::runif(S.per.LF, -1, 1)) * sqrt(comms[((i - 1) * S.per.LF + 1):(i * S.per.LF)])
  }

  # Focal trait:
  L[nrow(L), 1:n.LSF] <- sqrt((1 - Y.psi) / n.LSF)

  # Genetic correlation matrix:
  Rg <- L %*% phi %*% t(L) + psi

  # Check:
  stopifnot(`Diagonal entries of simulated Rg are not all 1!` = all.equal(as.numeric(diag(Rg)), rep(1, n.LF * S.per.LF + 1)))

  # Genetic and residual standard deviations of the secondary features and focal trait:
  SD.G <- diag(x = c(rep(sqrt(S.sg2), n.LF * S.per.LF), sqrt(Y.sg2)))
  SD.E <- diag(x = c(rep(sqrt(S.se2), n.LF * S.per.LF), sqrt(Y.se2)))

  Sg <- outer(diag(SD.G), diag(SD.G)) * Rg
  if (resCors) {
    Re <- cov2cor(clusterGeneration::genPositiveDefMat(dim = n.LF * S.per.LF + 1, covMethod = "eigen")$Sigma)
    Se <- outer(diag(SD.E), diag(SD.E)) * Re
  } else {
    Se <- diag(x = 1 - diag(Sg))
    Re <- stats::cov2cor(Se)
  }

  # G = LX + E
  # Simulating X (factor scores):
  X <- simMVData(rowCov = K, colCov = phi)
  LX <- X %*% t(L)

  # Simulating E:
  E <- simMVData(rowCov = K, colCov = psi)

  G <- LX + E

  # Y = G + e
  # Y = LX + E + e
  # Simulating e on correlation scale:
  if (resCors) {
    e <- simMVData(rowCov = diag(nrow(K) * r), colCov = Re)
  } else {
    e <- simMVData(rowCov = diag(nrow(K) * r), colCov = diag(n.LF * S.per.LF + 1))
  }

  # Scaling everything to the covariance scale and adding it all up:
  Y <- kronecker((G %*% SD.G), matrix(rep(1, r), ncol = 1)) + e %*% SD.E

  gens <- rep(paste0("g", 1:nrow(K)), each = r)
  rownames(G) <- unique(gens)
  data.real <- data.frame(G = gens, Y)
  names(data.real)[2:ncol(data.real)] <- rownames(L)
  data.real$G <- as.factor(data.real$G)

  data.benchmark <- data.frame(G = unique(gens), cbind(X, G[, ncol(G)]))
  names(data.benchmark)[2:ncol(data.benchmark)] <- c(colnames(L), "Y")
  data.benchmark$G <- as.factor(data.benchmark$G)

  pred.target <- G[, ncol(G)]

  Rg.benchmark <- rbind(cbind(phi, matrix(0, nrow(phi))),
                        matrix(0, 1, (ncol(phi) + 1)))
  Rg.benchmark[nrow(Rg.benchmark), ncol(Rg.benchmark)] <- 1
  Rg.benchmark[nrow(Rg.benchmark), 1:n.LSF] <- Rg.benchmark[1:n.LSF, ncol(Rg.benchmark)] <- sqrt((1 - Y.psi) / n.LSF)
  Sg.benchmark <- outer(c(rep(1, n.LF), sqrt(Y.sg2)), c(rep(1, n.LF), sqrt(Y.sg2))) * Rg.benchmark

  Se.benchmark <- matrix(0, n.LF + 1, n.LF + 1)
  Se.benchmark[nrow(Se.benchmark), ncol(Se.benchmark)] <- Y.se2

  return(list(Sg.real = Sg,
              Se.real = Se,
              Sg.bm = Sg.benchmark,
              Se.bm = Se.benchmark,
              data.real = data.real,
              data.bm = data.benchmark,
              pred.target = pred.target,
              G = G))
}
