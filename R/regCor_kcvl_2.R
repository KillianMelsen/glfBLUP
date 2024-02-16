regCor_kcvl_2 <- function(penalty, vecs) {

  # Number of folds
  n.folds <- length(vecs)

  # Start the cross-validated LL at 0
  cvLL <- 0

  # Increase cvLL for each fold
  for (i in 1:n.folds){

    cvLL <- cvLL + vecs[[i]]$nf * regCor_ll_2(penalty = penalty,
                                              R.eigenvalues = vecs[[i]]$R.EV,
                                              A.diag = vecs[[i]]$A.diag)
  }

  # Return the mean cross-validated LL
  return(cvLL / n.folds)

}
