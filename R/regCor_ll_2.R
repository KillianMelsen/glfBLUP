regCor_ll_2 <- function(penalty, R.eigenvalues, A.diag) {

  ll <- log(prod((1 - penalty) * R.eigenvalues + penalty)) + sum(A.diag / ((1 - penalty) * R.eigenvalues + penalty))
  return(ll)

}
