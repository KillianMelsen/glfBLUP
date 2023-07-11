#' regCor_ll
#'
#' \code{.regCor_ll} calculates the negative log-likelihood assuming an identity
#' matrix regularization target.
#'
#' @param S A correlation matrix based on the left out fold.
#' @param R A correlation matrix based on the other folds.
#'
#' @return The negative log-likelihood.
#' @keywords internal
#'
regCor_ll <- function(S, R) {
    return(log(det(R)) + sum(S*solve(R)))
}
