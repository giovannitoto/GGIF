# ---------------------------------------------------------------------------- #

#' Check if all the required variables are a list
#'
#' @param list_to_check A list.
#' @param required_variables A vector containing the names of the variables which must be in \code{list_to_check}.
#'
#' @return It returns \code{TRUE} if \code{list_to_check} contains all the required variables specified in \code{required_variables}.
check_list <- function(list_to_check, required_variables) {
  if(all(required_variables %in% names(list_to_check))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# ---------------------------------------------------------------------------- #

#' Inversion of omega
#'
#' This function is called in 'posterior_mean' and 'lposterior_function'.
#'
#' @param lambda A pxk matrix of factorial weights.
#' @param sigmai A p-dimensional vector of \code{1/sigma^2}.
#'
#' @return A pxp matrix.
omega_inversion <- function(lambda, sigmai) {
  p <- dim(lambda)[1]
  k <- dim(lambda)[2]
  # Woodbury matrix identity
  lambda_sigmai <- crossprod(lambda, diag(sigmai))
  capac <- diag(k) + lambda_sigmai %*% lambda
  omega_inv <- diag(sigmai) - diag(sigmai) %*% lambda %*% solve(capac) %*% lambda_sigmai
  return(omega_inv)
}

# ---------------------------------------------------------------------------- #
