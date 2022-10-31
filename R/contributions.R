# ---------------------------------------------------------------------------- #

#' Contributions
#'
#' Compute the sorted contributions \eqn{C^*_1,\ldots,C^*_{k_{max}}} for each iteration.
#'
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{AGS_SIS}}.
#' @param reference An integer specifying the reference iteration used to sort the contributions of each iteration; if \code{reference=NULL}, the contributions are sorted using the Frobenius norm. Default is \code{NULL}.
#' @param verbose Logical: if \code{TRUE}, print the number of active factors every 50 iterations. Default is \code{TRUE}.
#'
#' @return A list containing draws from the posterior distribution of \eqn{C^*_1,\ldots,C^*_{k_{max}}}.
#'
#' @seealso This function is applied to an output of \code{\link{AGS_SIS}}.
#'
#' @export
contributions <- function(out_MCMC, reference = NULL, verbose = TRUE) {
  # -------------------------------------------------------------------------- #
  # check whether all necessary variables are available in out_MCMC
  required_variables <- c("numFactors", "eta", "lambda")
  if(check_list(out_MCMC, required_variables) == FALSE) {
    stop(paste(paste(c(required_variables), collapse = ", "), "must be stored in the MCMC output."))
  }
  # -------------------------------------------------------------------------- #
  if(!is.null(reference)) {
    if((length(reference) != 1) || !is.numeric(reference) || (reference != round(reference)) || reference > length(out_MCMC$numFactors)) {
      stop("'reference' not valid: it must be NULL or an integer between 1 and the number of iterations..")
    }
    # number of active factors of reference iteration
    kstar_ref <- out_MCMC$numFactors[reference]
    # compute contributions NOT in order of reference iteration
    C_list <- list()
    for (h in 1:kstar_ref) {
      C_list[[paste0("C", h)]] <- tcrossprod(out_MCMC$eta[[reference]][, h], out_MCMC$lambda[[reference]][, h])
    }
    # compute order of the contributions using Frobenius norm
    C_order <- order(sapply(C_list, function(C) norm(C, type = "F")), decreasing = TRUE)
    # order contributions of reference iteration
    C_ref <- list()
    for (h in 1:kstar_ref) {
      C_ref[[paste0("C", h)]] <- C_list[[paste0("C", C_order[h])]]
    }
  }
  # -------------------------------------------------------------------------- #
  n <- nrow(out_MCMC$eta[[1]])
  p <- nrow(out_MCMC$lambda[[1]])
  kstar_max <- max(out_MCMC$numFactors)
  # -------------------------------------------------------------------------- #
  output <- list()
  # output <- list("reference" = reference)
  for (h in 1:kstar_max) {
    output[[paste0("C", h)]] <- vector("list", length = length(out_MCMC$numFactors))
  }
  # -------------------------------------------------------------------------- #
  for (it in 1:length(out_MCMC$numFactors)) {
    # number of active factors
    kstar <- out_MCMC$numFactors[it]
    # compute contributions NOT in order
    C_list <- list()
    for (h in 1:kstar) {
      C_list[[paste0("C", h)]] <- tcrossprod(out_MCMC$eta[[it]][, h], out_MCMC$lambda[[it]][, h])
    }
    # add matrices of zeros if kstar is less than the maximum kstar observed
    if (kstar < kstar_max) {
      for (h in (kstar+1):kstar_max) C_list[[paste0("C", h)]] <- matrix(0, nrow = n, ncol = p)
    }
    # compute order using Frobenius norm or reference iteration
    if(is.null(reference)) {
      # compute order of the contributions using Frobenius norm
      C_order <- order(sapply(C_list, function(C) norm(C, type = "F")), decreasing = TRUE)
      # save contributions
      for (h in 1:kstar_max) {
        output[[paste0("C", h)]][[it]] <- C_list[[paste0("C", C_order[h])]]
      }
    } else {
      # compute order of the contributions using reference iteration
      for (h in 1:kstar_ref) {
        min_idx <- names(which.min(sapply(C_list, function(C) norm(C - C_ref[[paste0("C", h)]], type = "F"))))
        output[[paste0("C", h)]][[it]] <- C_list[[min_idx]]
        C_list[[min_idx]] <- NULL
      }
      if(kstar_ref < kstar_max) {
        # the remaining contributions are sorted using Frobenius norm
        for (h in (kstar_ref+1):kstar_max) {
          min_idx <- names(which.min(sapply(C_list, function(C) norm(C, type = "F"))))
          output[[paste0("C", h)]][[it]] <- C_list[[paste0("C", min_idx)]]
          C_list[[paste0("C", min_idx)]] <- NULL
        }
      }
    }
    if(verbose == TRUE && it %% 50 == 0) cat(it, ":", kstar, "active factors\n")
  }
  # -------------------------------------------------------------------------- #
  return(output)
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
