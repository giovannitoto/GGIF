# ---------------------------------------------------------------------------- #

#' Compute the posterior mean of one or more parameters
#'
#' @param out_MCMC A list containing the results of the Adaptive Gibbs Sampler.
#' @param parameters A vector containing the names of the parameters for which you want to compute the posterior mean. The possible valid strings are \code{mu}, \code{bmu}, \code{sigmacol}, \code{omega} and \code{omega_inv}. Default is \code{"all"}, which is equivalent to writing \code{c("mu", "bmu", "sigmacol","omega","omega_inv")}.
#' @param columns Default is \code{"k"}.
#'
#' @return A list containing the posterior means of the parameters specified in \code{parameters}.
posterior_mean <- function(out_MCMC, parameters = "all", columns = "k") {
  # remove invalid strings from 'parameters'
  valid_parameters <- c("mu", "bmu", "sigmacol", "omega", "omega_inv")
  if("all" %in% parameters) {
    parameters <- valid_parameters
  } else {
    parameters <- intersect(parameters, valid_parameters)
  }
  # start time
  t0 <- proc.time()
  # prepare output
  output <- list()
  # check whether all necessary variables are available in out_MCMC
  for (par in intersect(parameters, valid_parameters[1:3])) {
    if(par %in% names(out_MCMC)) {
      # compute the posterior mean of the parameter 'par'
      output[[par]] <- Reduce("+", out_MCMC[[par]]) / length(out_MCMC[[par]])
    } else {
      stop(paste("'", par, "' not stored in the MCMC output.", sep=""))
    }
  }
  # compute the posterior mean of 'omega' and/or 'omega_inv'
  if (any(c("omega", "omega_inv") %in% parameters)) {
    # check the value of 'columns'
    if(!(columns %in% c("k", "kstar"))) {
      stop("'columns' not valid: it must be 'k' or 'kstar'.")
    }
    if(check_list(out_MCMC, c("lambda", "numFactors", "sigmacol")) == FALSE) {
      stop("'lambda', 'numFactors' and 'sigmacol' must be stored in the MCMC output.")
    } else {
      # number of iterations
      t <- length(out_MCMC$sigmacol)
      p <- nrow(out_MCMC$lambda[[1]])
      out_MCMC$lambda <- rescale_parameter(out_MCMC = out_MCMC, parameters = "lambda", columns = columns)$lambda
      if ("omega" %in% parameters) {
        output[["omega"]] <- 0
        for (i in 1:t) {
          output[["omega"]] <- output[["omega"]] + omega_inversion(out_MCMC$lambda[[i]], out_MCMC$sigmacol[[i]])
        }
        output[["omega"]] <- output[["omega"]] / t
      }
      if ("omega_inv" %in% parameters) {
        output[["omega_inv"]] <- 0
        for (i in 1:t) {
          output[["omega_inv"]] <- output[["omega_inv"]] + tcrossprod(out_MCMC$lambda[[i]]) + diag(1 / out_MCMC$sigmacol[[i]])
        }
        output[["omega_inv"]] <- output[["omega_inv"]] / t
      }
    }
  }
  # stop if no valid string is provided as input
  if(length(output) == 0) {
    stop("No valid parameter provided as input.")
  }
  # print runtime
  print(proc.time() - t0, "\n")
  # return a list
  return(output)
}

# ---------------------------------------------------------------------------- #
