% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lposterior.R
\name{lposterior}
\alias{lposterior}
\title{Log-posterior probability}
\usage{
lposterior(
  out_MCMC,
  frac_sampled = 1,
  samples = NULL,
  columns = "k",
  parameters = "all",
  seed = 28
)
}
\arguments{
\item{out_MCMC}{A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{AGS_SIS}}.}

\item{frac_sampled}{A value in (0,1] specifying the fraction of the iterations to evaluate. Default is 1, that is all the iterations.}

\item{samples}{A vector of integers between 1 and the length of the MCMC chain which specifies the iterations to evaluate. Default is all the iterations.}

\item{columns}{A string specifying whether to consider only active factors (\code{"kstar"}) or both active and inactive factors (\code{"k"}). Default is \code{"k"}.}

\item{parameters}{A vector containing the names of the parameters for which you want to compute the posterior mean. The possible valid strings are \code{"beta"}, \code{"eta"}, \code{"lambda"}, \code{"mu"} and \code{"bmu"}. Default is \code{all}, which is equivalent to writing \code{c("beta", "eta", "lambda", "mu", "bmu")}.}

\item{seed}{Seed. Default is 28.}
}
\value{
A list with the following elements:
\itemize{
\item \code{sampled}: indices of the iterations evaluated by the function.
\item \code{lposterior}: log-posterior probabilities of the iterations evaluated by the function.
\item \code{lposterior_max}: maximum log-posterior probability.
\item \code{iteration_max}: index of the iteration in which the log-posterior probability reaches its maximum.
\item \code{beta_max}: value of \eqn{\beta} at \code{iteration_max}.
\item \code{eta_max}: value of \eqn{\eta} at \code{iteration_max}.
\item \code{lambda_max}: value of \eqn{\Lambda} at \code{iteration_max}.
\item \code{mu}: value of \eqn{\mu} at \code{iteration_max}.
\item \code{bmu}: value of \eqn{b_{\mu}} at \code{iteration_max}.
}
If \code{columns="k"}, \code{beta_max}, \code{eta_max} and \code{lambda_max} contain only columns linked to active factors; if \code{columns="kstar"}, they contain columns linked to both active and inactive factors
}
\description{
Compute the log-posterior probabilities of part or all the MCMC iterations.
To speed up computations, it is possible to evaluate only a random fraction of the iterations, using the \code{frac_sampled} argument, or to define the iterations of interest, using the \code{samples} argument.
}
\seealso{
This function is applied to an output of \code{\link{AGS_SIS}}.
}
