# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Adaptive Gibbs Sampler (AGS)
#' Implementation in C++ of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
#'
#' @param alpha \eqn{Beta(1, \code{alpha})} for pi.
#' @param as,bs \eqn{\sigma^{-2} \sim Gamma(\code{as}, \code{bs})}.
#' @param a_y,a_yp1 Two nxp matrices containing the lower and upper endpoints of the truncated normal distribution.
#' @param a_theta,b_theta DA SCRIVERE
#' @param Beta A qxk matrix of shrinkage coefficients.
#' @param b0,b1 Positive constants for the adaptive probability.
#' @param burn An integer number of burning iterations (number of iterations to discard).
#' @param b_mu A qxc matrix of coefficients.
#' @param c An integer number of covariates.
#' @param d DA SCRIVERE
#' @param eta A nxk matrix of coefficients.
#' @param k An integer number of factors.
#' @param kmax Maximum number of latent factors.
#' @param kstar An integer number of active factors.
#' @param Lambda A pxk matrix of coefficients.
#' @param Lambda_star A pxk matrix of DA SCRIVERE
#' @param logit DA SCRIVERE
#' @param mu A pxc matrix of coefficients.
#' @param n An integer number of observations.
#' @param nrun An integer number of iterations. Default is 100.
#' @param out List of variables to save during the Adaptive Gibbs Sampler.
#' @param p An integer number of variables.
#' @param Phi DA SCRIVERE
#' @param Plam DA SCRIVERE
#' @param prec_b DA SCRIVERE
#' @param prec_mu DA SCRIVERE
#' @param pred DA SCRIVERE
#' @param prob DA SCRIVERE
#' @param ps DA SCRIVERE
#' @param p_constant Factor probability constant.
#' @param q_mean An integer number of meta-covariates.
#' @param q_cov An integer number of meta-covariates.
#' @param rho DA SCRIVERE
#' @param sd_b Standard deviation for \eqn{b}.
#' @param sd_beta Standard deviation for \eqn{\beta}.
#' @param sd_mu Standard deviation for \eqn{\mu}.
#' @param sp DA SCRIVERE
#' @param start_adapt An integer number of iterations before adaptation.
#' @param thin An integer thinning value (number of iterations to skip between saving iterations).
#' @param uu DA SCRIVERE
#' @param v DA SCRIVERE
#' @param verbose Boolean: if \code{true}, print the number of active factors every 10 iterations.
#' @param w DA SCRIVERE
#' @param W A nxc matrix of covariates.
#' @param Wnull Boolean: if \code{true}, the matrix \code{W} of covariates is not considered.
#' @param X_mean A pxq_mean matrix of meta-covariates.
#' @param X_cov A pxq_cov matrix of meta-covariates.
#'
#' @return A list.
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
Rcpp_AGS_SIS <- function(alpha, as, a_y, a_yp1, a_theta, Beta, bs, b0, b1, burn, b_mu, b_theta, c, d, eta, k, kmax, kstar, Lambda, Lambda_star, logit, mu, n, nrun, out, p, Phi, Plam, prec_b, prec_mu, pred, prob, ps, p_constant, q_mean, q_cov, rho, sd_b, sd_beta, sd_mu, sp, start_adapt, thin, uu, v, verbose, w, W, Wnull, X_mean, X_cov) {
    .Call(`_GGIF_Rcpp_AGS_SIS`, alpha, as, a_y, a_yp1, a_theta, Beta, bs, b0, b1, burn, b_mu, b_theta, c, d, eta, k, kmax, kstar, Lambda, Lambda_star, logit, mu, n, nrun, out, p, Phi, Plam, prec_b, prec_mu, pred, prob, ps, p_constant, q_mean, q_cov, rho, sd_b, sd_beta, sd_mu, sp, start_adapt, thin, uu, v, verbose, w, W, Wnull, X_mean, X_cov)
}

#' Update b_mu in the Adaptive Gibbs Sampler
#'
#' @param X A pxq matrix.
#' @param prec_mu A positive number.
#' @param prec_b A positive number.
#' @param mu A pxc matrix.
#' @param q An integer.
#' @param c An integer.
#'
#' @return A qxc matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_bmu <- function(X, prec_mu, prec_b, mu, q, c) {
    .Call(`_GGIF_update_bmu`, X, prec_mu, prec_b, mu, q, c)
}

#' Update the hth row of mu in the Adaptive Gibbs Sampler
#'
#' @param j An integer.
#' @param Qbet A cxc matrix.
#' @param W A nxc matrix.
#' @param Z_res A nxp matrix.
#' @param ps A p-dimensional vector.
#' @param b_mu A cxq matrix.
#' @param Xcov A pxq matrix.
#' @param c An integer.
#'
#' @return A 1xc matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_mu <- function(j, Qbet, W, Z_res, ps, b_mu, Xcov, c) {
    .Call(`_GGIF_update_mu`, j, Qbet, W, Z_res, ps, b_mu, Xcov, c)
}

#' Update eta in the Adaptive Gibbs Sampler
#'
#' @param Lambda A pxk matrix.
#' @param ps A p-dimensional vector.
#' @param k An integer.
#' @param Z A nxp matrix.
#' @param n An integer.
#'
#' @return A nxk matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_eta <- function(Lambda, ps, k, Z, n) {
    .Call(`_GGIF_update_eta`, Lambda, ps, k, Z, n)
}

#' Update the hth column of beta in the Adaptive Gibbs Sampler
#'
#' @param h An integer.
#' @param Xcov A pxq matrix.
#' @param Dt A pxk matrix.
#' @param Bh_1 A qxq matrix.
#' @param Phi_L A pxk matrix.
#' @param q An integer.
#'
#' @return A qx1 matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_beta <- function(h, Xcov, Dt, Bh_1, Phi_L, q) {
    .Call(`_GGIF_update_beta`, h, Xcov, Dt, Bh_1, Phi_L, q)
}

#' Update the jth row of Lambda_star in the Adaptive Gibbs Sampler
#'
#' @param j An integer.
#' @param etarho A kxn matrix.
#' @param Phi A pxk matrix.
#' @param Plam A kxk matrix.
#' @param ps A p-dimensional vector;
#' @param Z A nxp matrix.
#' @param k An integer.
#'
#' @return A 1xk matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_Lambda_star <- function(j, etarho, Phi, Plam, ps, Z, k) {
    .Call(`_GGIF_update_Lambda_star`, j, etarho, Phi, Plam, ps, Z, k)
}

#' Update the hth element of d in the Adaptive Gibbs Sampler
#'
#' @param h An integer.
#' @param Phi A pxk matrix.
#' @param p An integer.
#' @param n An integer.
#' @param rho A k-dimensional vector.
#' @param eta A nxk matrix.
#' @param lambdastar A pxk matrix.
#' @param Z A nxp matrix.
#' @param ps A p-dimensional vector.
#' @param k An integer.
#' @param w A k-dimensional vector.
#'
#' @return An integer in 1, ..., k.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_d <- function(h, Phi, p, n, rho, eta, lambdastar, Z, ps, k, w) {
    .Call(`_GGIF_update_d`, h, Phi, p, n, rho, eta, lambdastar, Z, ps, k, w)
}

#' Update Phi in the Adaptive Gibbs Sampler
#'
#' @param rho A k-dimensional vector.
#' @param logit A pxk matrix.
#' @param p_constant A number in (0,1).
#' @param p An integer.
#' @param n An integer.
#' @param eta A nxk matrix.
#' @param lambdastar A pxk matrix.
#' @param Phi A pxk matrix.
#' @param Z A nxp matrix.
#' @param ps A p-dimensional vector.
#' @param k An integer.
#'
#' @return A pxk matrix.
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
update_Phi <- function(rho, logit, p_constant, p, n, eta, lambdastar, Phi, Z, ps, k) {
    .Call(`_GGIF_update_Phi`, rho, logit, p_constant, p, n, eta, lambdastar, Phi, Z, ps, k)
}

rcpp_pgdraw <- function(b, c) {
    .Call(`_GGIF_rcpp_pgdraw`, b, c)
}

samplepg <- function(z) {
    .Call(`_GGIF_samplepg`, z)
}

#' Sample from a truncated normal distribution. Samples are drawn
#' componentwise, so each component of the vector is allowed its own
#' mean, standard deviation, and upper and lower limits. The components
#' are assumed to be independent.
#'
#' @param y_lower \code{n x p} matrix of lower endpoints
#' @param y_upper \code{n x p} matrix of upper endpoints
#' @param mu \code{n x p} matrix of conditional expectations
#' @param sigma \code{p x 1} vector of conditional standard deviations
#' @param u_rand \code{n x p} matrix of uniform random variables
#'
#' @return z_star \code{n x p} draw from the truncated normal distribution
#'
#' @note This function uses \code{Rcpp} for computational efficiency.
#'
truncnorm_lg <- function(y_lower, y_upper, mu, sigma, u_rand) {
    .Call(`_GGIF_truncnorm_lg`, y_lower, y_upper, mu, sigma, u_rand)
}

