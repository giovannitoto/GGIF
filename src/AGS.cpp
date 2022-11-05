#include <RcppArmadillo.h>
#include "AGS_update.h"
#include "helper_functions.h"
#include "truncnorm_lg.h"
#include "rcpp_pgdraw.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

//' Adaptive Gibbs Sampler (AGS)
//' Implementation in C++ of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
//'
//' @param alpha \eqn{Beta(1, \code{alpha})} for pi.
//' @param as,bs \eqn{\sigma^{-2} \sim Gamma(\code{as}, \code{bs})}.
//' @param a_y,a_yp1 Two nxp matrices containing the lower and upper endpoints of the truncated normal distribution.
//' @param a_theta,b_theta DA SCRIVERE
//' @param Beta A qxk matrix of shrinkage coefficients.
//' @param b0,b1 Positive constants for the adaptive probability.
//' @param burn An integer number of burning iterations (number of iterations to discard).
//' @param b_mu A qxc matrix of coefficients.
//' @param c An integer number of covariates.
//' @param d DA SCRIVERE
//' @param eta A nxk matrix of coefficients.
//' @param k An integer number of factors.
//' @param kmax Maximum number of latent factors.
//' @param kstar An integer number of active factors.
//' @param Lambda A pxk matrix of coefficients.
//' @param Lambda_star A pxk matrix of DA SCRIVERE
//' @param logit DA SCRIVERE
//' @param mu A pxc matrix of coefficients.
//' @param n An integer number of observations.
//' @param nrun An integer number of iterations. Default is 100.
//' @param out List of variables to save during the Adaptive Gibbs Sampler.
//' @param p An integer number of variables.
//' @param Phi DA SCRIVERE
//' @param Plam DA SCRIVERE
//' @param prec_b DA SCRIVERE
//' @param prec_mu DA SCRIVERE
//' @param pred DA SCRIVERE
//' @param prob DA SCRIVERE
//' @param ps DA SCRIVERE
//' @param p_constant Factor probability constant.
//' @param q_mean An integer number of meta-covariates.
//' @param q_cov An integer number of meta-covariates.
//' @param rho DA SCRIVERE
//' @param sd_b Standard deviation for \eqn{b}.
//' @param sd_beta Standard deviation for \eqn{\beta}.
//' @param sd_mu Standard deviation for \eqn{\mu}.
//' @param sp DA SCRIVERE
//' @param start_adapt An integer number of iterations before adaptation.
//' @param thin An integer thinning value (number of iterations to skip between saving iterations).
//' @param uu DA SCRIVERE
//' @param v DA SCRIVERE
//' @param verbose Boolean: if \code{true}, print the number of active factors every 10 iterations.
//' @param w DA SCRIVERE
//' @param W A nxc matrix of covariates.
//' @param Wnull Boolean: if \code{true}, the matrix \code{W} of covariates is not considered.
//' @param X_mean A pxq_mean matrix of meta-covariates.
//' @param X_cov A pxq_cov matrix of meta-covariates.
//'
//' @return A list.
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
Rcpp::List Rcpp_AGS_SIS(double alpha, double as, arma::mat a_y, arma::mat a_yp1, double a_theta,
          arma::mat Beta, double bs, double b0, double b1, int burn, arma::mat b_mu, double b_theta,
          int c,
          arma::vec d,
          arma::mat eta,
          int k, int kmax, int kstar,
          arma::mat Lambda, arma::mat Lambda_star, arma::mat logit,
          arma::mat mu,
          int n, int nrun,
          Rcpp::List out,
          int p, arma::mat Phi, arma::mat Plam, double prec_b, double prec_mu, arma::mat pred, arma::vec prob, arma::vec ps, double p_constant,
          int q_mean, int q_cov,
          arma::vec rho,
          double sd_b, double sd_beta, double sd_mu, int sp, int start_adapt,
          int thin,
          arma::vec uu,
          arma::vec v, bool verbose,
          arma::vec w, arma::mat W, bool Wnull,
          arma::mat X_mean, arma::mat X_cov) {
  // ---------------------------------------------------------------------------
  // output
  Rcpp::List MU(sp);
  Rcpp::List BMU(sp);
  Rcpp::List BETA(sp);
  Rcpp::List ETA(sp);
  Rcpp::List LAMBDA(sp);
  Rcpp::List SIG(sp);
  arma::vec K(sp);
  // ---------------------------------------------------------------------------
  int it, i, j, h;
  int ind = 0;
  for (it = 0; it < nrun; it++) {
    if(verbose && (it + 1) % 50 == 0) {
      Rcout << it + 1 << " : " << k << " active factors\n";
    }
    // -------------------------------------------------------------------------
    //Rcout << "1  ";
    // 1 - update Z
    arma::mat Zmean(n, p);
    if(Wnull) {
      Zmean = eta * trans(Lambda);
    } else {
      Zmean = eta * trans(Lambda) + W * trans(mu);
    }
    arma::mat n_unif = runif_mat(n, p, 0, 1);
    arma::mat Z = truncnorm_lg(log(a_y), log(a_yp1), Zmean, sqrt(1 / ps), n_unif);
    // -------------------------------------------------------------------------
    if(!Wnull) {
      // -----------------------------------------------------------------------
      //Rcout << "2  ";
      // 2 - update b_mu
      b_mu = update_bmu(X_mean, prec_mu, prec_b, mu, q_mean, c);
      // -----------------------------------------------------------------------
      //Rcout << "3  ";
      // 3 - update mu
      arma::mat Z_res = Z - eta * trans(Lambda);
      arma::mat Qbet(c, c);
      if(c > 1) {
        Qbet = arma::diagmat(prec_mu * arma::ones(c)) + trans(W) * W;
      } else {
        Qbet = prec_mu + trans(W) * W;
      }
      if (!Qbet.is_symmetric()) {
        Rcout << "update_mu" << "\t";
        Qbet = arma::symmatu(Qbet);
      }
      for (j = 0; j < p; j++) {
        mu.row(j) = update_mu(j, Qbet, W, Z_res, ps, b_mu, X_mean, c);
      }
      Z = Z - W * trans(mu);
      // -----------------------------------------------------------------------
    }
    //Rcout << "4  ";
    // 4 update eta
    eta = update_eta(Lambda, ps, k, Z, n);
    // -------------------------------------------------------------------------
    //Rcout << "5  ";
    // 5 - update Sigma
    arma::mat Z_res = Z - eta * trans(Lambda);
    for (j = 0; j < p; j++) {
      ps(j) = R::rgamma(as + 0.5 * n, 1 / (bs + 0.5 * arma::accu(arma::pow(Z_res.col(j), 2))));
    }
    // -------------------------------------------------------------------------
    //Rcout << "6  ";
    // 6 - update beta
    pred = X_cov * Beta;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    // 6.1 Update phi_L
    arma::mat Phi_L = arma::ones(p, k);
    arma::uvec Phi0 = arma::find(Phi == 0);
    arma::vec logit_phi0 = logit.elem(Phi0);
    arma::uvec which_zero = arma::randu(logit_phi0.n_elem) < (1 - logit_phi0) / (1 - logit_phi0 * p_constant);
    Phi_L.elem(Phi0.elem(arma::find(which_zero))) -= 1;
    // 6.2 Polya gamma
    arma::mat Dt(p, k);
    for (j = 0; j < k; j++) {
      for (i = 0; i < p; i++) {
        // https://cran.r-project.org/web/packages/pgdraw/index.html
        Dt(i, j) = samplepg(pred(i, j));
      }
    }
    // 6.3 Update beta_h
    arma::mat Bh_1 = arma::diagmat(arma::ones(q_cov) / pow(sd_beta, 2));
    for (h = 0; h < k; h++) {
      Beta.col(h) = trans(update_beta(h, X_cov, Dt, Bh_1, Phi_L, q_cov));
    }
    // -------------------------------------------------------------------------
    //Rcout << "7  ";
    // 7 - update Lambda and Lambda_star
    arma::mat etarho = trans(eta.each_row() % trans(rho));
    for (j = 0; j < p; j++) {
      Lambda_star.row(j) = update_Lambda_star(j, etarho, Phi, Plam, ps, Z, k);
    }
    Lambda = (Lambda_star.each_row() % trans(sqrt(rho))) % Phi;
    // -------------------------------------------------------------------------
    //Rcout << "8  ";
    // 8.1 - update d
    for (h = 0; h < k; h++) {
      d(h) = update_d(h, Phi, p, n, rho, eta, Lambda_star, Z, ps, k, w);
    }
    rho = arma::ones(k);
    rho.elem(find(d <= arma::linspace<arma::vec>(0, k - 1, k))) -= 1;
    // 8.2
    arma::vec Plam_diag(k);
    for (h = 0; h < k; h++) {
      Plam_diag(h) = R::rgamma(a_theta + 0.5 * p, 1 / (b_theta + 0.5 * arma::accu(arma::pow(Lambda_star.col(h), 2))));
    }
    Plam = arma::diagmat(Plam_diag);
    // 8.3
    for (h = 0; h < k - 1; h++) {
      v(h) = R::rbeta(1 + arma::accu(d == h), alpha + arma::accu(d > h));
    }
    v(k - 1) = 1;
    w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
    // -------------------------------------------------------------------------
    //Rcout << "9  ";
    // 9 - update Phi
    pred = X_cov * Beta;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    Phi = update_Phi(rho, logit, p_constant, p, n, eta, Lambda_star, Phi, Z, ps, k);
    // -------------------------------------------------------------------------
    //Rcout << "save  ";
    // save sampled values (after burn-in period)
    if ((it + 1) % thin == 0 && (it + 1) > burn) {
      if (!Wnull) {
        if(out.containsElementNamed("mu")) { MU[ind] = mu; }
        if(out.containsElementNamed("bmu")) { BMU[ind] = b_mu; }
      }
      if(out.containsElementNamed("beta")) { BETA[ind] = Beta; }
      if(out.containsElementNamed("eta")) { ETA[ind] = eta; }
      if(out.containsElementNamed("lambda")) { LAMBDA[ind] = Lambda; }
      if(out.containsElementNamed("sigmacol")) { SIG[ind] = ps; }
      if(out.containsElementNamed("numFactors")) { K[ind] = kstar; }
      ind += 1;
    }
    // -------------------------------------------------------------------------
    //Rcout << "adapt\n";
    // Adaptation
    if (uu(it) < prob(it) && (it + 1) > start_adapt) {
      arma::uvec active = find(d > arma::linspace<arma::vec>(0, k - 1, k));
      int kstar_new = active.n_elem;
      kstar = kstar_new;
      //Rcout << it+1 << " : "<< (kstar < k - 1) << " - " << (k < kmax) << "\n";
      if (kstar < k - 1) {
        // set truncation to kstar and subset all variables, keeping only active columns
        k = kstar + 1;
        eta = arma::join_rows(eta.cols(active), rnorm_vec(n, 0, 1));
        double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
        Plam_diag = join_elem(Plam_diag.elem(active), vartheta_k);
        Plam = arma::diagmat(Plam_diag);
        Lambda_star = arma::join_rows(Lambda_star.cols(active), rnorm_vec(p, 0, sqrt(vartheta_k)));
        Phi = arma::join_rows(Phi.cols(active), rbinom_vec(p, 1, p_constant));
        rho = join_elem(rho.elem(active), 1);
        Lambda = arma::join_rows(Lambda.cols(active), Lambda_star.col(k - 1) % Phi.col(k - 1));
        Beta = arma::join_rows(Beta.cols(active), rnorm_vec(q_cov, 0, sqrt(sd_beta)));
        w = join_elem(w.elem(active), 1 - sum(w.elem(active)));
        v = join_elem(v.elem(active), 1);
        d = join_elem(d.elem(active), k - 1);
      } else if (k < kmax) {
      // increase truncation by 1 and extend all variables, sampling from the prior/model
      k += 1;
      eta = arma::join_rows(eta, rnorm_vec(n, 0, 1));
      double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
      Plam_diag = join_elem(Plam_diag, vartheta_k);
      Plam = arma::diagmat(Plam_diag);
      Lambda_star = arma::join_rows(Lambda_star, rnorm_vec(p, 0, sqrt(vartheta_k)));
      Phi = arma::join_rows(Phi, rbinom_vec(p, 1, p_constant));
      rho = join_elem(rho, 1);
      Lambda = arma::join_rows(Lambda, Lambda_star.col(k - 1) % Phi.col(k - 1));
      Beta = arma::join_rows(Beta, rnorm_vec(q_cov, 0, sqrt(sd_beta)));
      v(k - 2) = R::rbeta(1, alpha);
      v = join_elem(v, 1);
      w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
      d = join_elem(d, k - 1);
      }
    }
    // -------------------------------------------------------------------------
  }
  // ---------------------------------------------------------------------------
  if (!Wnull) {
    if(out.containsElementNamed("mu")) { out["mu"] = MU; }
    if(out.containsElementNamed("bmu")) { out["bmu"] = BMU; }
  }
  if(out.containsElementNamed("beta")) { out["beta"] = BETA; }
  if(out.containsElementNamed("eta")) { out["eta"] = ETA; }
  if(out.containsElementNamed("lambda")) { out["lambda"] = LAMBDA; }
  if(out.containsElementNamed("sigmacol")) { out["sigmacol"] = SIG; }
  if(out.containsElementNamed("numFactors")) { out["numFactors"] = K; }
  // ---------------------------------------------------------------------------
  return out;
  // ---------------------------------------------------------------------------
}
