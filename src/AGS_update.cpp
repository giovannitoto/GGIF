#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/* -------------------------------------------------------------------------- */
// Functions in this file:
//  - update_bmu
//  - update_mu
//  - update_eta
//  - update_beta
//  - update_Lambda_star
//  - update_d
//  - update_Phi
/* -------------------------------------------------------------------------- */

//' Update b_mu in the Adaptive Gibbs Sampler
//'
//' @param X A pxq matrix.
//' @param prec_mu A positive number.
//' @param prec_b A positive number.
//' @param mu A pxc matrix.
//' @param q An integer.
//' @param c An integer.
//'
//' @return A qxc matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_bmu(arma::mat X, double prec_mu, double prec_b,
                           arma::mat mu, int q, int c) {
  int i, j;
  arma::mat Trsg = X * prec_mu;
  arma::mat Vgmbeta1(q, q);
  if (q > 1) {
    Vgmbeta1 = arma::diagmat(prec_b * arma::ones(q)) + trans(Trsg) * X;
  } else {
    Vgmbeta1 = prec_b + trans(Trsg) * X;
  }
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Vgmbeta1);
  eigval = arma::reverse(eigval);
  eigvec = arma::fliplr(eigvec);
  // ---------------------------------------------------------------------------
  // multiply an eigenvector by -1 if its first element is negative
  for (j = 0; j < q; j++) {
    eigvec.col(j) *= arma::sign(eigvec(0, j));
  }
  // ---------------------------------------------------------------------------
  arma::mat Tmat(q, q);
  if (arma::min(eigval) > 0.000001) {
    Tmat = trans(eigvec.each_row() % trans(sqrt(eigval)));
  } else {
    Tmat = arma::chol(Vgmbeta1);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q, R, Tmat);
  // ---------------------------------------------------------------------------
  // https://www.mathworks.com/matlabcentral/answers/83798-sign-differences-in-qr-decomposition
  // enforce positive diagonals of R
  R = arma::diagmat(arma::sign(R.diag())) * R;
  // ---------------------------------------------------------------------------
  arma::mat S = arma::inv(R);
  arma::mat M(c, q);
  for (j = 0; j < q; j++) {
    for (i = 0; i < c; i++) {
      M(i, j) = R::rnorm(0, 1);
    }
  }
  arma::mat b_mu = trans(trans(mu) * Trsg * S * trans(S) + M * trans(S));
  return b_mu;
}

/* -------------------------------------------------------------------------- */

//' Update the hth row of mu in the Adaptive Gibbs Sampler
//'
//' @param j An integer.
//' @param Qbet A cxc matrix.
//' @param W A nxc matrix.
//' @param Z_res A nxp matrix.
//' @param ps A p-dimensional vector.
//' @param b_mu A cxq matrix.
//' @param Xcov A pxq matrix.
//' @param c An integer.
//'
//' @return A 1xc matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_mu(int j, arma::mat Qbet, arma::mat W, arma::mat Z_res,
                    arma::vec ps, arma::mat b_mu, arma::mat Xcov, int c) {
  arma::mat Lbet = trans(arma::chol(Qbet));
  arma::vec bbet = trans(W) * Z_res.col(j) + ps(j) * trans(b_mu) * trans(Xcov.row(j));
  // mean
  arma::mat vbet = solve(arma::trimatl(Lbet), bbet);
  arma::mat mbet = solve(arma::trimatu(trans(Lbet)), vbet);
  // var
  arma::vec zbet(c);
  for (int i = 0; i < c; i++) {
    zbet(i) = R::rnorm(0, 1);
  }
  arma::mat ybet = solve(trimatu(trans(Lbet)), zbet);
  arma::mat muj = trans(ybet + mbet);
  return muj;
}

/* -------------------------------------------------------------------------- */

//' Update eta in the Adaptive Gibbs Sampler
//'
//' @param Lambda A pxk matrix.
//' @param ps A p-dimensional vector.
//' @param k An integer.
//' @param Z A nxp matrix.
//' @param n An integer.
//'
//' @return A nxk matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_eta(arma::mat Lambda, arma::vec ps, int k, arma::mat Z, int n) {
  int i, j;
  arma::mat Lmsg = Lambda.each_col() % ps;
  arma::mat Veta1 = arma::diagmat(arma::ones(k)) + trans(Lmsg) * Lambda;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Veta1);
  eigval = arma::reverse(eigval);
  eigvec = arma::fliplr(eigvec);
  // ---------------------------------------------------------------------------
  // multiply an eigenvector by -1 if its first element is negative
  for (j = 0; j < k; j++) {
    eigvec.col(j) *= arma::sign(eigvec(0, j));
  }
  // ---------------------------------------------------------------------------
  arma::mat Tmat(k, k);
  if (arma::min(eigval) > 0.000001) {
    Tmat = trans(eigvec.each_row() % trans(sqrt(eigval)));
  } else {
    Tmat = arma::chol(Veta1);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q, R, Tmat);
  // ---------------------------------------------------------------------------
  // https://www.mathworks.com/matlabcentral/answers/83798-sign-differences-in-qr-decomposition
  // enforce positive diagonals of R
  R = arma::diagmat(arma::sign(R.diag())) * R;
  // ---------------------------------------------------------------------------
  arma::mat S = arma::inv(R);
  arma::mat M(n, k);
  for (j = 0; j < k; j++) {
    for (i = 0; i < n; i++) {
      M(i, j) = R::rnorm(0, 1);
    }
  }
  arma::mat eta = Z * Lmsg * S * trans(S) + M * trans(S);
  return eta;
}

/* -------------------------------------------------------------------------- */

//' Update the hth column of beta in the Adaptive Gibbs Sampler
//'
//' @param h An integer.
//' @param Xcov A pxq matrix.
//' @param Dt A pxk matrix.
//' @param Bh_1 A qxq matrix.
//' @param Phi_L A pxk matrix.
//' @param q An integer.
//'
//' @return A qx1 matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_beta(int h, arma::mat Xcov, arma::mat Dt,
                      arma::mat Bh_1, arma::mat Phi_L, int q) {
  arma::mat Qbeta = trans(Xcov.each_col() % Dt.col(h)) * Xcov + Bh_1;
  arma::mat Lbeta = trans(arma::chol(Qbeta));
  arma::vec bbeta = trans(Xcov) * (Phi_L.col(h) - 0.5);
  // mean
  arma::mat vbeta = solve(trimatl(Lbeta), bbeta);
  arma::mat mbeta = solve(trimatu(trans(Lbeta)), vbeta);
  // var
  arma::vec zbeta(q);
  for (int i = 0; i < q; i++) {
    zbeta(i) = R::rnorm(0, 1);
  }
  arma::mat ybeta = solve(trimatu(trans(Lbeta)), zbeta);
  arma::mat betah = trans(ybeta + mbeta);
  return betah;
}

/* -------------------------------------------------------------------------- */

//' Update the jth row of Lambda_star in the Adaptive Gibbs Sampler
//'
//' @param j An integer.
//' @param etarho A kxn matrix.
//' @param Phi A pxk matrix.
//' @param Plam A kxk matrix.
//' @param ps A p-dimensional vector;
//' @param Z A nxp matrix.
//' @param k An integer.
//'
//' @return A 1xk matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_Lambda_star(int j, arma::mat etarho, arma::mat Phi,
                             arma::mat Plam, arma::vec ps, arma::mat Z, int k) {
  arma::mat etaj = trans(etarho.each_col() % trans(Phi.row(j)));
  arma::mat Qlam = Plam + ps(j) * trans(etaj) * etaj;
  arma::mat Llam = trans(arma::chol(Qlam));
  arma::vec blam = ps(j) * trans(etaj) * Z.col(j);
  // mean
  arma::mat vlam = solve(trimatl(Llam), blam);
  arma::mat mlam = solve(trimatu(trans(Llam)), vlam);
  // var
  arma::vec zlam(k);
  for (int i = 0; i < k; i++) {
    zlam(i) = R::rnorm(0, 1);
  }
  arma::mat ylam = solve(trimatu(trans(Llam)), zlam);
  arma::mat lambda_starj = trans(ylam + mlam);
  return lambda_starj;
}

/* -------------------------------------------------------------------------- */

//' Update the hth element of d in the Adaptive Gibbs Sampler
//'
//' @param h An integer.
//' @param Phi A pxk matrix.
//' @param p An integer.
//' @param n An integer.
//' @param rho A k-dimensional vector.
//' @param eta A nxk matrix.
//' @param lambdastar A pxk matrix.
//' @param Z A nxp matrix.
//' @param sdy A nxp matrix.
//' @param k An integer.
//' @param w A k-dimensional vector.
//'
//' @return An integer in 1, ..., k.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
int update_d(int h, arma::mat Phi, int p, int n, arma::vec rho,
             arma::mat eta, arma::mat lambdastar, arma::mat Z,
             arma::mat sdy, int k, arma::vec w) {
  int i, j, l;
  double lnorm0 = 0.0, lnorm1 = 0.0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      // initialize mean
      double muijh = -1 * sqrt(rho(h) * Phi(j, h)) * lambdastar(j, h) * eta(i, h);
      for (l = 0; l < k; l++) {
        muijh += sqrt(rho(l) * Phi(j, l)) * lambdastar(j, l) * eta(i, l);
      }
      lnorm0 += std::log(R::dnorm(Z(i, j), muijh, sdy(i, j), FALSE));
      // update mean
      muijh += Phi(j, h) * lambdastar(j, h) * eta(i, h);
      lnorm1 += std::log(R::dnorm(Z(i, j), muijh, sdy(i, j), FALSE));
    }
  }
  // adjust the scale
  double mlnorm = std::max(lnorm0, lnorm1);
  lnorm0 -= mlnorm;
  lnorm1 -= mlnorm;
  Rcpp::NumericVector prob_h = as<NumericVector>(wrap(log(w)));
  for (i = 0; i <= h; i++) {
    prob_h(i) += lnorm0;
  }
  for (i = h + 1; i < k; i++) {
    prob_h(i) += lnorm1;
  }
  prob_h = exp(prob_h);
  prob_h = prob_h / sum(prob_h);
  if (sum(prob_h) == 0) {
    prob_h = rep(0, k);
    prob_h[k] = 1;
  }
  // draw d from a multinomial distribution
  Rcpp::IntegerVector d(k);
  R::rmultinom(1, prob_h.begin(), k, d.begin());
  return Rcpp::which_max(d);
}

/* -------------------------------------------------------------------------- */

//' Update Phi in the Adaptive Gibbs Sampler
//'
//' @param rho A k-dimensional vector.
//' @param logit A pxk matrix.
//' @param p_constant A number in (0,1).
//' @param p An integer.
//' @param n An integer.
//' @param eta A nxk matrix.
//' @param lambdastar A pxk matrix.
//' @param Phi A pxk matrix.
//' @param Z A nxp matrix.
//' @param sdy A nxp matrix.
//' @param k An integer.
//'
//' @return A pxk matrix.
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat update_Phi(arma::vec rho, arma::mat logit, double p_constant,
                     int p, int n, arma::mat eta, arma::mat lambdastar,
                     arma::mat Phi, arma::mat Z, arma::mat sdy, int k) {
  // define variables
  int i, j, l, h;
  arma::uword f;
  arma::vec p_phi0, p_phi1, p_phi_sum, lnorm0, lnorm1;
  arma::uvec wr0 = find(rho == 0);
  arma::uvec wr1 = find(rho == 1);
  // update for inactive factors
  for (f = 0; f < wr0.n_elem; f++) {
    h = wr0(f);
    p_phi0 = 1 - logit.col(h) * p_constant;
    p_phi1 = logit.col(h) * p_constant;
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
  }
  // update for active factors
  int ih = 0;
  for (f = 0; f < wr1.n_elem; f++) {
    h = wr1(f);
    lnorm0 = arma::zeros(p);
    lnorm1 = arma::zeros(p);
    for (j = 0; j < p; j++) {
      for (i = 0; i < n; i++) {
        // initialize mean
        double muijh = -1 * sqrt(rho(ih) * Phi(j, ih)) * lambdastar(j, ih) * eta(i, ih);
        for (l = 0; l < k; l++) {
            muijh += sqrt(rho(l) * Phi(j, l)) * lambdastar(j, l) * eta(i, l);
        }
        lnorm0(j) += std::log(R::dnorm(Z(i, j), muijh, sdy(i, j), FALSE));
        // update mean
        muijh += rho(h) * lambdastar(j, ih) * eta(i, ih);
        lnorm1(j) += std::log(R::dnorm(Z(i, j), muijh, sdy(i, j), FALSE));
      }
      // adjust the scale
      double mlnorm = std::max(lnorm0(j), lnorm1(j));
      lnorm0(j) -= mlnorm;
      lnorm1(j) -= mlnorm;
    }
    p_phi0 = exp(lnorm0 + log(1 - logit.col(h) * p_constant));
    p_phi1 = exp(lnorm1 + log(logit.col(h) * p_constant));
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
    ih += 1;
  }
  return Phi;
}

/* -------------------------------------------------------------------------- */
