#ifndef AGS_update_H
#define AGS_update_H

arma::mat update_bmu(const arma::mat& X, double prec_mu, double prec_b,
                     const arma::mat& mu, int q, int c);

arma::mat update_mu(int j, const arma::mat& Qbet, const arma::mat& W, const arma::mat& Z_res,
                    const arma::vec& ps, const arma::mat& b_mu, const arma::mat& Xcov, int c);

arma::mat update_eta(const arma::mat& Lambda, const arma::vec& ps, int k, const arma::mat& Z, int n);

arma::mat update_beta(int h, const arma::mat& Xcov, const arma::mat& Dt,
                      const arma::mat& Bh_1, const arma::mat& Phi_L, int q);

arma::mat update_Lambda_star(int j, const arma::mat& etarho, const arma::mat& Phi,
                             const arma::mat& Plam, const arma::vec& ps, const arma::mat& Z, int k);

int update_d(int h, const arma::mat& Phi, int p, int n, const arma::vec& rho,
             const arma::mat& eta, const arma::mat& lambdastar, const arma::mat& Z,
             const arma::mat& sdy, int k, const arma::vec& w);

arma::mat update_Phi(const arma::vec& rho, const arma::mat& logit, double p_constant,
                     int p, int n, const arma::mat& eta, const arma::mat& lambdastar,
                     arma::mat& Phi, const arma::mat& Z, const arma::mat& sdy, int k);

#endif
