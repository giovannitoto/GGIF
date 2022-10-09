#ifndef AGS_update_H
#define AGS_update_H

arma::mat update_bmu(arma::mat X, double prec_mu, double prec_b,
                     arma::mat mu, int q, int c);

arma::mat update_mu(int j, arma::mat Qbet, arma::mat W, arma::mat Z_res,
                    arma::vec ps, arma::mat b_mu, arma::mat Xcov, int c);

arma::mat update_eta(arma::mat Lambda, arma::vec ps, int k, arma::mat Z, int n);

arma::mat update_beta(int h, arma::mat Xcov, arma::mat Dt,
                      arma::mat Bh_1, arma::mat Phi_L, int q);

arma::mat update_Lambda_star(int j, arma::mat etarho, arma::mat Phi,
                             arma::mat Plam, arma::vec ps, arma::mat Z, int k);

int update_d(int h, arma::mat Phi, int p, int n, arma::vec rho,
             arma::mat eta, arma::mat lambdastar, arma::mat Z,
             arma::mat sdy, int k, arma::vec w);

arma::mat update_Phi(arma::vec rho, arma::mat logit, double p_constant,
                     int p, int n, arma::mat eta, arma::mat lambdastar,
                     arma::mat Phi, arma::mat Z, arma::mat sdy, int k);

#endif
