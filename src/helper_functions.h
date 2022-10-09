#ifndef helper_functions_H
#define helper_functions_H

arma::vec join_elem(arma::vec v1, int v2);
arma::vec join_elem(arma::vec v1, double v2);
arma::vec join_elem(int v1, arma::vec v2);
arma::vec join_elem(double v1, arma::vec v2);

arma::mat rbinom_vec(int len, int size, double prob);

arma::mat rnorm_mat(int rows, int cols, double mean, double sd);
arma::mat rnorm_vec(int len, double mean, double sd);

arma::mat runif_mat(int rows, int cols, double minVal, double maxVal);
arma::mat runif_vec(int len, double minVal, double maxVal);

#endif
