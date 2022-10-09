/* -------------------------------------------------------------------------- */

#include <RcppArmadillo.h>
using namespace Rcpp;

/* -------------------------------------------------------------------------- */

arma::vec join_elem(arma::vec v1, int v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}
arma::vec join_elem(arma::vec v1, double v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}
arma::vec join_elem(int v1, arma::vec v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}
arma::vec join_elem(double v1, arma::vec v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}

/* -------------------------------------------------------------------------- */

arma::mat rbinom_vec(int len, int size, double prob) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rbinom(size, prob);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat rnorm_mat(int rows, int cols, double mean, double sd) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::rnorm(mean, sd);
    }
  }
  return M;
}
arma::mat rnorm_vec(int len, double mean, double sd) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rnorm(mean, sd);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat runif_mat(int rows, int cols, double minVal, double maxVal) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::runif(minVal, maxVal);
    }
  }
  return M;
}
arma::mat runif_vec(int len, double minVal, double maxVal) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::runif(minVal, maxVal);
  }
  return v;
}

/* -------------------------------------------------------------------------- */
