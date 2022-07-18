// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// truncnorm_lg
arma::mat truncnorm_lg(arma::mat y_lower, arma::mat y_upper, arma::mat mu, arma::vec sigma, arma::mat u_rand);
RcppExport SEXP _GGIF_truncnorm_lg(SEXP y_lowerSEXP, SEXP y_upperSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP u_randSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_lower(y_lowerSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_upper(y_upperSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u_rand(u_randSEXP);
    rcpp_result_gen = Rcpp::wrap(truncnorm_lg(y_lower, y_upper, mu, sigma, u_rand));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GGIF_truncnorm_lg", (DL_FUNC) &_GGIF_truncnorm_lg, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_GGIF(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
