// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/RSSReQTL.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_nmsum
double calc_nmsum(const double m);
RcppExport SEXP RSSReQTL_calc_nmsum(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_nmsum(m));
    return rcpp_result_gen;
END_RCPP
}
// calc_dist_exp
Eigen::MatrixXd calc_dist_exp(arrayxd_external map);
RcppExport SEXP RSSReQTL_calc_dist_exp(SEXP mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arrayxd_external >::type map(mapSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_dist_exp(map));
    return rcpp_result_gen;
END_RCPP
}
// calc_cov_mkl_exp
Eigen::MatrixXd calc_cov_mkl_exp(Matrix_external mat);
RcppExport SEXP RSSReQTL_calc_cov_mkl_exp(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_cov_mkl_exp(mat));
    return rcpp_result_gen;
END_RCPP
}
// calc_cov_exp
Eigen::MatrixXd calc_cov_exp(Matrix_external mat);
RcppExport SEXP RSSReQTL_calc_cov_exp(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_cov_exp(mat));
    return rcpp_result_gen;
END_RCPP
}
// calc_variance_exp
Eigen::ArrayXd calc_variance_exp(Matrix_external mat);
RcppExport SEXP RSSReQTL_calc_variance_exp(SEXP matSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type mat(matSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_variance_exp(mat));
    return rcpp_result_gen;
END_RCPP
}
// cov_2_cor_exp
Eigen::MatrixXd cov_2_cor_exp(Matrix_external covmat, arrayxd_external rowvar);
RcppExport SEXP RSSReQTL_cov_2_cor_exp(SEXP covmatSEXP, SEXP rowvarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type covmat(covmatSEXP);
    Rcpp::traits::input_parameter< arrayxd_external >::type rowvar(rowvarSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_2_cor_exp(covmat, rowvar));
    return rcpp_result_gen;
END_RCPP
}
// calc_shrinkage
Eigen::MatrixXd calc_shrinkage(Matrix_external distmat, Matrix_external S, Matrix_external hmata, const double theta, const double m, const double Ne, const double cutoff);
RcppExport SEXP RSSReQTL_calc_shrinkage(SEXP distmatSEXP, SEXP SSEXP, SEXP hmataSEXP, SEXP thetaSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< Matrix_external >::type S(SSEXP);
    Rcpp::traits::input_parameter< Matrix_external >::type hmata(hmataSEXP);
    Rcpp::traits::input_parameter< const double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_shrinkage(distmat, S, hmata, theta, m, Ne, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// calcLD_exp
Eigen::MatrixXd calcLD_exp(Matrix_external hmata, arrayxd_external mapa, const double m, const double Ne, const double cutoff);
RcppExport SEXP RSSReQTL_calcLD_exp(SEXP hmataSEXP, SEXP mapaSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type hmata(hmataSEXP);
    Rcpp::traits::input_parameter< arrayxd_external >::type mapa(mapaSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(calcLD_exp(hmata, mapa, m, Ne, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// sp_calcLD_exp
Eigen::SparseMatrix<double> sp_calcLD_exp(Matrix_external hmata, arrayxd_external mapa, const double m, const double Ne, const double cutoff);
RcppExport SEXP RSSReQTL_sp_calcLD_exp(SEXP hmataSEXP, SEXP mapaSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type hmata(hmataSEXP);
    Rcpp::traits::input_parameter< arrayxd_external >::type mapa(mapaSEXP);
    Rcpp::traits::input_parameter< const double >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(sp_calcLD_exp(hmata, mapa, m, Ne, cutoff));
    return rcpp_result_gen;
END_RCPP
}
// orthogonalize_data_exp
Eigen::MatrixXd orthogonalize_data_exp(Matrix_external data, Matrix_external ortho_covar);
RcppExport SEXP RSSReQTL_orthogonalize_data_exp(SEXP dataSEXP, SEXP ortho_covarSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Matrix_external >::type ortho_covar(ortho_covarSEXP);
    rcpp_result_gen = Rcpp::wrap(orthogonalize_data_exp(data, ortho_covar));
    return rcpp_result_gen;
END_RCPP
}
// rssr_orthogonalize_covar
Eigen::MatrixXd rssr_orthogonalize_covar(Matrix_external covariates);
RcppExport SEXP RSSReQTL_rssr_orthogonalize_covar(SEXP covariatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type covariates(covariatesSEXP);
    rcpp_result_gen = Rcpp::wrap(rssr_orthogonalize_covar(covariates));
    return rcpp_result_gen;
END_RCPP
}
// map_eqtl_lm_exp
Rcpp::DataFrame map_eqtl_lm_exp(Matrix_external genotype, arrayxd_external expression);
RcppExport SEXP RSSReQTL_map_eqtl_lm_exp(SEXP genotypeSEXP, SEXP expressionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Matrix_external >::type genotype(genotypeSEXP);
    Rcpp::traits::input_parameter< arrayxd_external >::type expression(expressionSEXP);
    rcpp_result_gen = Rcpp::wrap(map_eqtl_lm_exp(genotype, expression));
    return rcpp_result_gen;
END_RCPP
}
