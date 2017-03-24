#ifndef LD_H
#define LD_H
#include <RcppEigen.h>
#include "daal.h"
#include "mkl.h"
#include "RSSReQTL_types.h"


double calc_nmsum(const double m);

Eigen::MatrixXd calc_dist(c_arrayxd_internal map);

Eigen::MatrixXd calc_cov( c_Matrix_internal mat);

Eigen::ArrayXd calc_variance(c_Matrix_internal mat);

void cov_2_cor(Matrix_internal covmat, arrayxd_internal rowvar);

Eigen::MatrixXd cov_2_cor_exp(Matrix_external covmat, arrayxd_external rowvar);

void compute_shrinkage(Matrix_internal distmat, c_Matrix_internal S,c_Matrix_internal hmata,const double theta, const double m,const double Ne,const double cutoff);

void compute_shrinkage_cor(Matrix_internal distmat, c_Matrix_internal S,c_Matrix_internal hmata,const double theta, const double m,const double Ne,const double cutoff);

Eigen::MatrixXd calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff);

Eigen::SparseMatrix<double> sp_calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff);
  
#endif