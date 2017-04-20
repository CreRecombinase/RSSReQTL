// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include <RSSReQTL.h>
// [[Rcpp::depends(RcppEigen)]]


Eigen::MatrixXd eqtl_orthogonalize_covar( c_Matrix_internal covariates){
  Eigen::HouseholderQR<Eigen::MatrixXd> qr(covariates);
  Eigen::MatrixXd Q=qr.householderQ();
  return(Q.leftCols(covariates.cols()));
}

void eqtl_orthogonalize_data(Matrix_internal data, c_Matrix_internal ortho_covar){
  data=data-((data.transpose()*ortho_covar)*ortho_covar.transpose()).transpose();
}


//[[Rcpp::export(name="orthogonalize_data")]]
Eigen::MatrixXd orthogonalize_data_exp(Matrix_external data, Matrix_external ortho_covar){
  Eigen::MatrixXd tdata=data;
  eqtl_orthogonalize_data(tdata,ortho_covar);
  return(tdata);
}

//[[Rcpp::export]]
Eigen::MatrixXd rssr_orthogonalize_covar( Matrix_external covariates){
  return(eqtl_orthogonalize_covar(covariates));
}


  
void fast_eqtl_lm(c_arrayxd_internal genotype, c_arrayxd_internal expression, const int n,arrayxd_internal tresid, arrayxd_internal retvec ){
  
  double xtx = genotype.matrix().dot(genotype.matrix());
  retvec(0)=expression.matrix().dot(genotype.matrix())/xtx;
  tresid=expression-genotype*retvec(0);
  double s2=(tresid.matrix().dot(tresid.matrix()))/(n-1);
  retvec(1)=sqrt(s2/xtx);
  
}

Eigen::MatrixXd map_eqtl_lm(Matrix_internal genotype, arrayxd_internal expression){
  int p=genotype.cols();
  int n=expression.size();
  Eigen::MatrixXd retmat(2,p);
  Eigen::ArrayXd retvec(2);
  Eigen::ArrayXd tresid(n);
  for(int i=0; i<p;i++){
    fast_eqtl_lm(genotype.col(i),expression,n,tresid,retmat.col(i));
  }
  retmat.transposeInPlace();
  return(retmat);
}

Eigen::MatrixXd map_beta(const Matrix_internal genotype,const Matrix_internal expression){
  
  Eigen::ArrayXd sx2 =genotype.array().square().colwise().sum();
  //  Rcpp::Rcout<<"sx2 has: "<<sx2.size()<<"elements"<<std::endl;
  return((genotype.transpose()*expression).array().colwise()/sx2);
  
}

Eigen::MatrixXd map_se(const Matrix_internal genotype,const Matrix_internal expression, const Matrix_internal betahat){
  size_t s=genotype.cols();
  size_t g=expression.cols();
  size_t n=genotype.rows();
  Eigen::MatrixXd semat(s,g);
  Eigen::ArrayXd sx2 =genotype.array().square().colwise().sum();
  Eigen::MatrixXd resmat(n,s);
  Eigen::ArrayXd resvec(s);
  Eigen::MatrixXd yh(n,s);
  for(int i=0;i<g;i++){
    yh = genotype.array().rowwise()*betahat.col(i).array().transpose();
    resmat =((yh.array().colwise()-expression.col(i).array()));
    resvec = (resmat.array().square()/n).colwise().sum().sqrt();
    semat.col(i)=resvec.array()*sx2.sqrt().inverse();
  }
  return(semat);
}

//[[Rcpp::export]]
Eigen::MatrixXd map_se_exp(const Matrix_external genotype,const Matrix_external expression,const Matrix_external betahat){
  return(map_se(genotype,expression,betahat));
}


//[[Rcpp::export]]
Eigen::MatrixXd map_beta_exp(const Matrix_external genotype,const Matrix_external expression){
  
  return(map_beta(genotype,expression));
  
}


//[[Rcpp::export(name="map_eqtl_lm")]]
Rcpp::DataFrame map_eqtl_lm_exp(Matrix_external genotype, arrayxd_external expression){
  using namespace Rcpp;
  Eigen::MatrixXd retmat=map_eqtl_lm(genotype,expression);
  return(Rcpp::DataFrame::create(_["betahat"]=Rcpp::wrap(retmat.col(0)),
                          _["se"]=Rcpp::wrap(retmat.col(1))));

}

