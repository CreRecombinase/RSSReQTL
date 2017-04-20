#include <RSSReQTL.h>
// [[Rcpp::depends(rssr)]]



Rcpp::DataFrame mt_rss_varbvsr(
    const c_Matrix_internal R,
    const c_arrayxd_internal sigma_beta,
    const c_arrayxd_internal logodds,
    const c_Matrix_internal  betahat_mat,
    const c_Matrix_internal  se_mat,
    const arrayxd_internal talpha0,
    const arrayxd_internal tmu0,
    double tolerance,
    int itermax,
    bool isVerbose,
    bool islnz_tol){
  
  using namespace Rcpp;
  size_t gridsize = sigma_beta.size();
  
  
  if(gridsize!=logodds.size()){
    Rcpp::stop("sigma_beta.size() !=logodds.size()");
  }
  size_t genes= betahat_mat.cols();
  size_t retsize=genes*gridsize;
  Rcpp::NumericVector lovec(retsize);
  Rcpp::NumericVector sigbvec(retsize);
  Rcpp::NumericVector nlzvec(retsize);
  Rcpp::IntegerVector fgeneid(retsize);
  size_t snps=betahat_mat.rows();
  for(int i=0;i<genes;i++){
    Eigen::ArrayXd se=se_mat.col(i);
    Eigen::ArrayXd betahat=betahat_mat.col(i);
    Eigen::MatrixXd srs= se.matrix().asDiagonal()*R*se.matrix().asDiagonal();
    for(int gridsize=0;i<gridsize;i++){
      
    }
  }
  
  
  
  return(Rcpp::DataFrame::create(_["logodds"]=lovec,
                                 _["sigb"]=sigbvec,
                                 _["lnZ"]=nlzvec,
                                 _["fgeneid"]=fgeneid));
}
  
  