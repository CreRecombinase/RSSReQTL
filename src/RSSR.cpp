#include <RSSReQTL.h>
#include <RcppParallel.h>
//[[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(rssr)]]



// 
// 
// 
// class src_body {
//   const c_Matrix_internal betahat_matrix;
//   const c_Matrix_internal se_matrix;
//   const int ngenes;
//   int genei;
// public:
//   src_body(const c_Matrix_internal bhm, const c_Matrix_internal sem, int ng) : betahat_matrix(bhm), genei(1),se_matrix(sem),ngenes(sem.cols()) {}
//   bool operator()( arrayxd_internal se,arrayxd_internal betahat) {
//     if ( genei <= ngenes ) {
//       se=se_matrix.col(genei-1);
//       betahat=betahat_matrix.col(genei-1);
//       genei++;
//       return true;
//     } else {
//       return false;
//     }
//   }
// };
//
// // int main() {
//   int sum = 0;
//   graph g;
//   function_node< int, int > squarer( g, unlimited, [](const int &v) {
//     return v*v;
//   } );
//   function_node< int, int > cuber( g, unlimited, [](const int &v) {
//     return v*v*v;
//   } );
//   function_node< int, int > summer( g, 1, [&](const int &v ) -> int {
//     return sum += v;
//   } );
//   make_edge( squarer, summer );
//   make_edge( cuber, summer );
//   source_node< int > src( g, src_body(10), false );
//   make_edge( src, squarer );
//   make_edge( src, cuber );
//   src.activate();
//   g.wait_for_all();
//   cout << "Sum is " << sum << "\n";
// // }
//
//
//
//
//
// using namespace tbb;
//
// Rcpp::DataFrame mt_rss_varbvsr(
//     const c_Matrix_internal R,
//     const c_arrayxd_internal sigma_beta,
//     const c_arrayxd_internal logodds,
//     const c_Matrix_internal  betahat_mat,
//     const c_Matrix_internal  se_mat,
//     const arrayxd_internal talpha0,
//     const arrayxd_internal tmu0,
//     double tolerance,
//     int itermax,
//     bool isVerbose,
//     bool islnz_tol){
//
//   using namespace Rcpp;
//   size_t gridsize = sigma_beta.size();
//
//
//   if(gridsize!=logodds.size()){
//     Rcpp::stop("sigma_beta.size() !=logodds.size()");
//   }
//   size_t genes= betahat_mat.cols();
//   size_t retsize=genes*gridsize;
//   Rcpp::NumericVector lovec(retsize);
//   Rcpp::NumericVector sigbvec(retsize);
//   Rcpp::NumericVector nlzvec(retsize);
//   Rcpp::IntegerVector fgeneid(retsize);
//   size_t snps=betahat_mat.rows();
//
//
//   using namespace Rcpp;
//   size_t sigb_size= sigma_beta.size();
//   size_t logodds_size=logodds.size();
//   size_t tot_size=sigb_size*logodds_size;
//
//   Rcpp::NumericVector nlzvec(tot_size);
//   Rcpp::NumericVector sigbvec(tot_size);
//   Rcpp::NumericVector lovec(tot_size);
//
//   Rcpp::LogicalVector verbose(1);
//   verbose(0)=isVerbose;
//   Rcpp::LogicalVector lnz_tol(1);
//   lnz_tol(0)=islnz_tol;
//
//   parallel_for(blocked_range<size_t>(0,genes),
//                [&](const blocked_range<size_t>& r){
//                  for(size_t t=r.begin(); t!=r.end(); t++){
//
//                    Eigen::ArrayXd se=se_mat.col(i);
//                    Eigen::ArrayXd betahat=betahat_mat.col(i);
//                    Eigen::MatrixXd srs= se.matrix().asDiagonal()*R*se.matrix().asDiagonal();
//                    for(int gridsize=0;i<gridsize;i++){
//                    size_t i=t%logodds_size;
//                    size_t j=t/logodds_size;
//                    sigbvec(t)=sigma_beta(j);
//                    lovec(t)=logodds(i);
//                    nlzvec(t)=rss_varbvsr_squarem_iter(SiRiS,
//                           sigma_beta(j),
//                           logodds(i),
//                           betahat,
//                           se,
//                           talpha0,
//                           tmu0,
//                           tSiRiSr0,
//                           tolerance,
//                           itermax,
//                           lnz_tol);}});
//
//
//
//
//   for(int i=0;i<genes;i++){
//
//
//     }
//   }
//
//
//
//   return(Rcpp::DataFrame::create(_["logodds"]=lovec,
//                                  _["sigb"]=sigbvec,
//                                  _["lnZ"]=nlzvec,
//                                  _["fgeneid"]=fgeneid));
// }
//   
//   