#include <RSSReQTL.h>

// [[Rcpp::depends(RcppEigen)]]




//[[Rcpp::export]]
double calc_nmsum(const double m){
  int msize=(2*(int)m-1);
  Eigen::ArrayXd tx(msize);
  tx.setLinSpaced(msize,1,(int)(2*m-1));
  return  (1/tx).sum();
}

//[[Rcpp::export]]
Rcpp::DataFrame ld2df_sp(const  sparseMatrix_external ldmat, Rcpp::StringVector rsid,const double r2cutoff=0.01){
  
  using namespace Eigen;
  
  using namespace Rcpp;
  size_t p=ldmat.rows();
  
  if(p!=ldmat.cols()){
    Rcpp::stop("ldmat is not square!");
  }
  if(p!=rsid.size()){
    Rcpp::stop("rsid must be of length p!");
  }
  size_t dfsize = ldmat.nonZeros()-p;
  
  
  
  std::vector<double>corv;
  corv.reserve(dfsize);
  std::vector<std::string> rowsnp;
  std::vector<std::string> colsnp;
  rowsnp.reserve(dfsize);
  colsnp.reserve(dfsize);
  size_t ssize = rsid[0].size();
  // Rcpp::Rcout<<"Generating DataFrame"<<std::endl;
  for(int k=0;k<ldmat.outerSize();++k){
    for(SparseMatrix<double>::InnerIterator it(ldmat,k);it;++it)
    {
      if(it.row()!=it.col()){
        double r2=it.value()*it.value();
        if(r2>r2cutoff){
          corv.push_back(r2);
          rowsnp.push_back(as<std::string>(rsid[it.row()]));
          colsnp.push_back(as<std::string>(rsid[it.col()]));
          
        }
        
      }
    }
  }
  // Rcpp::Rcout<<"Returning DataFrame"<<std::endl;
  
  return(DataFrame::create(_["rowsnp"]=wrap(rowsnp),_["colsnp"]=wrap(colsnp),_["r2"]=wrap(corv),_["stringsAsFactors"]=false));
}

//[[Rcpp::export]]
Rcpp::DataFrame ld2df(const Matrix_external ldmat, Rcpp::StringVector rsid,const double r2cutoff=0.01){
  using namespace Rcpp;
  size_t p=ldmat.rows();
  if(p!=ldmat.cols()){
    Rcpp::stop("ldmat is not square!");
  }
  if(p!=rsid.size()){
    Rcpp::stop("rsid must be of length p!");
  }
  size_t dfsize = (p*(p-1))/2;
  std::vector<double>corv;
  corv.reserve(dfsize);
  std::vector<std::string> rowsnp;
  std::vector<std::string> colsnp;
  rowsnp.reserve(dfsize);
  colsnp.reserve(dfsize);
  
  
  // Eigen::ArrayXd corv(dfsize);
  // Rcpp::StringVector rowsnp(dfsize);
  // Rcpp::StringVector colsnp(dfsize);
  // Rcpp::Rcout<<"Generating DataFrame"<<std::endl;
  size_t k=0;
  for(int i=0; i<p;i++){
    for(int j=i+1; j<p;j++ ){
      double r2=ldmat.coeff(i,j)*ldmat.coeff(i,j);
      if(r2>r2cutoff){
        corv.push_back(r2);
        rowsnp.push_back(as<std::string>(rsid[i]));
        colsnp.push_back(as<std::string>(rsid[j]));
      }
    }
  }
  // Rcpp::Rcout<<"Returning DataFrame"<<std::endl;
  
  return(DataFrame::create(_["rowsnp"]=wrap(rowsnp),_["colsnp"]=wrap(colsnp),_["r2"]=wrap(corv),_["stringsAsFactors"]=false));
}
  



Eigen::MatrixXd calc_dist(c_arrayxd_internal map){
  int p=map.size();
  Eigen::MatrixXd retmatrix(p,p);
  retmatrix.setZero();
  double tj=0;
  double ti=0;
  for(int i=0; i<p; i++){
    ti = map.coeff(i);
    for(int j=i+1;j<p; j++){
      tj=map.coeff(j);
      if(tj<0){
        Rcpp::stop("Negative map value encountered!");
      }
      retmatrix.coeffRef(i,j)=tj-ti;
    }
  }
 
  return(retmatrix);
}
//[[Rcpp::export(name="calc_dist")]]
Eigen::MatrixXd calc_dist_exp( arrayxd_external map){
  return(calc_dist(map));
}




Eigen::MatrixXd calc_cov( c_Matrix_internal mat){
  auto centered = mat.rowwise()-mat.colwise().mean();
  return (((centered.adjoint()*centered)/double(mat.rows()-1)));  
}



//[[Rcpp::export]]
Eigen::MatrixXd test_calc_shrink(arrayxd_external map,Matrix_external Hpanel,  const double Ne, const double m,const double cutoff){
  int p=map.size();
  Eigen::MatrixXd S=calc_cov(Hpanel);
  Eigen::MatrixXd retmatrix(p,p);
  retmatrix.setZero();
  double tj=0;
  double ti=0;
  for(int i=0; i<p; i++){
    ti = map.coeff(i);
    for(int j=i+1;j<p; j++){
      tj=map.coeff(j);
      double rho = exp(-(4*Ne*(tj-ti)/100)/(2*m));
      if(rho>=cutoff){
        retmatrix.coeffRef(i,j)=rho;
      }
     
    }
  }
  retmatrix.diagonal().setOnes(); 
  
  S=S.cwiseProduct(retmatrix);
  return(S);
}


// 
// Eigen::MatrixXd calc_cov_daal( Matrix_internal mat){
//   using namespace daal;
//   using namespace daal::data_management;
//   using namespace daal::algorithms;
//   using namespace daal::services;
//   
//   int *dim;
//   int nFeatures, nObservations;
//   double *mean;
// 
//   /* Convert input arguments to C types */
// 
//   double *x = mat.data();
//   
//   /* Get input matrix size */
//   
//   nObservations = mat.rows();
//   nFeatures     = mat.cols();
//   
//   /* Create structure-of-arrays (SOA) numeric table to store input data */
//   SharedPtr<SOANumericTable> dataTable(new SOANumericTable(nFeatures, nObservations));
//   for (int i = 0; i < nFeatures; i++)
//   {
//     dataTable->setArray(x + i*nObservations, i);
//   }
//   
//   /* Allocate memory to store results */
//   Eigen::MatrixXd cov(nFeatures,nFeatures);
//   
//   mean = new double[nFeatures];
//   
//   /* Create homogeneous numeric tables to store results */
//   SharedPtr<HomogenNumericTable<> > covarianceTable (new HomogenNumericTable<>(cov.data(), nFeatures, nFeatures));
//   SharedPtr<HomogenNumericTable<> > meanTable       (new HomogenNumericTable<>(mean, nFeatures, 1));
//   
//   /* Create algorithm to compute covariance matrix using default method */
//   covariance::Batch<> algorithm;
//   algorithm.input.set(covariance::data, dataTable);
//   
//   /* Create object to store the results of DAAL computations */
//   SharedPtr<covariance::Result> result(new covariance::Result());
//   
//   /* Provide memory for storing the results of DAAL computations */
//   result->set(covariance::covariance, covarianceTable);
//   result->set(covariance::mean,       meanTable);
//   
//   /* Register the object for storing results in DAAL algorithm */
//   algorithm.setResult(result);
//   
//   /* Compute covariance matrix */
//   algorithm.compute();
//   
//   delete [] mean;
//   
//   /* Return covariance matrix */
//   return cov;
//   
//   
// }

// 
// Eigen::MatrixXd calc_cov_mkl(Matrix_internal mat){
// 
//   VSLSSTaskPtr task;
// 
//   
// 
// 
//   int rows=mat.rows();
//   int cols=mat.cols();
//   Eigen::ArrayXd meana(cols);
//   Eigen::ArrayXd vara(cols);
//   
//   double *x=mat.data();
//   double *mean=meana.data(); 
//   double *var= vara.data(); 
//   double *w=0;
//   Eigen::MatrixXd cov(cols,cols);
//   MKL_INT p,n,xstorage,covstorage;
//   int status;
//   
//   p=cols;
//   n=rows;
//   xstorage=VSL_SS_MATRIX_STORAGE_COLS;
//   covstorage=VSL_SS_MATRIX_STORAGE_FULL;
//   
//   int errcode = vsldSSNewTask(&task,&p,&n,&xstorage,x,w,0);
//     
//   errcode = vsldSSEditCovCor(task,meana.data(),cov.data(),&covstorage,0,0);
//   int estimates = VSL_SS_COV;
//   status = vsldSSCompute(task,estimates,VSL_SS_METHOD_1PASS);
//   status = vslSSDeleteTask(&task);
//   return(cov);
// }


// Eigen::MatrixXd calc_cov_mkl_exp(Matrix_external mat){
//   return(calc_cov_mkl(mat));
// }

// 
// 
// Eigen::MatrixXd calc_cov_daal_exp( Matrix_external mat){
//   return(calc_cov_daal(mat));
// }

//[[Rcpp::export(name="calc_cov")]]
Eigen::MatrixXd calc_cov_exp(Matrix_external mat){
  return(calc_cov(mat));
}


Eigen::ArrayXd calc_variance(c_Matrix_internal mat){
  int n=mat.rows();
  return(((mat.rowwise()-(mat.colwise().mean())).array().square().colwise().sum())/(n-1));
}

//[[Rcpp::export(name="calc_variance")]]
Eigen::ArrayXd calc_variance_exp(Matrix_external mat){
  return(calc_variance(mat));
}


void cov_2_cor(Matrix_internal covmat, arrayxd_internal rowvar){
  rowvar = 1/rowvar.sqrt();
  covmat.array().colwise()*=rowvar;
  covmat.array().rowwise()*=rowvar.transpose();
}

//[[Rcpp::export(name="cov_2_cor")]]
Eigen::MatrixXd cov_2_cor_exp(Matrix_external covmat, arrayxd_external rowvar){
  
  Eigen::ArrayXd trowvar=rowvar;
  Eigen::MatrixXd tcovmat=covmat;
  cov_2_cor(tcovmat,trowvar);
  return(tcovmat);
  
}





void compute_shrinkage(Matrix_internal distmat, c_Matrix_internal S,c_Matrix_internal hmata,const double theta, const double m,const double Ne,const double cutoff){
  // auto tdv=distmat.selfadjointView<Eigen::Upper>();
  // tdv=(tdv*(4*Ne))/100;
  
  int n=hmata.rows();
  distmat=4*Ne*distmat/100;
  distmat =(-distmat/(2*m)).array().exp();
  int tot_size=distmat.rows()*distmat.cols();
  
  for(int i=0;i<tot_size;i++){
    double* temp=distmat.data()+i;
    if((*temp)<cutoff){
      temp=0;
    }
  }
  distmat.diagonal().setOnes();
  distmat=distmat.cwiseProduct(S);
   Eigen::ArrayXd vars=calc_variance(hmata);
  
  distmat.diagonal()=vars;
  Eigen::VectorXd diagvec(distmat.rows());
  diagvec.setConstant(0.5*theta*(1-0.5*theta));
  distmat = (1-theta)*(1-theta)*distmat;
  distmat.diagonal()=distmat.diagonal()+diagvec;
  vars=vars*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  
}


void compute_shrinkage_cor(Matrix_internal distmat, c_Matrix_internal S,c_Matrix_internal hmata,const double theta, const double m,const double Ne,const double cutoff){

  
  
  //Editing so if we detect genotype instead of haplotype, we'll simply multiply S by 0.5
  double dosage_max=hmata.maxCoeff();
  bool isGeno= dosage_max>1;
//  int n=hmata.rows();
  distmat=4*Ne*distmat/100;
  distmat =(-distmat/(2*m)).array().exp();
  int tot_size=distmat.rows()*distmat.cols();
  
  for(int i=0;i<tot_size;i++){
    double* temp=distmat.data()+i;
    if((*temp)<cutoff){
      temp=0;
    }
  }
  distmat=distmat.cwiseProduct(S);
  Eigen::ArrayXd vars=calc_variance(hmata);
   if(isGeno){
     // distmat=distmat*0.5;
     vars=vars*0.5;
   }
  
  
  distmat.diagonal()=vars;
  Eigen::VectorXd diagvec(distmat.rows());
  diagvec.setConstant(0.5*theta*(1-0.5*theta));
  distmat = (1-theta)*(1-theta)*distmat;
  distmat.diagonal()=distmat.diagonal()+diagvec;
  vars=vars*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);

  cov_2_cor(distmat,vars);
}



//[[Rcpp::export]]
Eigen::MatrixXd calc_shrinkage(Matrix_external distmat, Matrix_external S, Matrix_external hmata,const double theta, const double m,const double Ne,const double cutoff){
  
  Eigen::MatrixXd tdistmat=distmat;
  compute_shrinkage(tdistmat,S,hmata,theta,m,Ne,cutoff);
  return(tdistmat);
}


// Eigen::MatrixXd calc_LD_ns(c_Matrix_internal hmata){
//   Eigen::ArrayXd vars = calc_variance(hmata);
//   Eigen::ArrayXd means= hmata.colwise().mean();
// 
//   // return(((mat.rowwise()-(mat.colwise().mean())).array().square().colwise().sum())/(n-1));
//   //
//   // return cov_2_cor(calc_cov(hmata));
// }



Eigen::MatrixXd calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){

  double nmsum=calc_nmsum(m);
  double theta=(1/nmsum)/(2*m+(1/nmsum));
  Eigen::MatrixXd dist=calc_dist(mapa);
  Eigen::MatrixXd S=calc_cov(hmata);
  
  double dosage_max=hmata.maxCoeff();
  bool isGeno= dosage_max>1;
  if(isGeno){
    S=S*0.5;
  }
  
  compute_shrinkage_cor(dist,S,hmata,theta,m,Ne,cutoff);
  //  Eigen::SparseMatrix<double> retmat=dist.sparseView();
  return(dist);
}


//[[Rcpp::export(name="calcLD")]]
Eigen::MatrixXd calcLD_exp(Matrix_external hmata,arrayxd_external mapa,const double m, const double Ne, const double cutoff){
  return(calcLD(hmata,mapa,m,Ne,cutoff));
}



Eigen::SparseMatrix<double> sp_calcLD(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){
  
  // compute_shrinkage_cor(dist,S,hmata,theta,m,Ne,cutoff);
  
  Eigen::MatrixXd dist = calcLD(hmata,mapa,m,Ne,cutoff);
  // dist.triangularView<Eigen::Lower>()=dist.transpose();
  Eigen::SparseMatrix<double> retmat=dist.sparseView();
  return(retmat);
}

//[[Rcpp::export(name="sp_calcLD")]]
Eigen::SparseMatrix<double> sp_calcLD_exp(Matrix_external hmata,arrayxd_external mapa,const double m, const double Ne, const double cutoff){
  return(sp_calcLD(hmata,mapa,m,Ne,cutoff));
}



Eigen::SparseMatrix<double> sp_calcLD_symm(c_Matrix_internal hmata,c_arrayxd_internal mapa,const double m, const double Ne, const double cutoff){
  
  // compute_shrinkage_cor(dist,S,hmata,theta,m,Ne,cutoff);
  
  Eigen::MatrixXd dist = calcLD(hmata,mapa,m,Ne,cutoff);
  // dist.triangularView<Eigen::Lower>()=dist.transpose();
  Eigen::SparseMatrix<double> retmat=Eigen::MatrixXd(dist.triangularView<Eigen::Upper>()).sparseView();
  return(retmat);
}

//[[Rcpp::export(name="sp_calcLD_symm")]]
Eigen::SparseMatrix<double> sp_calcLD_symm_exp(Matrix_external hmata,arrayxd_external mapa,const double m, const double Ne, const double cutoff){
  return(sp_calcLD_symm(hmata,mapa,m,Ne,cutoff));
}
