

varbvs_mean <- function(X,ymat,sigb,log10odds,tolerance=1e-3,itermax=100,fgeneid=NULL,tlog10odds=NULL,tsigb=NULL,tresid=NULL,verbose=F){
  library(dplyr)
  library(progress)
  library(varbvs)
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(ymat)
  }
  
  ng <- ncol(ymat)
  retl <- list()
  p <- ncol(X)
  n <- nrow(X)

  pb <- progress_bar$new(total=ng)
  if(!is.null(tsigb)){
    stopifnot(length(tsigb)==ng)

  }
  if(!is.null(tlog10odds)){
    stopifnot(length(tlog10odds)==ng)

  }
  if(!is.null(tresid)){
    stopifnot(length(tresid)==ng)
  }
  
  mean_logodds <- numeric(ng)
  mean_sigb <- numeric(ng)
  mean_alpha <- numeric(ng)
  mean_pi <- numeric(ng)

  
  
  for(i in 1:ng){
    pb$tick()
    if(!is.null(tlog10odds)){
      glog10odds <- tlog10odds[i]
    }else{
      glog10odds <- log10odds
    }
    if(!is.null(tsigb)){
      gsigb <- tsigb[i]
    }else{
      gsigb <- sigb^2
    }
    if(!is.null(tresid)){
      gsigma <- tresid[i]^2
      gsigb <- gsigb/gsigma
    }else{
      gsigma <- NULL
    }
    if(!is.null(gsigma)){
      tvb <- varbvs(X,Z = NULL,y = ymat[,i],sigma = gsigma,
                    sa = gsigb,
                    logodds = glog10odds,
                    update.sa = F,
                    update.sigma = F,
                    maxiter = itermax,tol=tolerance,verbose=verbose)
      
    }else{
      tvb <- varbvs(X,Z = NULL,y = ymat[,i],
                    sa = gsigb,
                    logodds = glog10odds,
                    update.sa = F,
                    update.sigma = F,
                    maxiter = itermax,tol=tolerance,verbose=verbose)
    }

    # tvb<- varbvs::varbvs(X = X,Z = NULL,y = ymat[,i],sa = sigb,logodds = log10odds,update.sigma = F,verbose=F,tol = 1e-3)
    tvb_w <- normalizelogweights(tvb$logw)
    tvb_mean_pi <- sum(tvb_w*10^(tvb$logodds)/(1+10^(tvb$logodds)))
    tvb_mean_alpha <- mean(c(tvb$alpha%*%tvb_w))
    tvb_mean_log10odds <- sum(tvb_w*tvb$logodds)
    tvb_mean_logodds <- sum(tvb_w*log(10)*tvb$logodds)
    tvb_mean_sigb <- sum(tvb_w*sqrt(tvb$sa*tvb$sigma))
    mean_logodds[i] <- tvb_mean_logodds
    mean_sigb[i] <- tvb_mean_sigb
    mean_alpha[i] <- tvb_mean_alpha
    mean_pi <- tvb_mean_pi
    mean_pve <- mean(tvb$model.pve)
    # retmat[i,1] <- tvb_mean_logodds
    # retmat[i,2] <- tvb_mean_sigb
    # retmat[i,3] <- tvb_mean_alpha
    # retmat[i,4] <- tvb_mean_pi
    # retmat[i,5] <- as.integer(fgeneid[i])
  }
  return(data_frame(fgeneid=fgeneid,
                    mean_logodds=mean_logodds,
                    mean_sigb=mean_sigb,
                    mean_alpha=mean_alpha,
                    mean_pi=mean_pi,
                    mean_pve=mean_pve,
                    use_tsigb=!is.null(tsigb),
                    use_tlogodds=!is.null(tlogodds),
                    use_tresid=!is.null(tresid)))
  
}







glmnet_max <- function(X,ymat,glmnet_alpha,fgeneid){
  library(glmnet)
  library(dplyr)
  library(progress)
  retmat <- matrix(0,ncol(ymat),3)
  colnames(retmat) <- c("sigb","logodds","fgeneid")
  pb <- progress_bar$new(total=ncol(ymat))
  for(i in 1:ncol(ymat)){
      pb$tick()
    tcv <- cv.glmnet(X,ymat[,i],alpha=0.95)
    betas <- coef(tcv,s="lambda.1se")[-1,]
    retpi <- mean(betas!=0)
    retsigb <- sd(betas[betas!=0])
    retmat[i,1] <- retsigb
    retmat[i,2] <- retpi
    retmat[i,3] <- i
  }
  return(as_data_frame(retmat))
}


rssr_max <- function(R,betahat_mat,se_mat,sigb,logodds,tolerance=1e-3,lnz_tol=T,itermax=100,fgeneid=NULL){
  library(dplyr)
  library(progress)
  library(purrr)
  paramdf <- list(sigb=sigb,logodds=logodds) %>% cross_d()
  sigb <- paramdf$sigb
  logodds <- paramdf$logodds
  stopifnot(ncol(R)==nrow(betahat_mat),
            ncol(betahat_mat)==ncol(se_mat),
            length(logodds)==length(sigb))
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(betahat_mat)
  }
  
  retl <- list()
  p <- nrow(betahat_mat)
  retmat <- matrix(0,ncol(betahat_mat),4)
  colnames(retmat) <- c("lnZ","sigb","logodds","fgeneid")
  pb <- progress_bar$new(total=ncol(betahat_mat))
  for(i in 1:ncol(betahat_mat)){
    pb$tick()
    SiRiS <- SiRSi_d(R,Si=1/se_mat[,i])
    
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- (SiRiS%*%(alpha0*mu0))
    retdf <- grid_search_rss_varbvsr(SiRiS = SiRiS,
                                     sigma_beta = sigb,
                                     logodds = logodds,
                                     betahat = betahat_mat[,i],
                                     se = se_mat[,i],talpha0 = alpha0,
                                     tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                     tolerance = tolerance,itermax=100,verbose = F,
                                     lnz_tol = lnz_tol) %>% filter(lnZ==max(lnZ)) %>% slice(1)
    retmat[i,1] <- retdf$lnZ
    retmat[i,2] <- retdf$sigb
    retmat[i,3] <- retdf$logodds
    retmat[i,4] <- i
  }
  return(as_data_frame(retmat))
  
}



rssr_mean <- function(R,betahat_mat,se_mat,sigb,logodds,tolerance=1e-3,lnz_tol=T,itermax=100,fgeneid=NULL,tlogodds=NULL,tsigb=NULL){
  library(dplyr)
  library(rssr)
  library(progress)
  library(purrr)
  
  paramdf <- list(sigb=sigb,logodds=logodds) %>% cross_d()
  sigb <- paramdf$sigb
  logodds <- paramdf$logodds
  
  
  stopifnot(ncol(R)==nrow(betahat_mat),
            ncol(betahat_mat)==ncol(se_mat),
            length(logodds)==length(sigb))
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(betahat_mat)
  }
  
  
  
  
  
  
  ng <- ncol(betahat_mat)
  retl <- list()
  p <- nrow(betahat_mat)
  retmat <- matrix(0,ng,5)
  colnames(retmat) <- c("mean_logodds","mean_sigma","mean_alpha","mean_pi", "fgeneid")
  pb <- progress_bar$new(total=ng)
  if(!is.null(tsigb)){
    stopifnot(length(tsigb)==ng)
    # sigb <- tsigb
  }
  if(!is.null(tlogodds)){
    stopifnot(length(tlogodds)==ng)
    #  logodds <- tlogodds
  }
  
  for(i in 1:ng){
    pb$tick()
    SiRiS <- SiRSi_d(R,Si=1/se_mat[,i])
    
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- (SiRiS%*%(alpha0*mu0))
    
    if(!is.null(tlogodds)){
      logodds <- tlogodds[i]
    }
    if(!is.null(tsigb)){
      sigb <- tsigb[i]
    }
    retdf <- grid_search_rss_varbvsr(SiRiS = SiRiS,
                                     sigma_beta = sigb,
                                     logodds = logodds,
                                     betahat = betahat_mat[,i],
                                     se = se_mat[,i],talpha0 = alpha0,
                                     tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                     tolerance = tolerance,itermax=100,verbose = F,
                                     lnz_tol = lnz_tol) %>% 
      mutate(pi=exp(logodds)/(exp(logodds)+1),w=normalizeLogWeights(lnZ,na.rm=T)) %>% 
      summarise(mean_pi=sum(w*pi),mean_sigma=sum(w*sigb),mean_logodds=sum(w*logodds),mean_alpha=sum(w*(alpha_mean)))
    retmat[i,1] <- retdf$mean_logodds
    retmat[i,2] <- retdf$mean_sigma
    retmat[i,3] <- retdf$mean_alpha
    retmat[i,4] <- retdf$mean_pi
    retmat[i,5] <- i
  }
  return(as_data_frame(retmat))
}


rssr_pve <- function(R,betahat_mat,se_mat,sigb,logodds,tolerance=1e-3,lnz_tol=T,itermax=100,fgeneid=NULL,tlogodds=NULL,tsigb=NULL){
  library(dplyr)
  library(rssr)
  library(progress)
  
  
  stopifnot(ncol(R)==nrow(betahat_mat),
            ncol(betahat_mat)==ncol(se_mat),
            length(logodds)==length(sigb))
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(betahat_mat)
  }
  
  ng <- ncol(betahat_mat)
  retl <- list()
  p <- nrow(betahat_mat)
  retmat <- matrix(0,ng,5)
  colnames(retmat) <- c("mean_logodds","mean_sigma","mean_alpha","mean_pi", "fgeneid")
  pb <- progress_bar$new(total=ng)
  if(!is.null(tsigb)){
    stopifnot(length(tsigb)==ng)
    # sigb <- tsigb
  }
  if(!is.null(tlogodds)){
    stopifnot(length(tlogodds)==ng)
    #  logodds <- tlogodds
  }
  
  for(i in 1:ng){
    pb$tick()
    SiRiS <- SiRSi_d(R,Si=1/se_mat[,i])
    
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- (SiRiS%*%(alpha0*mu0))
    
    if(!is.null(tlogodds)){
      logodds <- tlogodds[i]
    }
    if(!is.null(tsigb)){
      sigb <- tsigb[i]
    }
    
    
    retdf <- grid_search_rss_varbvsr(SiRiS = SiRiS,
                                     sigma_beta = sigb,
                                     logodds = logodds,
                                     betahat = betahat_mat[,i],
                                     se = se_mat[,i],talpha0 = alpha0,
                                     tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                     tolerance = tolerance,itermax=100,verbose = F,
                                     lnz_tol = lnz_tol) %>% 
      mutate(pi=exp(logodds)/(exp(logodds)+1),w=normalizeLogWeights(lnZ,na.rm=T)) %>% 
      summarise(mean_pi=sum(w*pi),mean_sigma=sum(w*sigb),mean_logodds=sum(w*logodds),mean_alpha=sum(w*(alpha_mean)))
    retmat[i,1] <- retdf$mean_logodds
    retmat[i,2] <- retdf$mean_sigma
    retmat[i,3] <- retdf$mean_alpha
    retmat[i,4] <- retdf$mean_pi
    retmat[i,5] <- i
    
  }
  retdf <- as_data_frame(retmat)
  
  
  
  return(retdf)
}




tspve_all <- function(R,betahat_mat,se_mat,betamat,n,fgeneid=NULL){
  require(dplyr)
  ng <- ncol(betahat_mat)
  stopifnot(ng==ncol(betamat))
  if(is.null(fgeneid)){
    fgeneid=1:ng()
  }
  stopifnot(length(fgeneid)==ng)
  pvevec <- numeric(ng)
  for(j in 1:ng){
    bh <- betahat_mat[,j]
    seh <- se_mat[,j]
    b <- betamat[,j]
    siris <- SiRSi_d(R,Si = 1/seh)  
    spve <- t(b)%*%siris%*%b
    pvevec[j] <- c(spve)/n
  }
  return(data_frame(fgeneid=fgeneid,tspve=pvevec))
}




rssr_all <- function(R,betahat_mat,se_mat,sigb,logodds,tolerance=1e-3,lnz_tol=T,itermax=200,fgeneid=NULL,n=1,use_squarem=T){
  
  library(dplyr)
  library(rssr)
  library(progress)
  library(BBmisc)
  library(purrr)

  # paramdf <- list(sigb=sigb,logodds=logodds) %>% cross_d() %>% distinct()
  # sigb <- paramdf$sigb
  # logodds <- paramdf$logodds  
  
  stopifnot(ncol(R)==nrow(betahat_mat),
            ncol(betahat_mat)==ncol(se_mat),
            length(sigb)==length(logodds))
  
  if(is.null(fgeneid)){
    fgeneid <- 1:ncol(betahat_mat)
  }

  
  num_sigb <- length(sigb)
  num_logodds <- length(logodds)
  ng <- ncol(betahat_mat)

  retl <- list()
  p <- nrow(betahat_mat)
  
  
  pb <- progress_bar$new(total=ng)
  retdfl <- list()
  nind=n
  for(i in 1:ng){

    SiRiS <- SiRSi_d(R,Si=1/se_mat[,i])
    
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- (SiRiS%*%(alpha0*mu0))
    mfgeneid <- fgeneid[i]
     if(use_squarem){
       retdfl[[i]] <- grid_search_rss_varbvsr(SiRiS = SiRiS,
                                              sigma_beta = sigb,
                                              logodds = logodds,
                                              betahat = betahat_mat[,i],
                                              se = se_mat[,i],talpha0 = alpha0,
                                              tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                              tolerance = tolerance,itermax=itermax,verbose = F,
                                              lnz_tol = lnz_tol) %>% 
         mutate(pi=exp(logodds)/(exp(logodds)+1),pve=pve/nind,fgeneid=mfgeneid)
     }else{
       retdfl[[i]] <- grid_search_rss_varbvsr_naive(SiRiS = SiRiS,
                                              sigma_beta = sigb,
                                              logodds = logodds,
                                              betahat = betahat_mat[,i],
                                              se = se_mat[,i],talpha0 = alpha0,
                                              tmu0 = mu0,tSiRiSr0 = SiRiSr,
                                              tolerance = tolerance,itermax=itermax,verbose = F,
                                              lnz_tol = lnz_tol) %>% 
         mutate(pi=exp(logodds)/(exp(logodds)+1),pve=pve/nind,fgeneid=mfgeneid)
       
       
     }
    pb$tick()
  }
  return(bind_rows(retdfl))
#    data_frame(lnZ=lnZvec,pi=pivec,alpha_mean=alpha_meanvec,sigb=sigbvec,fgeneid=fgeneidvec,pve=pvevec))
}





sim_eff <- function(p,neff,tsigb){

  stopifnot(neff<=p,neff>=1)
  Z <- sample(1:p,neff,replace=F)
  
  beta <- numeric(p)
  beta[Z] <- rnorm(length(Z),mean=0,sd=tsigb)
  return(beta)
}

calc_resid <- function(X,beta,pve){
  part_1 = (1-pve) / pve;
  n <- nrow(X)
  xb     = X%*%beta
  stopifnot(all(!is.na(part_1)))
  part_2 = c(crossprod(xb) / n)
  stopifnot(all(!is.na(part_2)))
  ret <- sqrt(part_1*part_2)
  stopifnot(!is.na(ret))
  return(ret)
}

sim_y_pve <- function(X,beta,sigma){
  n  <- nrow(X)
  return(X%*%beta+rnorm(p,0,sigma))  
}





sim_betamat <- function(tparamdf,p){
  require(progress)
  np <- nrow(tparamdf)
  betamat <- matrix(0,p,np)
  pb <- progress_bar$new(total=np)
  tparamdf <- mutate(tparamdf,neff=ceiling(tpi*p))
  for(i in 1:ncol(betamat)){
    pb$tick()
    betamat[,i] <- sim_eff(p=p,neff = tparamdf$neff[i],tsigb=tparamdf$tsigb[i])
  }
  return(betamat)
}


calc_residvec <- function(tparamdf,SNP,betamat){
  require(progress)
  n <- nrow(SNP)
  np <- nrow(tparamdf)
  p <- ncol(betamat)
  residvec <- numeric(np)
  pb <- progress_bar$new(total=np)
  for(i in 1:np){
    pb$tick()
    residvec[i] <- calc_resid(SNP,betamat[,i],tparamdf$tpve[i])
  }
  return(residvec)
}

sim_residmat <- function(n,residvec){
  require(progress)
  np <- length(residvec)
  residmat <- matrix(0,n,np)
  pb <- progress_bar$new(total=np)
  for(i in 1:np){
    residmat[,i] <- rnorm(n,mean = 0,sd = residvec[i])
  }
  return(residmat)
}
# 
#   stopifnot(file.exists(eqtl_snpfile))
#   snpA <-read_2d_index_h5(eqtl_snpfile,"SNPdata","genotype",eqtl_snplist)
#   return( c(snpA%*%beta))
#   
# }

sim_ys <- function(eqtl_snpfile,eqtl_snplist,betamat){
  
  stopifnot(file.exists(eqtl_snpfile))
  snpA <-read_2d_index_h5(eqtl_snpfile,"SNPdata","genotype",eqtl_snplist)
  return( c(snpA%*%beta))
  
}



gen_LD <- function(ld_snpfile,ld_snplist,m=85,Ne=11490.672741,cutoff=1e-3,mapA=NULL,snpA=NULL){
  #  ld_ldf<-read_attr(param_snpfile,as.character(snp_chunk),"1kg_filepath")
  stopifnot(file.exists(ld_snpfile))
  #  ld_snplist <- read_vec(param_snpfile,paste0("/",snp_chunk,"/1kg"))
  stopifnot(!is.unsorted(ld_snplist))
  if(is.null(snpA)){
    snpA <-read_2d_index_h5(ld_snpfile,"SNPdata","genotype",ld_snplist)
  }
  if(is.null(mapA)){
    mapA <- read_1d_index_h5(ld_snpfile,"SNPinfo","map",ld_snplist)
    # mapA <-read_vec(ld_snpfile,"/SNPinfo/map")[ld_snplist]
  }
  stopifnot(!is.unsorted(mapA))
  
  stopifnot(ncol(snpA)==length(ld_snplist),length(mapA)==length(ld_snplist))
  
  sp_R <- calcLD(snpA,mapA,m,Ne,cutoff)
  diag(sp_R) <- 1
  return(sp_R)
}


gen_LD_sp <- function(ld_snpfile,ld_snplist,m=85,Ne=11490.672741,cutoff=1e-3,mapA=NULL,snpA=NULL){
  #  ld_ldf<-read_attr(param_snpfile,as.character(snp_chunk),"1kg_filepath")
  stopifnot(file.exists(ld_snpfile))
  #  ld_snplist <- read_vec(param_snpfile,paste0("/",snp_chunk,"/1kg"))
  
  if(is.null(snpA)){
    snpA <-read_2d_index_h5(ld_snpfile,"SNPdata","genotype",ld_snplist)
  }
  if(is.null(mapA)){
    mapA <- read_1d_index_h5(ld_snpfile,"SNPinfo","map",ld_snplist)
    # mapA <-read_vec(ld_snpfile,"/SNPinfo/map")[ld_snplist]
  }
  
  stopifnot(ncol(snpA)==length(ld_snplist),length(mapA)==length(ld_snplist))
  
  sp_R <- sp_calcLD(snpA,mapA,m,Ne,cutoff)
  return(sp_R)
}


gen_LD_symm <- function(ld_snpfile,ld_snplist,m=85,Ne=11490.672741,cutoff=1e-3,mapA=NULL,snpA=NULL){
  #  ld_ldf<-read_attr(param_snpfile,as.character(snp_chunk),"1kg_filepath")
  stopifnot(file.exists(ld_snpfile))
  #  ld_snplist <- read_vec(param_snpfile,paste0("/",snp_chunk,"/1kg"))
  
  if(is.null(snpA)){
    snpA <-read_2d_index_h5(ld_snpfile,"SNPdata","genotype",ld_snplist)
  }
  if(is.null(mapA)){

    mapA <-read_1d_index_h5(ld_snpfile,groupname = "SNPinfo",dataname = "map",indvec = ld_snplist)
  }
  
  stopifnot(ncol(snpA)==length(ld_snplist),length(mapA)==length(ld_snplist))
  
  sp_R <- sp_calcLD_symm(snpA,mapA,m,Ne,cutoff)
  return(sp_R)
}





gen_eQTL <- function(eqtl_snpfile,eqtl_expfile,eqtl_explist,eqtl_snplist,scale_ortho_exp=T){
  
  library(dplyr)
  stopifnot(file.exists(eqtl_expfile))
  exp_fgeneid <- read_vec(eqtl_expfile,"EXPinfo/fgeneid")[eqtl_explist]
  eqtl_snpdat <- read_2d_index_h5(eqtl_snpfile,"SNPdata","orthogenotype",eqtl_snplist)
  eqtl_expdat <- read_2d_index_h5(eqtl_expfile,"EXPdata","orthoexpression",eqtl_explist)

  eqtl_covdat <- read_2d_mat_h5(eqtl_expfile,"Covardat","covariates")
  ortho_covdat <- rssr_orthogonalize_covar(cbind(1,eqtl_covdat))
  
  
  # teqtl_df <- list()
  ortho_genodat <- scale(orthogonalize_data(eqtl_snpdat,ortho_covdat),center=T)
  ortho_expdat <- scale(orthogonalize_data(eqtl_expdat,ortho_covdat),center=T)
  
  betahat_mat <- map_beta_exp(ortho_genodat,ortho_expdat)
  semat <- map_se_exp(ortho_genodat,expression = ortho_expdat,betahat = betahat_mat)
  return(list("betahat"=betahat_mat,"se"=semat))
}
  
  
run_RSS <- function(param_snpfile,param_expfile,snp_chunk,exp_chunk,tissue,
                    scale_ortho_exp=F,m=85,Ne=11490.672741,cutoff=1e-3,
                    logodds=NULL,sigb=NULL,exp_chunksize=1,sparseLD=FALSE){
  require(Matrix)
  require(dplyr)
  require(rssr)
  require(RcppEigenH5)
  library(h5)
  library(BBmisc)
  stopifnot(file.exists(param_snpfile),
            file.exists(param_expfile))
  
  
  eqtl_snpfile <-read_attr(param_snpfile,as.character(snp_chunk),paste0(tissue,"_filepath"))
  eqtl_snp_chromosome <-read_attr(param_snpfile,as.character(snp_chunk),"chromosome")
  stopifnot(file.exists(eqtl_snpfile))
  
  eqtl_snplist <-read_vec(param_snpfile,paste0(snp_chunk,"/",tissue))
  eqtl_snpdat <- read_2d_index_h5(eqtl_snpfile,"SNPdata","genotype",eqtl_snplist)
  
  
  if(exp_chunksize==1){
    eqtl_explist <-read_vec(param_expfile,paste0("all/",tissue))[exp_chunk]
    exp_fgeneid <- read_vec(param_expfile,"all/fgeneid")[exp_chunk]
    eqtl_expfile <- read_attr(param_expfile,"all",paste0(tissue,"_filepath"))
  }else{
    if(!exists_group(param_expfile,as.character(exp_chunk))){
      eqtl_explist <-chunk(read_vec(param_expfile,paste0("all/",tissue)),chunk.size = exp_chunksize)[[exp_chunk]]
      eqtl_fgeneid <-chunk(read_vec(param_expfile,paste0("all/fgeneid")),chunk.size = exp_chunksize)[[exp_chunk]]
      eqtl_expfile <- read_attr(param_expfile,"all",paste0(tissue,"_filepath"))
    }else{
      eqtl_explist <-read_vec(param_expfile,paste0(exp_chunk,"/",tissue))
      exp_fgeneid <- read_vec(param_expfile,paste0(exp_chunk,"/fgeneid"))
      eqtl_expfile <- read_attr(param_expfile,as.character(exp_chunk),paste0(tissue,"_filepath"))
    }
  }
  
  stopifnot(file.exists(eqtl_expfile))
  eqtl_expdat <- read_2d_index_h5(eqtl_expfile,"EXPdata","expression",eqtl_explist)
  
  
  eqtl_covdat <- read_2d_mat_h5(eqtl_expfile,"Covardat","covariates")
  ortho_covdat <- rssr_orthogonalize_covar(cbind(1,eqtl_covdat))
  
  
  teqtl_df <- list()
  ortho_genodat <- orthogonalize_data(eqtl_snpdat,ortho_covdat)
  ortho_expdat <- scale(orthogonalize_data(eqtl_expdat,ortho_covdat),center=scale_ortho_exp,scale=scale_ortho_exp)
  
  
  for(i in 1:ncol(eqtl_expdat)){
    teqtl_df[[i]] <-map_eqtl_lm(ortho_genodat,ortho_expdat[,i]) %>% mutate(fgeneid=exp_fgeneid[i])
  }
  
  summary_stats <- bind_rows(teqtl_df)
  
  ld_ldf<-read_attr(param_snpfile,as.character(snp_chunk),"1kg_filepath")
  stopifnot(file.exists(ld_ldf))
  ld_snplist <- read_vec(param_snpfile,paste0("/",snp_chunk,"/1kg"))
  
  
  snpA <-read_2d_index_h5(ld_ldf,"SNPdata","genotype",ld_snplist)
  mapA <-read_vec(ld_ldf,"/SNPinfo/map")[ld_snplist]
  
  stopifnot(ncol(snpA)==length(ld_snplist),length(mapA)==length(ld_snplist))
  
  
  if(sparseLD){
    R <- sp_calcLD(snpA,mapA,m,Ne,cutoff)
  }else{
    R <- calcLD(snpA,mapA,m,Ne,cutoff)
  }
  
  # R <- RColumbo::gen_sparsemat(RColumbo::calcLD(hmata=snpA,hmatb = snpA,mapa = mapA,mapb = mapA,m = m,Ne=Ne,cutoff = cutoff,isDiag = T),istart=1,jstart=1,nSNPs = ncol(snpA),makeSymmetric = T)
  
  
  
  rm(snpA,mapA,eqtl_covdat)
  
  fit_dfl <- list()
  summ_statl <- split(summary_stats,summary_stats$fgeneid)
  for(i in 1:length(summ_statl)){
    cat(i,"of ",length(summ_statl),"\n")
    se <- summ_statl[[i]][["se"]]
    betahat <- summ_statl[[i]][["betahat"]]
    if(sparseLD){
      SiRiS <- SiRSi(R,Si = 1/se)
    }else{
      SiRiS <- SiRSi_d(R,Si=1/se)
    }
    
    p <- length(betahat)
    alpha0 <- ralpha(p = p)
    mu0 <-rmu(p)
    SiRiSr <- SiRiS%*%(alpha0*mu0)
    if(sparseLD){
      fit_df <- grid_search_rss_varbvsr_sp(SiRiS = SiRiS,sigma_beta = sigb,logodds=logodds,betahat = betahat,
                                        se = se,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = SiRiSr@x,tolerance = 1e-3,itermax = 100,verbose = T,lnz_tol = T)
    }else{
      fit_df <- grid_search_rss_varbvsr(SiRiS = SiRiS,sigma_beta = sigb,logodds=logodds,betahat = betahat,
                                        se = se,talpha0 = alpha0,tmu0 = mu0,tSiRiSr0 = SiRiSr@x,tolerance = 1e-3,itermax = 100,verbose = T,lnz_tol = T)
    }
    
    fit_dfl[[i]] <- fit_df %>% mutate(fgeneid=exp_fgeneid[i])
    gc()
  }
  fit_df <- bind_rows(fit_dfl) %>% mutate(snp_chunk=snp_chunk,exp_chunk=exp_chunk,snp_chromosome=eqtl_snp_chromosome)
  return(fit_df)
}


libpath <- function() {
  cat(sprintf(
    "%s/RSSReQTL/libs/RSSReQTL%s",
    installed.packages()["RSSReQTL","LibPath"][1],
    .Platform$dynlib.ext
  ))
}