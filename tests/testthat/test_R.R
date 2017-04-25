context("compare to R")


test_that("QR decomposition works as in RColumbo",{
  
  rcovarmat <- cbind(1,matrix(rnorm(500*10),500,10))
  tqr <- qr(rcovarmat)
  RQ <- qr.Q(tqr)
  rss_covar <- rssr_orthogonalize_covar(rcovarmat)
  expect_equal(rss_covar,RQ)
})

test_ortho <- function(X,covx){
  return(X-t((t(X)%*%covx)%*%t(covx)))
}

test_that("orthogonalization of data works as in RColumbo",{
  rcovarmat <- cbind(1,matrix(rnorm(500*10),500,10))
  rdatamat <- matrix(rnorm(500*20),500,20)
  rss_covar <- rssr_orthogonalize_covar(rcovarmat)
  R_ortho <- test_ortho(rdatamat,rss_covar)
  rss_ortho <- orthogonalize_data(rdatamat,rss_covar)
  expect_equal(rss_ortho,rc_ortho)
  
})


test_that("eQTL mapping works as in RColumbo",{
  tgeno <- runif(500)
  texp <- runif(500)
  
  rc_eqtl <- as.numeric(RColumbo::map_eqtl_lm(t(t(tgeno)),t(t(texp)),scale_ortho_exp = F)[1,])
  rss_eqtl <- map_eqtl_lm(tgeno,texp)
  expect_equal(rss_eqtl,rc_eqtl,tolerance=1e-4)
})



test_that("eQTL mapping (for multiple traits & multiple SNPs) works as in R",{
  
  
  R_lm <- function(X,Y){
    n <- nrow(X)
    s <- ncol(X)
    g <- ncol(Y)
    
    bmat <- matrix(0,s,g)
    for(i in 1:s){
      for(j in 1:g){
        bmat[i,j] <- c(coef(lm(Y[,j]~X[,i]+0)))
      }
    }
    return(bmat)
  }
  R_se <- function(X,Y){
    n <- nrow(X)
    s <- ncol(X)
    g <- ncol(Y)
    
    semat <- matrix(0,s,g)
    for(i in 1:s){
      for(j in 1:g){
        semat[i,j] <- coef(summary(lm(Y[,j]~X[,i]+0)))[1,2]
      }
    }
    return(semat)
  }
  R_resid <- function(X,Y){
    n <- nrow(X)
    s <- ncol(X)
    g <- ncol(Y)
    
    res_array <- array(0,c(n,s,g))
    for(i in 1:s){
      for(j in 1:g){
        tlm <- lm(Y[,j]~X[,i]+0)
        res_array[,i,j] <- resid(tlm)
      }
    }
    return(res_array)
  } 
  
  
  n <- 500
  g <- 50
  s <- 100
  tgeno <- matrix(runif(n*s),n,s)
  texp <- matrix(runif(n*g),n,g)
  
  mbhat <- map_beta_exp(tgeno,texp)
  mse <- map_se_exp(tgeno,texp,mbhat)
  r_bmat <- R_lm(tgeno,texp)
  r_semat <- R_se(tgeno,texp)
  expect_equal(r_bmat,mbhat)
  
  expect_equal(mse,r_semat,tolerance=1e-4)
})



test_that("calculation of ntheta works as in RColumbo",{
  m <- 85
  rc_theta <- RColumbo::calc_nmsum(m)
  rss_theta <- calc_nmsum(m)
  expect_equal(rss_theta,rc_theta)
})

test_that("calculation of distamce matrix works as in RColumbo",{
  tmap <- cumsum(runif(100))
  rc_dist <- RColumbo::wrap_ip_dist(tmap,tmap,F)
  rss_dist <- calc_dist(tmap)
  expect_equal(rss_dist[upper.tri(rss_dist)],rc_dist[upper.tri(rc_dist)])
})

test_that("calculation of covariance matrix works as expected",{
  tdat <- matrix(runif(9*8),8,9)
  R_cov <- cov(tdat)
  R_cov[lower.tri(R_cov,diag=T)] <- 0
  diag(R_cov) <- 1
  rss_cov <- calc_cov(tdat)
  expect_equal(rss_cov,R_cov)
  
})

test_that("Calculation of variance works as expected",{
  tdat <- matrix(runif(9*8),9,8)
  R_vars <- apply(tdat,2,var)
  rss_vars <- calc_variance(tdat)
  expect_equal(rss_vars,R_vars)
})

test_that("Shrinkage works as expected",{
  m=85
  tmap <- cumsum(runif(100))
  thmat <- matrix(runif(100*m),m,100) 
  
  Ne=11490.672741
  nmsum <- calc_nmsum(m)
  theta =(1/nmsum)/(2*m+1/nmsum)
  cutoff=1e-3
  thc <- cov(thmat)
  thc[lower.tri(thc)]<-0
  diag(thc) <- 1
  tdist <- calc_dist(tmap)
  Rshrink <- RColumbo::wrap_compute_shrinkage(tdist,thc,thmat,theta = theta,m = m,Ne=Ne,cutoff = cutoff,isDiag = T)
  rss_shrink <- calc_shrinkage(distmat = tdist,S = thc,hmata = thmat,theta = theta,m = m,Ne = Ne,cutoff = cutoff)
  expect_equal(rss_shrink,Rshrink,tolerance=1e-4)
})

test_that("cov_2_cor works as in RColumbo",{
  thmat <- scale(matrix(runif(100*m),m,100),center=T,scale = F)
  
  tcmat <- cov(thmat)
  covs <- apply(thmat,2,var)
  rc_c2c <- RColumbo::cov_2_cor(tcmat,t(covs))
  R_c2c <- cov2cor(tcmat)
  expect_equal(R_c2c,rc_c2c)
  rss_c2c <- cov_2_cor(tcmat,covs)
  expect_equal(R_c2c,rss_c2c)
})

test_that("calcLD works as in RColumbo",{
  library(dplyr)
  m=85
  tmap <- cumsum(runif(100))
  thmat <- matrix(runif(100*m),m,100) 
  Ne=11490.672741
  cutoff=1e-3
  
  rc_ld <- RColumbo::calcLD(thmat,thmat,tmap,tmap,m,Ne,cutoff,T)
  # rc_ld[is.na(rc_ld)] <- 0
  rss_ld <- calcLD(thmat,tmap,m,Ne,cutoff)
  rss_ld[lower.tri(rss_ld)] <- 0
  gen_rsid <- paste0("rs",1:100)
  lddf <- ld2df(rss_ld,gen_rsid)
  ldcols <- gen_rsid[col(rss_ld)[upper.tri(rss_ld)]]
  ldrows <- gen_rsid[row(rss_ld)[upper.tri(rss_ld)]]
  ldv <- rss_ld[upper.tri(rss_ld)]
  rlddf <- data_frame(rowsnp=ldrows,colsnp=ldcols,r=ldv)
  expect_equal(rlddf,lddf)
  expect_equal(rss_ld[upper.tri(rss_ld)],rc_ld[upper.tri(rc_ld)],tolerance=1e-4)
})


test_that("eqtl mapping (with covariates) works as in RColumbo",{
  rcovarmat <- matrix(rnorm(500*10),500,10)
  rgenodat <- matrix(rnorm(500*20),500,20)
  rexpdat <- matrix(rnorm(500),500,1)
  
  rss_covar <- rssr_orthogonalize_covar(cbind(1,rcovarmat))
  rc_eqtl <- RColumbo::map_eqtl_lm(rgenodat,rexpdat,rcovarmat,scale_ortho_exp=T)
  ortho_geno <- orthogonalize_data(rgenodat,rss_covar)
  ortho_exp <- scale(orthogonalize_data(rexpdat,rss_covar),center=T,scale=T)
  
  rss_eqtl <- map_eqtl_lm(ortho_geno,ortho_exp)
  expect_equal(rss_eqtl$betahat,rc_eqtl$betahat)
  expect_equal(rss_eqtl$se,rc_eqtl$se,tolerance=1e-4)
})