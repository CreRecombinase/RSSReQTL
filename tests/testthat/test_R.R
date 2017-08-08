context("compare to R")


test_that("QR decomposition works as I think it does",{
  
  rcovarmat <- cbind(1,matrix(rnorm(500*10),500,10))
  tqr <- qr(rcovarmat)
  RQ <- qr.Q(tqr)
  rss_covar <- rssr_orthogonalize_covar(rcovarmat)
  expect_equal(rss_covar,RQ)
})

test_ortho <- function(X,covx){
  return(X-t((t(X)%*%covx)%*%t(covx)))
}

test_that("orthogonalization of data works as I think it does",{
  rcovarmat <- cbind(1,matrix(rnorm(500*10),500,10))
  rdatamat <- matrix(rnorm(500*20),500,20)
  rss_covar <- rssr_orthogonalize_covar(rcovarmat)
  R_ortho <- test_ortho(rdatamat,rss_covar)
  rss_ortho <- orthogonalize_data(rdatamat,rss_covar)
  expect_equal(rss_ortho,R_ortho)
  
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

test_that("new way of computing see works",{
  
  
  test_se <-   function(genotype,expression,betahat){
  
    yh <- genotype%*%betahat
    byh <- kronecker(betahat,genotype)
    bresid 
    
    
    
  }
})


test_that("calculation of covariance matrix works as expected",{
  tdat <- matrix(runif(9*8),8,9)
  R_cov <- cov(tdat)
  rss_cov <- calc_cov(tdat)
  expect_equal(rss_cov,R_cov)
  
})

test_that("Calculation of variance works as expected",{
  tdat <- matrix(runif(9*8),9,8)
  R_vars <- apply(tdat,2,var)
  rss_vars <- calc_variance(tdat)
  expect_equal(rss_vars,R_vars)
})

test_that("LD shrinkage estimator throws error when the cummap isn't sorted",{
  m <- 100
  Ne <- 10000
  n <- 100
  p <- 500
  cutoff <- 1e-3
  tmap <- sample(cumsum(runif(p)/10))
  
  Hpanel <- matrix(sample(c(0,1),n*2*p,replace=T),n*2,p)
  # mfile <- system.file("m_files/run_install.m",package="rssr")
  expect_error(calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff),"Recombination map must be non-decreasing\n")
})



test_that("LD shrinkage estimator doesn't throw error when the cummap is non-decreasing (but not strictly increasing)",{
  m <- 100
  Ne <- 10000
  n <- 100
  p <- 500
  cutoff <- 1e-3
  tmap <- cumsum(runif(p)/10)
  tmap[1] <- tmap[2]
  
  Hpanel <- matrix(sample(c(0,1),n*2*p,replace=T),n*2,p)
  # mfile <- system.file("m_files/run_install.m",package="rssr")
  mLD <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
})





