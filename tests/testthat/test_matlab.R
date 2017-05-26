context("RSS")

test_that("theta is computed correctly",{
  m <- 100
  nmsum = sum(1 / (1:(2*m-1)))
  theta = (1/nmsum) / (2*m + 1/nmsum)
  expect_equal(calc_theta(m),theta)
})


test_that("LD shrinkage estimators work as expected on simulated data",{
  m <- 100
  Ne <- 10000
  n <- 100
  p <- 500
  cutoff <- 1e-3
  tmap <- cumsum(runif(p)/10)
  
  Hpanel <- matrix(sample(c(0,1),n*2*p,replace=T),n*2,p)
  # mfile <- system.file("m_files/run_install.m",package="rssr")
  mdir <- system.file("m_files",package="RSSReQTL")
  
  #change to the directory with the .m files in Octave
  library(RcppOctave)
  .CallOctave('cd',mdir)
  msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
  
  Rmsig <- cov2cor(msig)
  # Rmsig[lower.tri(Rmsig)] <- 0
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  
  
  # Rsig[lower.tri(Rsig)] <- 0
  expect_equal(Rsig,Rmsig)
})

test_that("LD shrinkage estimators work the same real data ",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("haplomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  mdir <- system.file("m_files",package="RSSReQTL")  
  #change to the directory with the .m files in Octave
  library(RcppOctave)
  .CallOctave('cd',mdir)
  msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
  Rmsig <- cov2cor(msig)
  # Rmsig[lower.tri(Rmsig)] <- 0
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  # Rsig[lower.tri(Rsig)] <- 0
  expect_equal(Rsig,Rmsig)
})





test_that("LD shrinkage estimators work the same real (larger) data ",{
  
  m=85
  Ne=1490.672741
  cutoff=1e-3
  data("bighaplo")
  data("bigmap")
  Hpanel <- bighaplo
  tmap <- bigmap
  
  
  mdir <- system.file("m_files",package="RSSReQTL")  
  #change to the directory with the .m files in Octave
  library(RcppOctave)
  .CallOctave('cd',mdir)
  msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
  Rmsig <- cov2cor(msig)
  
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  
  expect_equal(Rsig,Rmsig)
})


test_that("LD shrinkage estimators work the same real (larger) data and a big cutoff ",{
  
  m=85
  Ne=1490.672741

  data("bighaplo")
  data("bigmap")
  Hpanel <- bighaplo
  tmap <- bigmap
  
  
  mdir <- system.file("m_files",package="RSSReQTL")  
  #change to the directory with the .m files in Octave
  library(RcppOctave)
  .CallOctave('cd',mdir)
  msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
  Rmsig <- cov2cor(msig)
  
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  
  expect_equal(Rsig,Rmsig)
  cutoff=.5
  msig <- .CallOctave('shrink_cov',m,Ne,tmap,Hpanel,cutoff)
  Rmsig <- cov2cor(msig)
  
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  expect_equal(Rsig,Rmsig)
  evd <- eigen(Rsig)
  
})





test_that("LD shrinkage estimators give similar results for sparse and dense data",{
  
  m=85
  Ne=11490.672741
  cutoff=1e-3
  data("mapdat")
  data("haplomat")
  Hpanel <- haplomat
  tmap <- mapdat
  
  
  Rsig_h_d <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig_h_s <- as.matrix(sp_calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff))
  expect_equivalent(Rsig_h_d,Rsig_h_s)
  
  
  
})




