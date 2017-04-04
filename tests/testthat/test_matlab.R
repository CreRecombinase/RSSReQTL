context("RSS")





test_that("LD shrinkage estimators work as expected",{
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
  Rmsig[lower.tri(Rmsig)] <- 0
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig[lower.tri(Rsig)] <- 0
  expect_equal(Rsig,Rmsig,tolerance=1e-5)
})



test_that("LD shrinkage estimators work as expected",{
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
  Rmsig[lower.tri(Rmsig)] <- 0
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig[lower.tri(Rsig)] <- 0
  expect_equal(Rsig,Rmsig,tolerance=1e-5)
})

test_that("LD shrinkage estimators work the same real data ",{
  
  m=85
  Ne=11490.672741
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
  Rmsig[lower.tri(Rmsig)] <- 0
  Rsig <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig[lower.tri(Rsig)] <- 0
  expect_equal(Rsig,Rmsig,tolerance=1e-5)
})



test_that("LD shrinkage estimators give similar results for genotype and haplotype info",{
  
  m=85
  Ne=11490.672741
  cutoff=1e-3
  data("genomat")
  data("mapdat")
  Hpanel <- haplomat
  tmap <- mapdat
  

  Rsig_h <- calcLD(hmata = Hpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig_h[lower.tri(Rsig_h)] <- 0
  data("genomat")
  Gpanel <- genomat
  Rsig_g <- calcLD(hmata = Gpanel,mapa = tmap,m = m,Ne = Ne,cutoff = cutoff)
  Rsig_g[lower.tri(Rsig_g)] <- 0

  expect_equal(Rsig_h,Rsig_g)
  
  
})




