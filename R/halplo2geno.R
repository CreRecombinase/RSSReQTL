

# haplofile <- '/media/nwknoblauch/Data/GTEx/1kg_SNP_H5/EUR.chr1_1kg.h5'
# haplomat <- read_2d_index_h5(haplofile,"SNPdata","genotype",1:1000)
# mapdat <- read_1d_index_h5(haplofile,"SNPinfo","map",1:1000)
# devtools::use_data(mapdat)


haplo2geno <- function(haplomat){
  nhaps <- nrow(haplomat)
  nind <- nhaps/2
  p <- ncol(haplomat)
  genomat <- matrix(0,nind,p)
  ind <- gl(n = nind,k=2,length = nhaps)
  return(as.matrix(aggregate(haplomat,by=list(ind),FUN=sum)[,-1]))
  # use_data(haplomat,compress = "gzip")
  # use_data(genomat,compress="gzip")
  
  
  
}