library(biomaRt)
ensembl=useMart("ensembl")

listDatasets(ensembl)
###################
###Append description
###################
ensembl_yeast=useDataset("scerevisiae_gene_ensembl", mart=ensembl)

yeast_omim=read.csv(file="learning/biomaRt/Yeast_homo_in_any1.csv",sep=";",header=T)
yeast_omim=as.matrix(yeast_omim)
getBM(attributes=c('ensembl_gene_id','description'),filters='ensembl_gene_id',values=ensembl_id,mart=ensembl_yeast)
description=getBM(attributes=c('ensembl_gene_id','description'),filters='ensembl_gene_id',values=ensembl_id,mart=ensembl_yeast)

description=getBM(attributes=c('ensembl_gene_id','description'),filters='ensembl_gene_id',values=ensembl_id[i],mart=ensembl_yeast)
#yeast_omim_description=cbind(yeast_omim,apply(yeast_omim,2,function(x) getBM(attributes=c('ensembl_gene_id','description'),filters='ensembl_gene_id',values=yeast_omim[x,2],mart=ensembl_yeast)))
yeast_omim_description=c()
for(i in 1:671) {
  yeast_omim_description=c(yeast_omim_description,getBM(attributes=c('ensembl_gene_id','description'),filters='ensembl_gene_id',values=yeast_omim[i,2],mart=ensembl_yeast))
}
