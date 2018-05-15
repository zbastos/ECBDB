source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(“limma”)
library(Biobase)
library(GEOquery)
gds5047 <- getGEO('GDS5047', destdir=".")

tabgds = Table(gds5047) #tabela do dataset
colnames(tabgds) #colunas
dim(tabgds) #dimensão
tabgds[1,1:5] #primeiras 5 colunas da primeira linha

Meta(gds5047)$sample_organism
Meta(gds5047)$title

#conversão para expressionSet
expset = GDS2eSet(gds5047,do.log2=FALSE,GPL=NULL,AnnotGPL=TRUE,getGPL=TRUE)
expset
dim(expset)
featureNames(expset)[1:5]
sampleNames(expset)[1:5]
varMetadata(expset) #metadados

expset$sample
expset$agent
expset$individual
expset$description