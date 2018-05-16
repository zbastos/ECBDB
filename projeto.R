#source("http://bioconductor.org/biocLite.R")
#biocLite("GEOquery")
library(Biobase)
library(GEOquery)
library(limma)
library(hgug4112a.db)
library(AnnotationDbi)
library(genefilter)

#Download GDS file, put it in the current directory, and load it:
gds5047 <- getGEO('GDS5047', destdir=".")

#Organismo
Meta(gds5047)$sample_organism
#Experiencia
Meta(gds5047)$title

#Expression set
eset <- GDS2eSet(gds5047, do.log2=TRUE)
dim(eset)

#Matriz de dados de expressao
exp=exprs(eset)
dim(exp)
class(exp)

#Atributos dos metadatos
varMetadata(eset)

#informacao sobre a experiencia
abstract(eset)


#Remover NA
eset_NA= eset[complete.cases(exp),]

##############################################Funcoes filtragem e normalizacao #########################
#filtro flat patterns sd > 2*mediana do sd
filter_sd<-function(dataset){
  exp2=exprs(dataset)
  sds=rowSds(exp2)
  m=median(sds)
  hist(sds,breaks=50,col="mistyrose")
  abline(v=m,col="blue",lwd=4,lty=2)
  abline(v=m*2,col="red",lwd=4,lty=2)
  esetr = dataset[sds >= 3*median(sds),]
  return(esetr)
}

#filtro flat patterns racio max/min > 2 ---> Resulta em 0 features logo nao usamos
filter_max<-function(dataset){
  exp2=exprs(dataset)
  maximos=apply(exp2,1,max)
  minimos=apply(exp2,1,min)
  vl=maximos/minimos > 2
  esetmax=dataset[vl,]
  return(esetmax)
}

normalizar<-function(dataset){
  exp_m2 = exprs(dataset)
  return(scale(exp_m2))
}
##################################################################################################

esetr=filter_sd(eset_NA)
exprs(esetr)=normalizar(esetr)

######################################## Funcoes Expressao diferencial ####################################
tstudentagent<-function(dataset,nelements){
  table(dataset$agent)
  tt=rowttests(dataset,"agent")
  rank=order(tt$p.value)
  p20=rank[1:nelements]
  tt$p.value[p20]
  return(featureNames(dataset[p20]))
}

lmfitagent<-function(dataset,n){
  design= model.matrix(~dataset$agent)
  fit=lmFit(dataset,design)
  fit2=eBayes(fit)
  diff=topTable(fit2,coef = 2,n)
  return(diff)
}

####################################################################################################

#tstudentrows
#results=tstudentagent(esetr,20)
#results

#LmFit e Bayes ----> resultados mais semelhantes ao do artigo
results=lmfitagent(esetr,20)
results[,c(3,23,24,25,26,27,28)]

##############################################Funcao clustering######################################

clusterhierarch<-function(dataset,res,metodo){
    cl=dataset[res[,1]]
    eucD=dist(exprs(cl))
    cl.hier <- hclust(eucD,method=metodo)
    plot(cl.hier)
}

heatmapTop<-function(dataset,res){
  cl=dataset[res[,1]]
  heatmap(exprs(cl),labCol = F)
}

##########################################################################################################

#Cluster plot
clusterhierarch(esetr,results,"average")
#Heatmap
heatmapTop(esetr,results)

train = t(exprs(esetr[,11:50]))
test=t(exprs(esetr[,c(1:10,51:60)]))

library(class)
valores.previstos= knn(train,test,esetr$agent[11:50])
valores.reais= esetr$agent[c(1:10,51:60)]
sum(valores.previstos == valores.reais)/ length(esetr$agent[c(1:10,51:60)])
table(valores.previstos,valores.reais)
valores.previstos= knn(train,test,esetr$agent[11:50],k=3)
table(valores.previstos,valores.reais)

library(nnet)
ann = nnet(esetr$agent[11:50]~.,data.frame(train),size=3,MaxNWts=10000)
valores.prev.ann = predict(ann,data.frame(test),type="class")
table(valores.prev.ann,valores.reais)
sum(valores.prev.ann == valores.reais)/length(valores.reais)

library(rpart)
arv = rpart(esetr$agent[11:50]~.,data.frame((train)))
plot(arv,uniform=T,branch=0.4,margin=0.1,compress=T)
text(arv,use.n=T,cex=0.9)
classes.previstas.arv = predict(arv,data.frame(test),type="class")
sum(classes.previstas.arv == valores.reais)/length(valores.reais)