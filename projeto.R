source("http://bioconductor.org/biocLite.R")
library(Biobase)
library(gtools)
library(GEOquery)
library(limma)
library(AnnotationDbi)
library(genefilter)
library(gplots)
library(e1071)
library(caret)
library(class)
library(nnet)
library(rpart)

#Download GDS file, put it in the current directory, and load it:
gds5047 <- getGEO('GDS5047', destdir=".")

#Expression set
eset <- GDS2eSet(gds5047, do.log2=TRUE)
dim(eset)
eset

#Matriz de dados de expressao
exp=exprs(eset)
dim(exp)

#informacao sobre a experiencia
abstract(eset)

#MetaDados
vars=pData(eset)
names(vars)

levels(vars$agent)
levels(vars$individual)


#Frequencia absoluta das amostras para as variaveis Agent,Individual
ftable(vars$agent, vars$individual, dnn = c("Agent", "Individual"))


####################################################Pre-processamento##################################
#Remover NA
eset_NA= eset[complete.cases(exp),]
dim(eset_NA)

##############################################Funcoes filtragem e normalizacao #########################
# Filtrar dados pela mediana
filter_sd<-function(dataset){
  exp2=exprs(dataset)
  sds=rowSds(exp2)
  m=median(sds)
  hist(sds,breaks=50,col="mistyrose",main = "Histograma dos desvios padrões",
       xlab = "Desvios padrões",
       ylab = "Frequência")
  abline(v=m,col="blue",lwd=4,lty=2)
  abline(v=m*2,col="red",lwd=4,lty=2)
  legend("topright",legend=c("mediana","dobro da mediana"), bty="n",fill = c("blue","red"))
  esetr = dataset[sds >= 3*median(sds),]
  return(esetr)
}

#Filtro racio max-min ---> Resulta em 0 features logo nao usamos
filter_max<-function(dataset){
  exp2=exprs(dataset)
  maximos=apply(exp2,1,max)
  minimos=apply(exp2,1,min)
  vl=maximos/minimos > 2
  esetmax=dataset[vl,]
  return(esetmax)
}

#Normalizacao
normalizar<-function(dataset){
  exp_m2 = exprs(dataset)
  return(scale(exp_m2))
}
##################################################################################################

esetr=filter_sd(eset_NA)
dim(esetr)
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

lmfitagent<-function(dataset){
  design= model.matrix(~dataset$agent)
  fit=lmFit(dataset,design)
  fit2=eBayes(fit)
  diff=topTable(fit2,coef = 2,1000)
  return(diff)
}

####################################################################################################

#tstudentrows
#results=tstudentagent(esetr,20)
#results

#LmFit e Bayes ----> resultados mais semelhantes ao do artigo

diff=lmfitagent(eset_NA)
threshold = foldchange2logratio(1.4)

#Genes
genes = diff[which(diff$adj.P.Val<0.05 & (diff$logFC > threshold | diff$logFC < -threshold)),]
genes=genes[,c(1,3,16,17,18,23:28)]
View(genes)

#Genes sobreexpressos:
genessobre = diff[which(diff$adj.P.Val<0.05 & (diff$logFC > threshold )),]
genessobre=genessobre[,c(1,3,16,17,18,23:28)]
View(genessobre)

#Genes subexpressos:
genessub = diff[which(diff$adj.P.Val<0.05 & (diff$logFC < -threshold )),]
genessub=genessub[,c(1,3,16,17,18,23:28)]
View(genessub)

##############################################Clustering######################################

clusterhierarch<-function(dataset,res){
    cl=exprs(dataset[res[,1],])
    corr=cor(cl,method="pearson")
    d=as.dist(1-corr)
    cl.hier <- hclust(d)
    plot(cl.hier, main = "Clustering Hier?rquico")
}

dist.fun = function(x) {
  return (as.dist (1 - cor(t (x), method = "pearson")))
}

clust.fun = function (x) {
  return (hclust (x))
}

color.map.tissue <- function(status) { if (status == "cocaine") "turquoise" else "chocolate1" }

heatmapTop<-function(dataset,res){
  cl=dataset[res[,1]]
  tissuecolorsh <- unlist(lapply(cl$agent, color.map.tissue))
  heatmap(exprs(cl),ColSideColors = tissuecolorsh, cexRow = 0.5, distfun = dist.fun, hclustfun = clust.fun, labRow = F, margins = c(7,7), 
          ylab = "Sondas", xlab = "Amostras", main = "Heatmap")
  legend("topright",legend=c("cocaine","control"), bty="n",fill = c("turquoise","chocolate1"))
}

##########################################################################################################

#Cluster plot
clusterhierarch(eset_NA,genes)
#Heatmap
genes2 = diff[which(diff$adj.P.Val<0.05),]
genes2 =genes2[,c(1,3,16,17,18,23:28)]
heatmapTop(eset_NA,genes2)

#############################################Analise preditiva###############################################

#KNN
model_knn = train(t(exprs(esetr)), esetr$agent, method = "knn", trControl=trainControl("cv", number = 5))
pred_knn = predict(model_knn, t(exprs(esetr)))
mk1=confusionMatrix(pred_knn, esetr$agent)
mk1$table; mk1$overal[1]

#Decision tree
model_tree = train(t(exprs(esetr)), esetr$agent, method = "rpart", trControl=trainControl("cv", number = 5))
pred_tree = predict(model_tree, t(exprs(esetr)))
mt2 = confusionMatrix(pred_tree, esetr$agent)
mt2$table; mt2$overal[1]

#SVM
model_svm = train(t(exprs(esetr)), esetr$agent, method = "svmLinear", trControl=trainControl("cv", number = 5))
pred_svm = predict(model_tree, t(exprs(esetr)))
ms3 = confusionMatrix(pred_tree, esetr$agent)
ms3$table; ms3$overal[1]


#ANN com particao
train = t(exprs(esetr[,11:50]))
test=t(exprs(esetr[,c(1:10,51:60)]))

ann = nnet(esetr$agent[11:50]~.,data.frame(train),size=3,MaxNWts=10000)
valores.prev.ann = predict(ann,data.frame(test),type="class")
valores.reais= esetr$agent[c(1:10,51:60)]
table(valores.prev.ann,valores.reais)
sum(valores.prev.ann == valores.reais)/length(valores.reais)
