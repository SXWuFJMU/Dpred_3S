# DPred_3S: Identifying dihydrouridine (D) modification on three species epitranscriptome based on multiple sequence-derived features 
Dihydrouridine is a conserved modification of tRNA among all three life domains. D modification enhances the flexibility of a single nucleotide base in spatial structure and is disease- and evolution-associated. Recent studies have also suggested the presence of dihydrouridine on mRNA. To identify D in epitranscirptome, we provided a prediction framework named “Dpred” based on machine learning approach. The optimal features were evaluated by F score and different features integration; our model achieved AUROC scores 0.933 and 0.912 for E. coli and S. pombe, respectively. The performances of different machine learning algorithms were also compared in this study. 

##Example for traning model of S.pombe D sites prediction<br/>
library(BSgenome)<br/>
library(dplyr)<br/>
library(e1071)<br/>
library(dplyr)<br/>
library(caret)<br/>
library(ROCR)<br/>
Pdata <- readDNAStringSet(paste0("/FA/Spombe_D_P.fa"))<br/>
Ndata <- readDNAStringSet(paste0("/FA/Spombe_D_N.fa"))<br/>

XXX1 <- c(138)<br/>
set.seed(XXX1)<br/>
Ndata_s <- Ndata[sample(1:length(Ndata))][1:(length(Pdata))]<br/>
set.seed(123)<br/>
Pdata <- Pdata[sample(1:length(Pdata))]<br/>

Ndata_s1 <- gsub("T","U",as.data.frame(Ndata_s)$x)<br/>
Pdata1 <- gsub("T","U",as.data.frame(Pdata)$x)<br/>

S_order_COMPOSI <- readRDS("/F1_order/S_order_COMPOSI.rds")<br/>
S_order_ChemicalProper <- readRDS("/F1_order/S_order_ChemicalProper.rds")<br/>
S_order_PseKNC <- readRDS("/F1_order/S_order_PseKNC.rds")<br/>

##feature encoding<br/>
source("/methods/composition.R")<br/>
source("/methods/PseKNC.R")<br/>
source("/methods/ChemicialProperty.R")<br/>
Ndata_F1 <- list()<br/>
Pdata_F1 <- list()<br/>
Pdata_F1[[1]] <- composition(Pdata1)[,S_order_COMPOSI]<br/>
Pdata_F1[[2]] <- ChemicalProperty(as.data.frame(Pdata)$x)[,S_order_ChemicalProper]<br/>
Pdata_F1[[3]] <- PseKNC(as.data.frame(Pdata)$x)[,S_order_PseKNC]<br/>
Ndata_F1[[1]] <- composition(Ndata_s1)[,S_order_COMPOSI]<br/>
Ndata_F1[[2]] <- ChemicalProperty(as.data.frame(Ndata_s)$x)[,S_order_ChemicalProper]<br/>
Ndata_F1[[3]] <- PseKNC(as.data.frame(Ndata_s)$x)[,S_order_PseKNC]<br/>


Pdata_F <- cbind(Pdata_F1[[1]],Pdata_F1[[2]],Pdata_F1[[3]])<br/>
Ndata_F <- cbind(Ndata_F1[[1]],Ndata_F1[[2]],Ndata_F1[[3]])<br/>

Train_nunber <- c(1:length(Pdata))[c(1:round(length(Pdata)/10*8))]<br/>
test_nunber <- c(1:length(Pdata))[-c(1:round(length(Pdata)/10*8))]<br/>


label_train <- c(rep("motif",length(Train_nunber)),rep("non",length(Train_nunber)))<br/>
label_train <- factor(label_train,labels=c("motif","non"))<br/>
testlabel <- c(rep(1,length(test_nunber)),rep(0,length(test_nunber)))<br/>

BIOmotifvsnontrain <- rbind(Pdata_F[Train_nunber,],Ndata_F[Train_nunber,])<br/>
BIOmotifvsnontest <- rbind(Pdata_F[test_nunber,],Ndata_F[test_nunber,])<br/>

##training model<br/>
BIOmotifvsnon_SVM <- svm(BIOmotifvsnontrain,label_train,cross=5, probability = TRUE)<br/>

pred <- predict(BIOmotifvsnon_SVM,BIOmotifvsnontest, probability = TRUE)<br/>
motif <- attr(pred, 'probabilities')[,1]<br/>
testppred <- prediction(motif,testlabel)<br/>
testpppred_auc <- performance(testppred,"auc")<br/>

#performance<br/>
attr(testpppred_auc,"y.values")[[1]]
