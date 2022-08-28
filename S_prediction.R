library(BSgenome)
library(dplyr)
library(e1071)
library(dplyr)
library(caret)
library(ROCR)

SpombS_D_N <- readRDS("/FA_RDS/Spombe_D_Negative_test.rds")
SpombS_D_P <- readRDS("/FA_RDS/Spombe_D_Positive_test.rds")

S_order_COMPOSI <- readRDS("/F1_order/S_order_COMPOSI.rds")
S_order_ChemicalProper <- readRDS("/F1_order/S_order_ChemicalProper.rds")
S_order_PseKNC <- readRDS("/F1_order/S_order_PseKNC.rds")

source("/methods/composition.R")
source("/methods/PseKNC.R")
source("/methods/ChemicialProperty.R")

Ndata_F1 <- list()
Pdata_F1 <- list()

Pdata_F1[[1]] <- composition(SpombS_D_P)[,S_order_COMPOSI]
Pdata_F1[[2]] <- ChemicalProperty(gsub("U","T",SpombS_D_P))[,S_order_ChemicalProper]
Pdata_F1[[3]] <- PseKNC(gsub("U","T",SpombS_D_P))[,S_order_PseKNC]
Ndata_F1[[1]] <- composition(SpombS_D_N)[,S_order_COMPOSI]
Ndata_F1[[2]] <- ChemicalProperty(gsub("U","T",SpombS_D_N))[,S_order_ChemicalProper]
Ndata_F1[[3]] <- PseKNC(gsub("U","T",SpombS_D_N))[,S_order_PseKNC]

Pdata_F <- cbind(Pdata_F1[[1]],Pdata_F1[[2]],Pdata_F1[[3]])
Ndata_F <- cbind(Ndata_F1[[1]],Ndata_F1[[2]],Ndata_F1[[3]])
BIOmotifvsnontest <- rbind(Pdata_F,Ndata_F)
testlabel <- c(rep(1,nrow(Pdata_F1[[1]])),rep(0,nrow(Pdata_F1[[1]])))

S_final_model <- readRDS("/model/S_final_model.rds")
pred <- predict(S_final_model,BIOmotifvsnontest, probability = TRUE)
motif <- attr(pred, 'probabilities')[,1]
testppred <- prediction(motif,testlabel)
testpppred_auc <- performance(testppred,"auc")

attr(testpppred_auc,"y.values")[[1]]
