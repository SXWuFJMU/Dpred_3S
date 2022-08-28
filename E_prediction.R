library(BSgenome)
library(dplyr)
library(e1071)
library(dplyr)
library(caret)
library(ROCR)

Ecoli_D_N <- readRDS("/FA_RDS/Ecoli_D_Negative_test.rds")
Ecoli_D_P <- readRDS("/FA_RDS/Ecoli_D_Positive_test.rds")

E_order_COMPOSI <- readRDS("/F1_order/E_order_COMPOSI.rds")
E_order_Frequency <- readRDS("/F1_order/E_order_Frequency.rds")
E_order_EIIP <- readRDS("/F1_order/E_order_EIIP.rds")

source("/methods/composition.R")
source("/methods/Frequency.R")
source("/methods/EIIP.R")
Ndata_F1 <- list()
Pdata_F1 <- list()
Pdata_F1[[1]] <- composition(Ecoli_D_P)[,E_order_COMPOSI]
Pdata_F1[[2]] <- Frequency(Ecoli_D_P)[,E_order_Frequency]
Pdata_F1[[3]] <- electronIonInteraction(gsub("U","T",Ecoli_D_P))[,E_order_EIIP]
Ndata_F1[[1]] <- composition(Ecoli_D_N)[,E_order_COMPOSI]
Ndata_F1[[2]] <- Frequency(Ecoli_D_N)[,E_order_Frequency]
Ndata_F1[[3]] <- electronIonInteraction(gsub("U","T",Ecoli_D_N))[,E_order_EIIP]

Pdata_F <- cbind(Pdata_F1[[1]],Pdata_F1[[2]],Pdata_F1[[3]])
Ndata_F <- cbind(Ndata_F1[[1]],Ndata_F1[[2]],Ndata_F1[[3]])
BIOmotifvsnontest <- rbind(Pdata_F,Ndata_F)
testlabel <- c(rep(1,nrow(Pdata_F1[[1]])),rep(0,nrow(Pdata_F1[[1]])))

E_final_model <- readRDS("/model/E_final_model.rds")
pred <- predict(E_final_model,BIOmotifvsnontest, probability = TRUE)
motif <- attr(pred, 'probabilities')[,1]
testppred <- prediction(motif,testlabel)
testpppred_auc <- performance(testppred,"auc")

attr(testpppred_auc,"y.values")[[1]]
