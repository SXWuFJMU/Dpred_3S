library(BSgenome)

NC1 <- function(data){
  for(i in 1:length(data) ){
    if(i==1){
      DCfirst <- unlist(as.vector(strsplit(data[1],"",fixed = TRUE)))
      DCsecond <- matrix(NA,nrow = length(data),ncol = length(DCfirst))
      DCsecond[1,] <-  DCfirst 
    }else{
      DCsecond[i,] <- unlist(as.vector(strsplit(data[i],"",fixed = TRUE)))
    }
  }
  return(DCsecond)
}

library(matrixStats)

Frequency <- function(Data,NTYPE="RNA") {
  sequences_M <- NC1(Data)
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  N = ncol(sequences_M)
  
  cumFreq_A <- rowCumsums(matrix(as.numeric( sequences_M == "A" ), ncol = N, byrow = F))
  cumFreq_T <- rowCumsums(matrix(as.numeric( sequences_M ==  U ), ncol = N, byrow = F))
  cumFreq_C <- rowCumsums(matrix(as.numeric( sequences_M == "C" ), ncol = N, byrow = F))
  cumFreq_G <- rowCumsums(matrix(as.numeric( sequences_M == "G" ), ncol = N, byrow = F))
  
  cumFreq_combined <- matrix(0,ncol = N, nrow = length(Data))
  cumFreq_combined[sequences_M == "A"] <- cumFreq_A[sequences_M == "A"]
  cumFreq_combined[sequences_M == U] <- cumFreq_T[sequences_M == U]
  cumFreq_combined[sequences_M == "C"] <- cumFreq_C[sequences_M == "C"]
  cumFreq_combined[sequences_M == "G"] <- cumFreq_G[sequences_M == "G"]
  
  cumFreq_combined <- t(t(cumFreq_combined) / seq_len(N))
  colnames(cumFreq_combined ) <-  paste0("cumFreq_",seq_len(N))
  return(cumFreq_combined)
}

