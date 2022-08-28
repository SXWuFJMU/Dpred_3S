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


composition <- function(Data,NI=3,NTYPE="RNA",Freq=2){
  # MA <- NC1(DATA) 
  if(NTYPE=="RNA"){
    U="U"
  }else{
    U="T"
  }
  for(i in 1:NI){
    if(i==1){
      NP <- c("A","G","C",U)
    }else{
      NP1 <- NULL
      for(j in c("A","G","C",U)){
        for(k in 1:length(NP)){
          NP1 <- c(NP1,paste0(NP[k],j))
        }
      }
      NP <- NP1
    }
  }
  MA2 <- matrix(NA,ncol = length(NP),nrow = length(Data))
  colnames(MA2) <- NP
  if(NTYPE=="RNA"){
    for(i in 1:length(NP)){
      MA2[,i]<- vcountPattern(NP[i],RNAStringSet(Data))
    }
  }else{
    for(i in 1:length(NP)){
      MA2[,i]<- vcountPattern(NP[i],DNAStringSet(Data))
    }
  }
  M3 <- MA2
  if(Freq==1){
    for(i in 1:nrow(M3)){
      M3[i,] <-  MA2[i,]/sum(MA2[i,])
    }
  }
  return(M3)
}
