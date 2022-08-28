ChemicalProperty <- function(Data){
  MA <- NC1(Data)
  MA2 <- matrix(NA,ncol =(3*ncol(MA)),nrow = nrow(MA))
  A <- c(1,1,1)%>%matrix(nrow=1)
  C <- c(0,1,0)%>%matrix(nrow=1)
  G <- c(1,0,0)%>%matrix(nrow=1)
  T <- c(0,0,1)%>%matrix(nrow=1)
  U <- c(0,0,1)%>%matrix(nrow=1)
  N <- c(0,0,0)%>%matrix(nrow=1)
  for(i in 1:ncol(MA)){
    i1 <- (i-1)*3 +1
    i2 <- (i-1)*3 +2
    i3 <- i*3
    for(j in c("A","C","G","T","U")){
      if(length(which(MA[,i]==j))>0){
        MA2[which(MA[,i]==j),i1] <- eval(parse(text=j))[1]
        MA2[which(MA[,i]==j),i2] <- eval(parse(text=j))[2]
        MA2[which(MA[,i]==j),i3] <- eval(parse(text=j))[3]
      }
      # print(eval(parse(text=j)))
    }
  }
  return(MA2)
}