stage1methy <- function(XX,methylationLUAD,myDMPLUADTSS1500){
  TSS1500 <- list()
  for (i in 1:length(XX)) {
    TSS1500[[i]] <- NULL
  }
  
  n <- length(methylationLUAD)
  methylationLUADTSS1500 <- list()
  for (i in 1:length(XX)) {
    methylationLUADTSS1500[[i]] <- NaN * seq(n)
  }
  
  for (i in 1:length(XX)) {
    TSS1500[[i]] <- as.character(myDMPLUADTSS1500$V1[myDMPLUADTSS1500$gene == XX[i]])
    for (j in 1:(length(TSS1500[[i]]))) {
      methylationLUADTSS1500[[i]] <- rbind(methylationLUADTSS1500[[i]], methylationLUAD[methylationLUAD[, 1] == TSS1500[[i]][j], ])
    }
    methylationLUADTSS1500[[i]] <- methylationLUADTSS1500[[i]][-1, ]
  }
  return(methylationLUADTSS1500)
}