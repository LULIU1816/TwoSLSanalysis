selectmethysample <- function(methylationLUAD){
  index <- substr(colnames(methylationLUAD), 14, 16)
  methylationLUAD1 <- methylationLUAD[, index == "01A"]
  methylationLUAD2 <- methylationLUAD[, index == "11A"]
  methylationLUAD <- cbind(methylationLUAD1, methylationLUAD2)
  return(methylationLUAD)
}