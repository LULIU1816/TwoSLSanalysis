selectexpressionsample <- function(gene_id,expressionLUAD){
  index <- substr(colnames(expressionLUAD), 14, 16)
  expressionLUAD1 <- expressionLUAD[, index == "01A"]
  expressionLUAD2 <- expressionLUAD[, index == "11A"]
  expressionLUAD <- cbind(expressionLUAD1, expressionLUAD2)
  expressionLUAD <- cbind(gene_id, expressionLUAD)
  return(expressionLUAD)
}