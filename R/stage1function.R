  stage1function <- function(i,XX,commonID,methylationLUAD,myDMPLUADTSS1500,expressionLUAD,clinical2){
    # 
    # TSS1500 <- list()
    # for (i in 1:length(XX)) {
    #   TSS1500[[i]] <- NULL
    # }
    # 
    n <- length(methylationLUAD)
    # methylationLUADTSS1500 <- list()
    # for (i in 1:length(XX)) {
    methylationLUADTSS1500 <- NaN * seq(n)
    # }
    # 
    # expressionTSS1500 <- list()
    # for (i in 1:length(XX)) {
    #   expressionTSS1500[[i]] <- NaN * seq(n)
    # }
    # w <- n-1
    # expressionTSS15003 <- list()
    # for (i in 1:length(XX)) {
    #   expressionTSS15003[[i]] <- NaN * seq(w)
    # }
    # 
    # methylationTSS15003 <- list()
    # for (i in 1:length(XX)) {
    #   methylationTSS15003[[i]] <- NaN * seq(w)
    # }
    # 
    # e <- w
    # X3 <- list()
    # for (i in 1:length(XX)) {
    #   X3[[i]] <- NaN * seq(e)
    # }
    # 
    # ZTSS1500 <- list()
    # for (i in 1:length(XX)) {
    #   ZTSS1500[[i]] <- NaN * seq(e)
    # }
    # 
    # XhTSS1500 <- list()
    # for (i in 1:length(XX)) {
    #   XhTSS1500[[i]] <- NaN
    # }
    
    TSS1500 <- as.character(myDMPLUADTSS1500$V1[myDMPLUADTSS1500$gene == XX[i]])
    for (j in 1:(length(TSS1500))) {
      methylationLUADTSS1500 <- rbind(methylationLUADTSS1500, methylationLUAD[methylationLUAD[, 1] == TSS1500[j], ])
    }
    methylationLUADTSS1500<- methylationLUADTSS1500[-1, ]
    
    expressionTSS1500 <- expressionLUAD[expressionLUAD[, 1] == XX[i], ]
    
    expressionTSS15003 <- expressionTSS1500[expressionLUAD[1,] %in% commonID]

    methylationTSS15003 <- methylationLUADTSS1500[expressionLUAD[1,] %in% commonID]
    
    X3 <- t(expressionTSS15003[1, ])
    X3 <- as.matrix(X3)
    X3 <- apply(X3, 2, as.numeric)
    
    ZTSS1500 <- t(methylationTSS15003)
    ZTSS1500<- data.frame(1, ZTSS1500)
    ZTSS1500 <- as.matrix(ZTSS1500)
    ZTSS1500 <- apply(ZTSS1500, 2, as.numeric)
    
    mod <- lm(X3 ~ clinical2$gender + clinical2$age_at_initial_pathologic_diagnosis + clinical2$tobacco_smoking_pack_years_smoked)
    e <- as.matrix(residuals(mod))
    ZTSS1500 <- t(na.omit(t(ZTSS1500)))
    if (det(t(ZTSS1500) %*% ZTSS1500) != 0) {
      beta <- solve(t(ZTSS1500) %*% ZTSS1500) %*% t(ZTSS1500) %*% e
      XhTSS1500 <- ZTSS1500 %*% beta
    } else {
      XhTSS1500 <- NA
    }
  return(XhTSS1500)
}
