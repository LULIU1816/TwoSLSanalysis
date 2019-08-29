  stage1function <- function(i,XX,commonID,methylationLUAD,myDMPLUADTSS1500,expressionLUAD,clinical2){
     
    TSS1500[[i]] <- as.character(myDMPLUADTSS1500$V1[myDMPLUADTSS1500$gene == XX[i]])
    for (j in 1:(length(TSS1500[[i]]))) {
      methylationLUADTSS1500[[i]] <- rbind(methylationLUADTSS1500[[i]], methylationLUAD[methylationLUAD[, 1] == TSS1500[[i]][j], ])
    }
    methylationLUADTSS1500[[i]] <- methylationLUADTSS1500[[i]][-1, ]
    
    expressionTSS1500[[i]] <- expressionLUAD[expressionLUAD[, 1] == XX[i], ]
    
    expressionTSS15003[[i]] <- expressionTSS1500[[i]][expressionLUAD[1,] %in% commonID]

    methylationTSS15003[[i]] <- methylationLUADTSS1500[[i]][expressionLUAD[1,] %in% commonID]
    
    X3[[i]] <- t(expressionTSS15003[[i]][1, ])
    X3[[i]] <- as.matrix(X3[[i]])
    X3[[i]] <- apply(X3[[i]], 2, as.numeric)
    
    ZTSS1500[[i]] <- t(methylationTSS15003[[i]])
    ZTSS1500[[i]] <- data.frame(1, ZTSS1500[[i]])
    ZTSS1500[[i]] <- as.matrix(ZTSS1500[[i]])
    ZTSS1500[[i]] <- apply(ZTSS1500[[i]], 2, as.numeric)
    
    mod <- lm(X3[[i]] ~ clinical2$gender + clinical2$age_at_initial_pathologic_diagnosis + clinical2$tobacco_smoking_pack_years_smoked)
    e <- as.matrix(residuals(mod))
    ZTSS1500[[i]] <- t(na.omit(t(ZTSS1500[[i]])))
    if (det(t(ZTSS1500[[i]]) %*% ZTSS1500[[i]]) != 0) {
      beta <- solve(t(ZTSS1500[[i]]) %*% ZTSS1500[[i]]) %*% t(ZTSS1500[[i]]) %*% e
      XhTSS1500[[i]] <- ZTSS1500[[i]] %*% beta
    } else {
      XhTSS1500[[i]] <- NA
    }
  return(XhTSS1500[[i]])
}
