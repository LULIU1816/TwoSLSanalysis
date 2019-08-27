#' TwoSLSanalysisLUSC
#' @param clinical input data
#' @param data_normalized input data
#' @param data_count input data
#' @param data_methy input data
#' @param testDir input teatDir
#'
#' @return data.frame:XX, the candidate genes; hTSS1500, beta value of 2SLS; ORTSS1500, HR of cox model; lower.95, upper.95: 95%CI; PTSS1500, P-valve; FDRH, FDR; test.ph, cox PH-test.
#'
#' @export
#'
TwoSLSanalysisLUSC <- function(clinical, data_normalized, data_count, data_methy, testDir) {

    # structure case-control data
    case_control_data <- casecontrolfunction(data_count,data_tumor)    
    # edgeR
    diff_gene_edgeR <- edgeRfunction(case_control_data)
    diff_gene <- diff_gene_edgeR[,1]
    # CHAMP
    myDMPLUAD <- ChAMPfunction(testDir)
    
    # select CpGs in TSS1500
    myDMPLUAD <- myDMPLUAD[myDMPLUAD$gene != "", ]
    myDMPLUADTSS1500 <- myDMPLUAD[myDMPLUAD$feature == "TSS1500", ]
    # corresponding genes(DMG)
    DMPgeneLUADTSS1500 <- myDMPLUADTSS1500$gene
    DMPgeneLUADTSS1500 <- unique(DMPgeneLUADTSS1500)
    # overlapping the DEG and DMG
    LUADgeneTSS1500OVERLAP <- intersect(DMPgeneLUADTSS1500, diff_gene)

    # sample matching between expression and methylation
    sampleID_methylationLUAD <- colnames(data_methy)[-1]
    sampleID_expressionLUAD <- colnames(data_normalized)[-1]
    ID_methylationLUAD <- substr(sampleID_methylationLUAD, 1, 16)
    ID_expressionLUAD <- substr(sampleID_expressionLUAD, 1, 16)
    commonID <- ID_expressionLUAD[ID_expressionLUAD %in% ID_methylationLUAD]
    ID_methylationLUAD0 <- c("ID_methylationLUAD", ID_methylationLUAD)
    ID_expressionLUAD0 <- c("ID_expressionLUAD", ID_expressionLUAD)
    methylationLUAD0 <- rbind(ID_methylationLUAD0, data_methy)
    expressionLUAD0 <- rbind(ID_expressionLUAD0, data_normalized)
    gene_id <- expressionLUAD0[, 1]
    CpG <- methylationLUAD0[, 1]
    methylationLUAD <- methylationLUAD0[, which(methylationLUAD0[1, ] %in% commonID)]
    methylationLUAD <- cbind(CpG, methylationLUAD)
    # select sample 01A,11A 
    methylationLUAD <- selectmethysample(methylationLUAD)

    expressionLUAD <- expressionLUAD0[, which(expressionLUAD0[1, ] %in% commonID)]
    # select sample 01A,11A 
    expressionLUAD <- selectexpressionsample(gene_id,expressionLUAD) 
    # take the average in overlapping samples(concrete analysis should be made according to concrete circumstance)
    if( dim(expressionLUAD)[2]!=dim(methylationLUAD)[2]+1){
      methylationLUAD <- overlapaverage(CpG,methylationLUAD)
    }
    # sample in order
    expressionLUAD <- expressionLUAD[, order(expressionLUAD[1, ])]
    methylationLUAD <- methylationLUAD0[, order(methylationLUAD0[1, ])]
    # stage one
    XX <- LUADgeneTSS1500OVERLAP
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
            methylationLUADTSS1500[[i]] <- rbind(methylationLUADTSS1500[[i]], methylationLUAD[methylationLUAD[, 1] == TSS1500[[i]][j],
                ])
        }
        methylationLUADTSS1500[[i]] <- methylationLUADTSS1500[[i]][-1, ]
    }
    m <- length(LUADexpression[1, ])
    expressionTSS1500 <- list()
    for (i in 1:length(XX)) {
        expressionTSS1500[[i]] <- NaN * seq(m)
    }
    for (i in 1:length(XX)) {
        expressionTSS1500[[i]] <- LUADexpression[LUADexpression[, 1] == XX[i], ]
    }

    commonID2 <- LUADexpression[2, -1]
    commonID2 <- t(commonID2)
    commonID2 <- c(commonID2)
    commonID22 <- commonID2[substr(commonID2, 14, 16) == "01A"]
    clinical2 <- clinicalfunction(clinical,commonID22)
    commonID <- clinical2$bcr_patient_barcode

    w <- length(expressionTSS1500[[1]][, -1])
    expressionTSS15003 <- list()
    for (i in 1:length(XX)) {
        expressionTSS15003[[i]] <- NaN * seq(w)
    }
    for (i in 1:length(XX)) {
        expressionTSS15003[[i]] <- expressionTSS1500[[i]][, -1]
        expressionTSS15003[[i]] <- expressionTSS15003[[i]][LUADexpression[2, -1] %in% commonID]
    }
    # methylation对应
    methylationTSS15003 <- list()
    for (i in 1:length(XX)) {
        methylationTSS15003[[i]] <- NaN * seq(w)
    }

    for (i in 1:length(XX)) {
        methylationTSS15003[[i]] <- methylationTSS1500[[i]][, -1]
        methylationTSS15003[[i]] <- methylationTSS15003[[i]][LUADexpression[2, -1] %in% commonID]
    }
    e <- length(expressionTSS15003[[1]][1, ])
    X3 <- list()
    for (i in 1:length(XX)) {
        X3[[i]] <- NaN * seq(e)
    }
    for (i in 1:length(XX)) {
        X3[[i]] <- t(expressionTSS15003[[i]][1, ])
        X3[[i]] <- as.matrix(X3[[i]])
        X3[[i]] <- apply(X3[[i]], 2, as.numeric)
    }
    ZTSS1500 <- list()
    for (i in 1:length(XX)) {
        ZTSS1500[[i]] <- NaN * seq(e)
    }
    for (i in 1:length(XX)) {
        ZTSS1500[[i]] <- t(methylationTSS15003[[i]])
        ZTSS1500[[i]] <- data.frame(rep(1, e), ZTSS1500[[i]])
        ZTSS1500[[i]] <- as.matrix(ZTSS1500[[i]])
        ZTSS1500[[i]] <- apply(ZTSS1500[[i]], 2, as.numeric)
    }

    XhTSS1500 <- list()
    for (i in 1:length(XX)) {
        XhTSS1500[[i]] <- NaN
    }
    for (i in 1:length(XX)) {
        mod <- lm(X3[[i]] ~ clinical2$gender + clinical2$age_at_initial_pathologic_diagnosis + clinical2$tobacco_smoking_pack_years_smoked)
        e <- as.matrix(residuals(mod))
        ZTSS1500[[i]] <- t(na.omit(t(ZTSS1500[[i]])))
        if (det(t(ZTSS1500[[i]]) %*% ZTSS1500[[i]]) != 0) {
            beta <- solve(t(ZTSS1500[[i]]) %*% ZTSS1500[[i]]) %*% t(ZTSS1500[[i]]) %*% e
            XhTSS1500[[i]] <- ZTSS1500[[i]] %*% beta
        } else {
            XhTSS1500[[i]] <- NA
        }
    }
    # stage two
    s <- Surv(time = clinical2$death_days_to, event = clinical2$vital_status)
    PTSS1500 <- NULL
    ORTSS1500 <- NULL
    hTSS1500 <- NULL
    lower.95 <- NULL
    upper.95 <- NULL
    test.ph <- NULL
    for (i in 1:length(XX)) {
        coxph.fit <- coxph(s ~ XhTSS1500[[i]] + clinical2$gender + clinical2$age_at_initial_pathologic_diagnosis + clinical2$tobacco_smoking_pack_years_smoked)
        PTSS1500[i] <- summary(coxph.fit)$coefficients[1, 5]
        ORTSS1500[i] <- summary(coxph.fit)$coefficients[1, 2]
        hTSS1500[i] <- summary(coxph.fit)$coefficients[1, 1]
        lower.95[i] <- summary(coxph.fit)$conf.int[1, 3]
        upper.95[i] <- summary(coxph.fit)$conf.int[1, 4]
        test.ph[i] <- cox.zph(coxph.fit)$table[4, 3]
    }
    H <- data.frame(XX, hTSS1500, ORTSS1500, PTSS1500, lower.95, upper.95)
    H <- H[!is.na(H$PTSS1500), ]
    FDRH <- p.adjust(H$PTSS1500, method = "fdr", n = length(H$PTSS1500))
    return(data.frame(XX, hTSS1500, ORTSS1500, lower.95, upper.95, PTSS1500, FDRH, test.ph))
}
