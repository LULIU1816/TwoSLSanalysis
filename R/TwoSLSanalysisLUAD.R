#' TwoSLSanalysisLUAD
#' @param clinical input data
#' @param data_normalized input data
#' @param data_count input data
#' @param data_methy input data
#' @param testDir input teatDir
#'
#' @return data.frame:XX, the candidate genes; hTSS1500, beta value of 2SLS; ORTSS1500, HR of cox model; lower.95, upper.95: 95%CI; PTSS1500, P-valve; FDRH, FDR; test.ph, cox PH-test.
#'
#' @export

TwoSLSanalysisLUAD <- function(clinical, data_normalized, data_count, data_methy, testDir) {

    clinical$gender <- ifelse(clinical$gender == "MALE", 1, 2)
    clinical$vital_status <- ifelse(clinical$vital_status == "Alive", 1, 2)
    clinical$tobacco_smoking_pack_years_smoked <- ifelse(clinical$tobacco_smoking_pack_years_smoked == "[Not Available]", NA, clinical$tobacco_smoking_pack_years_smoked)
    clinical$tobacco_smoking_pack_years_smoked <- as.numeric(clinical$tobacco_smoking_pack_years_smoked)
    for (i in 1:length(clinical$death_days_to)) {
        clinical$death_days_to[i] <- ifelse(clinical$death_days_to[i] == "[Not Available]", clinical$last_contact_days_to[i], clinical$death_days_to[i])
        clinical$death_days_to[i] <- ifelse(clinical$death_days_to[i] == "[Not Applicable]", clinical$last_contact_days_to[i], clinical$death_days_to[i])
    }
    clinical$death_days_to <- ifelse(clinical$death_days_to == "[Not Available]", NA, clinical$death_days_to)
    clinical$death_days_to <- ifelse(clinical$death_days_to == "[Discrepancy]", NA, clinical$death_days_to)
    clinical$death_days_to <- as.numeric(clinical$death_days_to)
    clinical1 <- clinical[!is.na(clinical$death_days_to), ]
    clinical1$bcr_patient_barcode <- paste(clinical1$bcr_patient_barcode, "01A", sep = "-")
    clinical <- clinical1

    # edgeR
    data_count <- data_count[, data_count[1, ] != "scaled_estimate"][-1, ]
    sampleID_count <- colnames(data_count)[-1]
    # leave sample 11B out
    index <- substr(sampleID_count, 14, 16)
    data_control <- data_count[, -1][, index == "11A"]
    data_tumor <- data_count[, -1][, index == "01A"]
    gene_id <- data_count[, 1]
    data_control <- cbind(gene_id, data_control)
    data_tumor <- cbind(gene_id, data_tumor)
    # leave gene ? out
    df <- data_control[, 1]
    n <- grep("\\?", df)
    data_control <- data_control[-n, ]
    dg <- data_tumor[, 1]
    m <- grep("\\?", dg)
    data_tumor <- data_tumor[-m, ]
    # sample matching
    sampleID_control <- colnames(data_control)[-1]
    sampleID_tumor <- colnames(data_tumor)[-1]
    ID_gene_normal <- substr(sampleID_control, 1, 12)
    ID_gene_tumor <- substr(sampleID_tumor, 1, 12)
    commonID <- ID_gene_tumor[ID_gene_tumor %in% ID_gene_normal]
    ID_gene_normal0 <- c("ID_gene_normal", ID_gene_normal)
    ID_gene_tumor0 <- c("ID_gene_tumor", ID_gene_tumor)
    data_tumor <- rbind(ID_gene_tumor0, data_tumor)
    data_control <- rbind(ID_gene_normal0, data_control)
    gene_id <- data_tumor[, 1]
    data_case0 <- data_tumor[, which(data_tumor[1, ] %in% commonID)]
    data_case1 <- cbind(gene_id, data_case0)
    gene_id <- data_control[, 1]
    data_control0 <- data_control[, which(data_control[1, ] %in% commonID)]
    data_control1 <- cbind(gene_id, data_control0)
    # merge the case-control table
    case_control_data <- merge(data_case1[-1, ], data_control1[-1, ], by = "gene_id")
    case_control_data[, 2:ncol(case_control_data)] <- lapply(case_control_data[, 2:ncol(case_control_data)], as.numeric)
    group <- factor(c(rep("H", (ncol(case_control_data) - 1)/2), rep("M", (ncol(case_control_data) - 1)/2)))
    # remove genes with low expression
    case_control_data1 <- case_control_data[rowSums(cpm(case_control_data[, 2:ncol(case_control_data)]) > 1) >= 2, ]
    c_c_data <- DGEList(counts = case_control_data1[, 2:ncol(case_control_data)], genes = case_control_data1[, 1], group = group)
    ## flitration and normalization
    c_c_data <- calcNormFactors(c_c_data)
    ## design matrix
    design <- model.matrix(~group)
    ## dispersion
    exprSet <- estimateCommonDisp(c_c_data, design)
    exprSet <- estimateTagwiseDisp(exprSet)
    ## different expression gene(DEG)
    et <- exactTest(exprSet)
    tTag <- topTags(et, n = nrow(exprSet))
    diff_gene_edgeR <- subset(tTag$table, FDR < 0.001 & (logFC > 2.5 | logFC < -2.5) & (logCPM > 2.5 | logCPM < -2.5))
    diff_gene <- diff_gene_edgeR[, 1]
    w <- gregexpr("\\|", diff_gene)
    p <- NULL
    for (i in 1:length(diff_gene)) {
        p[i] <- substr(diff_gene[i], 1, w[[i]][1] - 1)
    }
    diff_gene <- p
    diff_gene_edgeR[, 1] <- diff_gene

    # CHAMP
    myLoad <- champ.load(testDir, arraytype = "450K")
    # distrubution of CpGs CpG.GUI() QC QC.GUI(arraytype = '450K') champ.QC() Normalization
    myNorm <- champ.norm(beta = myLoad$beta, arraytype = "450K", cores = 5)
    # save write.csv(myNorm,file='./LUADNormalization Data.csv',quote=F,row.names = T) SVD champ.SVD(beta=myNorm,pd=myLoad$pd)
    # different methylaton position
    myDMP <- champ.DMP(beta = myNorm, pheno = myLoad$pd$Sample_Group)
    # DMP.GUI()
    myDMPLUAD <- myDMP[[1]]
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
    # leave sample 01B etc. out
    index <- substr(colnames(methylationLUAD), 14, 16)
    methylationLUAD1 <- methylationLUAD[, index == "01A"]
    methylationLUAD2 <- methylationLUAD[, index == "11A"]
    methylationLUAD <- cbind(methylationLUAD1, methylationLUAD2)
    # take the average in overlapping samples(concrete analysis should be made according to concrete circumstance)
    a <- as.numeric(methylationLUAD$`TCGA-44-5645-01A-01D-1626-05`)
    a <- ifelse(is.na(a), as.numeric(methylationLUAD$`TCGA-44-5645-01A-01D-A276-05`), a)
    b <- as.numeric(methylationLUAD$`TCGA-44-6146-01A-11D-1756-05`)
    b <- ifelse(is.na(b), as.numeric(methylationLUAD$`TCGA-44-6146-01A-11D-A276-05`), b)
    c <- as.numeric(methylationLUAD$`TCGA-44-6147-01A-11D-1756-05`)
    c <- ifelse(is.na(c), as.numeric(methylationLUAD$`TCGA-44-6147-01A-11D-A276-05`), c)
    d <- as.numeric(methylationLUAD$`TCGA-44-6775-01A-11D-1856-05`)
    d <- ifelse(is.na(d), as.numeric(methylationLUAD$`TCGA-44-6775-01A-11D-A276-05`), d)
    a <- (as.numeric(methylationLUAD$`TCGA-44-5645-01A-01D-1626-05`) + as.numeric(methylationLUAD$`TCGA-44-5645-01A-01D-A276-05`))/2
    b <- (as.numeric(methylationLUAD$`TCGA-44-6146-01A-11D-1756-05`) + as.numeric(methylationLUAD$`TCGA-44-6146-01A-11D-A276-05`))/2
    c <- (as.numeric(methylationLUAD$`TCGA-44-6147-01A-11D-1756-05`) + as.numeric(methylationLUAD$`TCGA-44-6147-01A-11D-A276-05`))/2
    d <- (as.numeric(methylationLUAD$`TCGA-44-6775-01A-11D-1856-05`) + as.numeric(methylationLUAD$`TCGA-44-6775-01A-11D-A276-05`))/2
    mydf <- subset(methylationLUAD, select = -c(`TCGA-44-5645-01A-01D-1626-05`, `TCGA-44-5645-01A-01D-A276-05`, `TCGA-44-6146-01A-11D-1756-05`,
        `TCGA-44-6146-01A-11D-A276-05`, `TCGA-44-6147-01A-11D-1756-05`, `TCGA-44-6147-01A-11D-A276-05`, `TCGA-44-6775-01A-11D-1856-05`,
        `TCGA-44-6775-01A-11D-A276-05`))
    methylationLUAD0 <- cbind(mydf, a)
    methylationLUAD0 <- cbind(methylationLUAD0, b)
    methylationLUAD0 <- cbind(methylationLUAD0, c)
    methylationLUAD0 <- cbind(methylationLUAD0, d)
    methylationLUAD0 <- cbind(CpG, methylationLUAD0)
    colnames(methylationLUAD0)[470:473] <- c("TCGA-44-5645-01A", "TCGA-44-6146-01A", "TCGA-44-6147-01A", "TCGA-44-6775-01A")
    methylationLUAD0[1, 470:473] <- c("TCGA-44-5645-01A", "TCGA-44-6146-01A", "TCGA-44-6147-01A", "TCGA-44-6775-01A")
    expressionLUAD <- expressionLUAD0[, which(expressionLUAD0[1, ] %in% commonID)]
    # leave sample 01B out
    index <- substr(colnames(expressionLUAD), 14, 16)
    expressionLUAD1 <- expressionLUAD[, index == "01A"]
    expressionLUAD2 <- expressionLUAD[, index == "11A"]
    expressionLUAD <- cbind(expressionLUAD1, expressionLUAD2)
    expressionLUAD <- cbind(gene_id, expressionLUAD)
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
    clinical0 <- clinical[clinical$bcr_patient_barcode %in% commonID22, ]
    commonID222 <- clinical0[, 3]
    clinical0$tobacco_smoking_pack_years_smoked <- ifelse(clinical0$tobacco_smoking_pack_years_smoked == "[Not Available]", NA, clinical0$tobacco_smoking_pack_years_smoked)
    clinical1 <- clinical0[!is.na(clinical0$tobacco_smoking_pack_years_smoked), ]
    clinical1 <- clinical1[clinical1$death_days_to != 0, ]
    clinical1$age_at_initial_pathologic_diagnosis <- as.numeric(clinical1$age_at_initial_pathologic_diagnosis)
    clinical1$followupage <- (-as.numeric(as.character(clinical1$birth_days_to)))/365
    clinical1$followupage <- floor(clinical1$followupage)
    clinical1$age <- ifelse(clinical1$followupage != clinical1$age_at_initial_pathologic_diagnosis, clinical1$followupage, clinical1$age_at_initial_pathologic_diagnosis)
    clinical1$age_at_initial_pathologic_diagnosis <- ifelse(is.na(clinical1$age), clinical1$age_at_initial_pathologic_diagnosis,
        clinical1$age)
    clinical2 <- clinical1[!is.na(clinical1$age_at_initial_pathologic_diagnosis), ]
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
