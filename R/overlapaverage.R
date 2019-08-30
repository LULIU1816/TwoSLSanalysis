overlapaverage <- function(CpG,methylationLUAD){
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
  return(methylationLUAD0)
  }