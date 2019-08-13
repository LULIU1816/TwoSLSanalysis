# TwoSLSanylsis

TwoSLSanylsis is an  two-stage least squares procedure for TCGA analysis. TwoSLSanylsis can provide the estimate of causal effect with the survival data download from TCGA.


# Installation

It is easy to install the development version of TwoSLSanylsis package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("LULU1816/TwoSLSanylsis")
```


# Usage

The main functions is *TwoSLSanylsisLUAD* for Lung adenocarcinoma (LUAD) individual level data, and *TwoSLSanylsisLUSC* for Lung squamous cell carcinoma (LUSC) individual level data. The difference between *TwoSLSanylsisLUAD* and *TwoSLSanylsisLUSC* is the exsitance of repeated measures data of LUAD methylation. Thus, concrete analysis should be made according to concrete circumstance. Moreover,  *TwoSLSanylsisLUSC* is the general function in other cancer analysis. You can find the instructions by '?TwoSLSanylsisLUAD' and '?TwoSLSanylsisLUSC'. 

library(TwoSLSanylsis)

?TwoSLSanylsisLUAD

?TwoSLSanylsisLUSC

# Example

#read clinical data
clinical<-fread("nationwidechildrens.org_clinical_patient_lusc.txt",data.table=F,head=T)
clinical <- as.data.frame(clinical)[-c(1,2),]

#read expression data
data_normalized<-fread("LUSC__gene.normalized_RNAseq__tissueTypeAll__20181022104424.txt",data.table=F,head=T)
data_count<-fread("LUSC__gene_RNAseq__tissueTypeAll__20181022105257.txt",data.table=F,head=T)
data_normalized <- as.data.frame(data_normalized)
data_count <-as.data.frame(data_count)

#read methylation data
data_methy<-fread("LUSC__methylation_450__tissueTypeAll__20181024055604.txt",data.table=F,head=T)
data_methy <- as.data.frame(data_methy)[,-c(2,3,4)]

#configure path of GDC data
testDir = "~/GDCdata"

LUSCresult <-TwoSLSanalysisLUSC(clinical,data_normalized,data_count,data_methy,testDir) 


# Development
This R package is developed by Lu Liu.
