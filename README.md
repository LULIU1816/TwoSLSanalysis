# TwoSLSanalysis

TwoSLSanalysis is an  two-stage least squares procedure for TCGA analysis. TwoSLSanalysis can provide the estimate of causal effect with the survival data download from TCGA.


# Installation

It is easy to install the development version of TwoSLSanalysis package using the 'devtools' package. The typical install time on a "normal" desktop computer is less than one minute.

```
# install.packages("devtools")
library(devtools)
install_github("LULIU1816/TwoSLSanalysis")
```


# Usage

The main function is *TwoSLSanalysis* for two stage analysis of cancer (eg.Lung adenocarcinoma (LUAD)) individual level data downloaded from TCGA, other functions *casecontrolfunction* which has a 1:1 case-control gene expression data; *edgeRfunction* which performs a differentially expressed genes analysis required edgeR package; *ChAMPfunction* which performs a differentially methylated genes analysis required ChAMP package; *selectmethysample*, *selectexpressionsample* which selects the tumor and normal tissues in methylation and expression respectively; *overlapaverage* which averages the repeated measures data of LUAD methylation if exisit, and concrete analysis should be made according to concrete circumstance; *clinicalfunction* which performs clinical data such as transforming one type of data to another; *stage1function* which obtains predictions of the first stage with the advantage of fast computation and memory release.  Moreover, you can find the instructions by '?TwoSLSanalysis','?casecontrolfunction' and so on. 
```R
library(TwoSLSanaylsis)

?TwoSLSanalysis

?casecontrolfunction

?edgeRfunction

?ChAMPfunction

?selectmethysample

?selectexpressionsample

?overlapaverage

?clinicalfunction

?stage1function
```

# Example
```

#load required R packages

library("TwoSLSanalysis")

library("data.table")

library("edgeR")

library("ChAMP")

library("survival")

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

result <-TwoSLSanalysis(clinical, data_normalized, data_count, data_methy, testDir) 


# Development
This R package is developed by Lu Liu.
```
