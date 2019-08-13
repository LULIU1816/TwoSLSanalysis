\name{TwoSLSanalysisLUSC}
\alias{TwoSLSanalysisLUSC}
\title{TwoSLSanalysis}
\usage{
TwoSLSanalysisLUSC(clinical,data_normalized,data_count,data_methy,testDir)
}
\description{
provides the estimate of causal effect with the survival data in LUAD download from TCGA.
}
\examples{
library("data.table")
library("edgeR")
library("ChAMP")
library("survival")

#read clinical data clinical<-fread("nationwidechildrens.org_clinical_patient_lusc.txt",data.table=F,head=T) clinical <- as.data.frame(clinical)[-c(1,2),]

#read expression data data_normalized<-fread("LUSC__gene.normalized_RNAseq__tissueTypeAll__20181022104424.txt",data.table=F,head=T) data_count<-fread("LUSC__gene_RNAseq__tissueTypeAll__20181022105257.txt",data.table=F,head=T) data_normalized <- as.data.frame(data_normalized) data_count <-as.data.frame(data_count)

#read methylation data data_methy<-fread("LUSC__methylation_450__tissueTypeAll__20181024055604.txt",data.table=F,head=T) data_methy <- as.data.frame(data_methy)[,-c(2,3,4)]

#configure path of GDC data
testDir = "~/GDCdata"

LUSCresult <-TwoSLSanalysisLUSC(clinical,data_normalized,data_count,data_methy,testDir)
}
