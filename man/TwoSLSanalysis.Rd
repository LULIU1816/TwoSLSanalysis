\name{TwoSLSanalysis}
\alias{TwoSLSanalysis}
\title{TwoSLSanalysis}
\description{
TwoSLSanylsis can provide the estimate of causal effect with the survival data download from TCGA.
}
\usage{
TwoSLSanalysisLUAD(clinical,data_normalized,data_count,data_methy,testDir,fdr=0.001,logfc=2.5,logcpm=2.5,cores=1)
}
\examples{
library("data.table")
library("edgeR")
library("ChAMP")
library("survival")

#read clinical data
clinical<-fread("nationwidechildrens.org_clinical_patient_luad.txt",data.table=F,head=T)
clinical <- as.data.frame(clinical)[-c(1,2),]
#read expression data
data_normalized<-fread("LUAD__gene.normalized_RNAseq__tissueTypeAll__20181022065043.txt",data.table=F,head=T)
data_count<-fread("LUAD__gene_RNAseq__tissueTypeAll__20181022065942.txt",data.table=F,head=T)
data_normalized <- as.data.frame(data_normalized)
data_count <-as.data.frame(data_count)
#read methylation data
data_methy<-fread("LUAD__methylation_450__tissueTypeAll__20181024043150.txt",data.table=F,head=T)
data_methy <- as.data.frame(data_methy)[,-c(2,3,4)]

#configure path of GDC data
testDir = "~/GDCdata"

LUADresult <-TwoSLSanalysis(clinical,data_normalized,data_count,data_methy,testDir)
}
