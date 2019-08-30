\name{casecontrolfunction}
\alias{casecontrolfunction}
\title{casecontrolfunction}
\usage{
casecontrolfunction(data_count,data_tumor)
}
\description{
structure 1:1 case-control gene expression data 
}
\examples{
library("data.table")

#read expression data 
data_normalized<-fread("LUSC__gene.normalized_RNAseq__tissueTypeAll__20181022104424.txt",data.table=F,head=T) 
data_count<-fread("LUSC__gene_RNAseq__tissueTypeAll__20181022105257.txt",data.table=F,head=T) 
data_normalized <- as.data.frame(data_normalized) 
data_count <-as.data.frame(data_count)

case_control_data <- casecontrolfunction(data_count,data_tumor)
}
