\name{edgeRfunction}
\alias{edgeRfunction}
\title{edgeRfunction}
\usage{
edgeRfunction(case_control_data,fdr=0.001,logfc=2.5,logcpm=2.5)
}
\description{
 perform a differentially expressed genes analysis required edgeR package
}
\examples{
library("data.table")
library("edgeR")

#read expression data 
data_normalized<-fread("LUSC__gene.normalized_RNAseq__tissueTypeAll__20181022104424.txt",data.table=F,head=T) 
data_count<-fread("LUSC__gene_RNAseq__tissueTypeAll__20181022105257.txt",data.table=F,head=T) 
data_normalized <- as.data.frame(data_normalized) 
data_count <-as.data.frame(data_count)

# structure case-control data
case_control_data <- casecontrolfunction(data_count,data_tumor)

# edgeR
diff_gene_edgeR <- edgeRfunction(case_control_data)
diff_gene <- diff_gene_edgeR[,1]
}
