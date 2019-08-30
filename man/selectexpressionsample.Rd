\name{selectexpressionsample}
\alias{selectexpressionsample}
\title{selectexpressionsample}
\usage{
selectexpressionsample(gene_id,expressionLUAD)
}
\description{
select the tumor and normal tissues in expression
}
\examples{
library("data.table")

#read expression data
data_normalized<-fread("LUAD__gene.normalized_RNAseq__tissueTypeAll__20181022065043.txt",data.table=F,head=T)
data_normalized <- as.data.frame(data_normalized)
gene_id <- data_normalized[, 1]

# select sample 01A,11A 
selectexpressionsample(gene_id,data_normalized)
}
