\name{overlapaverage}
\alias{overlapaverage}
\title{overlapaverage}
\usage{
overlapaverage(CpG,methylationLUAD)
}
\description{
average the repeated measures data of LUAD methylation if exisit, and concrete analysis should be made according to concrete circumstance in other cancer
}
\examples{
library("data.table")

#read methylation data 
data_methy<-fread("LUSC__methylation_450__tissueTypeAll__20181024055604.txt",data.table=F,head=T) data_methy <- as.data.frame(data_methy)[,-c(2,3,4)]

 methylationLUAD <- overlapaverage(CpG,methylationLUAD)
}
