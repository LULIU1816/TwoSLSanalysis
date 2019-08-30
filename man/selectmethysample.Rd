\name{selectmethysample}
\alias{selectmethysample}
\title{selectmethysample}
\usage{
selectmethysample(methylationLUAD)
}
\description{
select the tumor and normal tissues in DNA  methylation
}
\examples{
library("data.table")

#read methylation 
data data_methy<-fread("LUSC__methylation_450__tissueTypeAll__20181024055604.txt",data.table=F,head=T) data_methy <- as.data.frame(data_methy)[,-c(2,3,4)]

# select sample 01A,11A 
methylationLUAD <- selectmethysample(methylationLUAD)
}
