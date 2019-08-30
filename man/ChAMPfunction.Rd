\name{ChAMPfunction}
\alias{ChAMPfunction}
\title{ChAMPfunction}
\usage{
ChAMPfunction(testDir,core=1)
}
\description{
perform a differentially methylated genes analysis required ChAMP package
}
\examples{
library("data.table")
library("ChAMP")

#configure path of GDC data
testDir = "~/GDCdata"

myDMPLUAD <- ChAMPfunction(testDir)
}
