ChAMPfunction <- function(testDir,core=1){
  # CHAMP
myLoad <- champ.load(testDir, arraytype = "450K")
# distrubution of CpGs CpG.GUI() QC QC.GUI(arraytype = '450K') champ.QC() Normalization
myNorm <- champ.norm(beta = myLoad$beta, arraytype = "450K", cores = core)
# save write.csv(myNorm,file='./LUADNormalization Data.csv',quote=F,row.names = T) 
# SVD champ.SVD(beta=myNorm,pd=myLoad$pd)
# different methylaton position
myDMP <- champ.DMP(beta = myNorm, pheno = myLoad$pd$Sample_Group)
# DMP.GUI()
myDMPLUAD <- myDMP[[1]]
return(myDMPLUAD)
}
