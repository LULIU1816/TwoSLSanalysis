edgeRfunction <- function(case_control_data,fdr=0.001,logfc=2.5,logcpm=2.5){
  group <- factor(c(rep("H", (ncol(case_control_data) - 1)/2), rep("M", (ncol(case_control_data) - 1)/2)))
# remove genes with low expression
case_control_data1 <- case_control_data[rowSums(cpm(case_control_data[, 2:ncol(case_control_data)]) > 1) >= 2, ]
c_c_data <- DGEList(counts = case_control_data1[, 2:ncol(case_control_data)], genes = case_control_data1[, 1], group = group)
## flitration and normalization
c_c_data <- calcNormFactors(c_c_data)
## design matrix
design <- model.matrix(~group)
## dispersion
exprSet <- estimateCommonDisp(c_c_data, design)
exprSet <- estimateTagwiseDisp(exprSet)
## different expression gene(DEG)
et <- exactTest(exprSet)
tTag <- topTags(et, n = nrow(exprSet))
diff_gene_edgeR <- subset(tTag$table, FDR < fdr & (logFC > logfc | logFC < -logfc) & (logCPM > logcpm | logCPM < -logcpm))
diff_gene <- diff_gene_edgeR[, 1]
w <- gregexpr("\\|", diff_gene)
p <- NULL
for (i in 1:length(diff_gene)) {
  p[i] <- substr(diff_gene[i], 1, w[[i]][1] - 1)
}
diff_gene_edgeR[, 1] <- p
return(diff_gene_edgeR)
}
