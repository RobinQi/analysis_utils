library(edgeR)

data = read.table("/Users/robin/Desktop/software/deg/./raw_counts.txt", header=T, row.names=1, com='')
col_ordering = c(2,3,4,6,8,9,11,12,13,14,15,16,17,18,19,20,21,22,24,25,30,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58,60,62,65,66,67,68,69,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,87,88,89,90,91,92,93,94,95,96,99,100,101,102,103,104,106,107,108,109,110,111,113,115,116,117,118,119,120,122,123,124,125,127,128,129,130,132,133,134,135,136,137,138,139,140,142,143,144,145,147,148,149,150,151,152,153,154,155,156,157,158,159,160,1,5,7,10,23,26,27,28,29,31,46,59,61,63,64,70,86,97,98,105,112,114,121,126,131,141,146,161)
rnaseqMatrix = data[,col_ordering]
rnaseqMatrix = round(rnaseqMatrix)
rnaseqMatrix = rnaseqMatrix[rowSums(cpm(rnaseqMatrix) > 1) >= 2,]
conditions = factor(c(rep("G1", 133), rep("G2", 28)))

exp_study = DGEList(counts=rnaseqMatrix, group=conditions)
exp_study = calcNormFactors(exp_study)
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)
et = exactTest(exp_study, pair=c("G1", "G2"))
tTags = topTags(et,n=NULL)
result_table = tTags$table
result_table = data.frame(sampleA="G1", sampleB="G2", result_table)
result_table$logFC = -1 * result_table$logFC
write.table(result_table, file='raw_counts.txt.G1_vs_G2.edgeR.DE_results', sep='	', quote=F, row.names=T)
write.table(rnaseqMatrix, file='raw_counts.txt.G1_vs_G2.edgeR.count_matrix', sep='	', quote=F, row.names=T)
source("/Users/robin/Desktop/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/R/rnaseq_plot_funcs.R")
pdf("raw_counts.txt.G1_vs_G2.edgeR.DE_results.MA_n_Volcano.pdf")
plot_MA_and_Volcano(result_table$logCPM, result_table$logFC, result_table$FDR)
dev.off()
