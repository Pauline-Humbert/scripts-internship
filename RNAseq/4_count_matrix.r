#!/usr/bin/env Rscript
setwd("/home/phumbert/raw_data/RNAseq/")

##### Count # reads per gene for all genes => matrix

library(Rsubread)
library(limma)
library(edgeR)
library(openxlsx)

load("Count_reads_raw_data.RData") # fc_PE below from this file

# CREATE A DGEList (data.frame) => count matrix and information by genes  
x <- DGEList(counts=fc_PE$counts, genes=fc_PE$annotation)
dim(x)

# ORGANISING DATA
samplenames <- substr(colnames(x), 33, (nchar(colnames(x))-11)) # select just the name of gene in my directory
colnames(x) <- samplenames # replace the name of columns in x by just the name of my genes

# REMOVING GENES THAT ARE LOWLY EXPRESSED
# table(rowSums(x$counts==0)==ncol(x)) #de-comment if you want the result
#keep.exprs <- rowSums(x$counts)>=50 
#x <- x[keep.exprs,, keep.lib.sizes=FALSE]
x <- x[rowSums(x$counts>=10)>3,] # better way to filter genes lowly expressed than the commands in comments
dim(x)

# NORMALISING GENE EXPESSION DISTRIBUTIONS
x <- calcNormFactors(x)

matrix_after_normalising <- x$counts # to save it at the end

# DIFFERENTIEAL EXPRESSION ANALYSIS
### Creating a design matrix
######Creating samples
y <- colnames (x)
samples <- NULL
for (var in y) {
  if (grep(".*_.", var)) samples <- c(samples, (substr(var, 1, ((nchar(var))-2))))
  else print (var)
  }
  
######Design matrix
samples = factor(samples, levels = c("elba1", "elba2", "elba3", "insv", "yw"))  
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)

setwd("/home/phumbert/raw_data/RNAseq/results_DEG/") # to save results and plots in the right directory 

### Removing heteroscedascity from count data
######Transform RNA-Seq Data Ready for Linear Modelling
pdf("Rplots_DEA.pdf")
v <- voom(x,design,plot=TRUE)

######Linear model is fitted to the expression genes
vfit <- lmFit(v,design)

contr.matrix <- makeContrasts (
	elba1VSyw = elba1-yw, 
	elba2VSyw = elba2-yw, 
	elba3VSyw = elba3-yw, 
	insvVSyw = insv-yw, 
	levels = colnames(design))

vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

######To obtain more precise estimates of gene-wise variability
efit <- eBayes(vfit)

plotSA(efit, main="Mean_variance trend")


### Examining the number of DE genes
######For a quick look at DE levels
summary(decideTests(efit)) #related t-statistics as up, down or not significant 


######Extracting and writing results
write.fit(efit, decideTests(efit), file="results_numberDE_genes.txt") #not very important to save that but for a quick look it saves in .txt file 

### Examining individual DE genes from top to bottom
versus= c("elba1VSyw", "elba2VSyw", "elba3VSyw", "insvVSyw")
listDEGenes=list()
for (vs in versus) {
	listDEGenes[[vs]]= topTreat(efit, coef=vs, number=nrow(x$counts))
} #Extract a table of the top-ranked genes from a linear model fit

###### False discovery rate
listDEGenes_fdr05 <- lapply(listDEGenes, function(x) x[x[,"adj.P.Val"]<=0.05,]) # limit to FDR of 0.05
lapply( listDEGenes_fdr05, nrow)

save(matrix_after_normalising, design, v, contr.matrix, vfit, efit, listDEGenes, listDEGenes_fdr05, file = "matrix_DEG.RData")
write.xlsx(listDEGenes, file = "expression_genelevel.xlsx",  sep="",row.names=T, StartRow=1, StartCol=1) 

### Graphical representation of DE results

######Summarising results for all genes visually:
plotMD(efit, column=1, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main=colnames(efit)[1])

plotMD(efit, column=2, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main=colnames(efit)[2])

plotMD(efit, column=3, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main=colnames(efit)[3])

plotMD(efit, column=4, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main=colnames(efit)[4])                                                         

###### Plot counts against p-value for all genes visually: 
hist(listDEGenes[[1]]$P.Value, breaks=50, main = "elba1VSyw", xlab="p.value", ylab="counts", col="blue", border="green")
hist(listDEGenes[[2]]$P.Value, breaks=50, main = "elba2VSyw", xlab="p.value", ylab="counts", col="blue", border="green")
hist(listDEGenes[[3]]$P.Value, breaks=50, main = "elba3VSyw", xlab="p.value", ylab="counts", col="blue", border="green")
hist(listDEGenes[[4]]$P.Value, breaks=50, main = "insvVSyw", xlab="p.value", ylab="counts", col="blue", border="green")

dev.off()
