#!/usr/bin/env Rscript
setwd("/home/phumbert/raw_data/RNAseq/results_DEG")

library(limma)
library(stringr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

load("matrix_DEG.RData")


# PATHWAY ANALYSIS

setwd("/home/phumbert/genesets/")
path = "/home/phumbert/genesets/" # to save this path for retrieve files in this directory after

## using  ".NAME.RData" files only to create index
file.name <- dir(path, pattern = ".NAME.RData") # to list all ".Name.RData" file names in this directory
 
for (file in file.name) # loop to use different RData files as index in different commands (camera and mroast in the same loop)  
{
	load(file) # load index of each file
	indx <- ids2indices(geneset_name,id=rownames(v)) # save index in a variable

### Camera command: differentially expressed genes tend to be over-represented in the gene set, compared to all the other genes in the experiment 

	cam.elba1VSyw <- camera(v,index=indx,design=design,contrast=contr.matrix[,1], inter.gene.cor=0.05) 
	cam.elba2VSyw <- camera(v,index=indx,design=design,contrast=contr.matrix[,2], inter.gene.cor=0.05)
	cam.elba3VSyw <- camera(v,index=indx,design=design,contrast=contr.matrix[,3], inter.gene.cor=0.05)	
	cam.insvVSyw <- camera(v,index=indx,design=design,contrast=contr.matrix[,4], inter.gene.cor=0.05)

	namefile <- str_sub(file,1, str_length(file)-11) # cut the end of the file'sname (.NAME.RData) and namefile is used in the new file's name
	save(cam.elba1VSyw,cam.elba2VSyw,cam.elba3VSyw,cam.insvVSyw,file=paste("/home/phumbert/raw_data/RNAseq/pathway_analysis/cam_",namefile,".RData", sep="")) 

### mroast command: Are the genes in the set differentially expressed as a whole? (similar to camera but slower and more significant)

	rst.elba1VSyw <- mroast(v, index=indx, design=design, contrast=contr.matrix[,1], geneid=rownames(v), set.statistic = "mean", nrot = 999, approx.zscore = TRUE, adjust.method = "BH", midp = TRUE, sort = "directional")
	rst.elba2VSyw <- mroast(v, index=indx, design=design, contrast=contr.matrix[,2], geneid=rownames(v), set.statistic = "mean", nrot = 999, approx.zscore = TRUE, adjust.method = "BH", midp = TRUE, sort = "directional")
	rst.elba3VSyw <- mroast(v, index=indx, design=design, contrast=contr.matrix[,3], geneid=rownames(v), set.statistic = "mean", nrot = 999, approx.zscore = TRUE, adjust.method = "BH", midp = TRUE, sort = "directional")
	rst.insvVSyw <- mroast(v, index=indx, design=design, contrast=contr.matrix[,4], geneid=rownames(v), set.statistic = "mean", nrot = 999, approx.zscore = TRUE, adjust.method = "BH", midp = TRUE, sort = "directional")

	namefile <- str_sub(file,1, str_length(file)-11) # cut the end of the file'sname (.NAME.RData) and namefile is used in the new file's name
	save(rst.elba1VSyw,rst.elba2VSyw,rst.elba3VSyw,rst.insvVSyw,file=paste("/home/phumbert/raw_data/RNAseq/pathway_analysis/rst_",namefile,".RData", sep="")) 
}

## Select rows which include the word "pathway" just for results from msigdb geneset
setwd("/home/phumbert/raw_data/RNAseq/pathway_analysis/")

# For msigdb geneset from "camera"
load("cam_fly.msigdb.RData")

cam.elba1VSyw <-subset(cam.elba1VSyw, row.names(cam.elba1VSyw) %in% grep("pathway", row.names(cam.elba1VSyw), value=TRUE, ignore.case=TRUE))
cam.elba2VSyw <-subset(cam.elba2VSyw, row.names(cam.elba2VSyw) %in% grep("pathway", row.names(cam.elba2VSyw), value=TRUE, ignore.case=TRUE))
cam.elba3VSyw <-subset(cam.elba3VSyw, row.names(cam.elba3VSyw) %in% grep("pathway", row.names(cam.elba3VSyw), value=TRUE, ignore.case=TRUE))
cam.insvVSyw <-subset(cam.insvVSyw, row.names(cam.insvVSyw) %in% grep("pathway", row.names(cam.insvVSyw), value=TRUE, ignore.case=TRUE))

save(cam.elba1VSyw,cam.elba2VSyw,cam.elba3VSyw,cam.insvVSyw,file="/home/phumbert/raw_data/RNAseq/pathway_analysis/cam_fly.msigdb.RData")

# For msigdb geneset from "mroast"
load("rst_fly.msigdb.RData")

rst.elba1VSyw <-subset(rst.elba1VSyw, row.names(rst.elba1VSyw) %in% grep("pathway", row.names(rst.elba1VSyw), value=TRUE, ignore.case=TRUE))
rst.elba2VSyw <-subset(rst.elba2VSyw, row.names(rst.elba2VSyw) %in% grep("pathway", row.names(rst.elba2VSyw), value=TRUE, ignore.case=TRUE))
rst.elba3VSyw <-subset(rst.elba3VSyw, row.names(rst.elba3VSyw) %in% grep("pathway", row.names(rst.elba3VSyw), value=TRUE, ignore.case=TRUE))
rst.insvVSyw <-subset(rst.insvVSyw, row.names(rst.insvVSyw) %in% grep("pathway", row.names(rst.insvVSyw), value=TRUE, ignore.case=TRUE))

save(rst.elba1VSyw,rst.elba2VSyw,rst.elba3VSyw,rst.insvVSyw,file="/home/phumbert/raw_data/RNAseq/pathway_analysis/rst_fly.msigdb.RData")

## HEATMAP 

### From camera results:
pdf(file="Heatmaps(1).pdf", width=29, height=24)
path = "/home/phumbert/raw_data/RNAseq/pathway_analysis/" # to save this path for retrieve files in this directory after


all.files <- dir(path, pattern = "cam_") # to list all "cam_" file names in this directory

for (file in all.files) # loop to load all cam_"".RData files one by one to create for each one a table for heatmap
{
	load(file)
	# create a new data frame with only PValue and direction of elba1,2,3 and insv
	elba1 <- subset(cam.elba1VSyw, select=c("Direction","PValue"))
	elba2 <- subset(cam.elba2VSyw, select=c("Direction","PValue"))
	elba3 <- subset(cam.elba3VSyw, select=c("Direction","PValue"))
	insv <- subset(cam.insvVSyw, select=c("Direction","PValue"))

	all <- data.frame(elba1,elba2,elba3,insv)

	#select just data with PValue<=0.05 for elba1, or elba2 or elba3 or insv
	select_PValue <- subset(all, elba1$PValue<=0.05 | elba2$PValue<=0.05 | elba3$PValue<=0.05 | insv$PValue<=0.05)
	colnames(select_PValue)<-c("elba1_Direction","elba1_PValue","elba2_Direction","elba2_PValue","elba3_Direction","elba3_PValue","insv_Direction","insv_PValue")

	# rank rows according to PValues 
	select_rankPValue <- select_PValue[order(select_PValue$elba1_PValue, select_PValue$elba2_PValue, select_PValue$elba3_PValue, select_PValue$insv_PValue), , drop = FALSE] 
			# drop = FALSE for keeping row.names
	
	#select top 50 genes for heatmaps
	ifelse(nrow(select_rankPValue)<=50, top_50_genes <- select_rankPValue, top_50_genes <- select_rankPValue[1:50,])

	# Apply -log10() or log10() on each PValue in elba 1, 2, 3 and insv according to genes are up or down-regulated in elba 1, 2, 3 and insv
	# => distinct genes up and down regulated in heatmap
	top_50_genes$elba1_PValue <- ifelse(top_50_genes$elba1_Direction == "Up", -log10(top_50_genes$elba1_PValue), log10(top_50_genes$elba1_PValue))
	top_50_genes$elba2_PValue <- ifelse(top_50_genes$elba2_Direction == "Up", -log10(top_50_genes$elba2_PValue), log10(top_50_genes$elba2_PValue))
	top_50_genes$elba3_PValue <- ifelse(top_50_genes$elba3_Direction == "Up", -log10(top_50_genes$elba3_PValue), log10(top_50_genes$elba3_PValue))
	top_50_genes$insv_PValue <- ifelse(top_50_genes$insv_Direction == "Up", -log10(top_50_genes$insv_PValue), log10(top_50_genes$insv_PValue))

	# create table for heatmaps without columns Directions
	heatmap_table <- subset(top_50_genes, select=c("elba1_PValue","elba2_PValue","elba3_PValue","insv_PValue")) 
	colnames(heatmap_table)<- c("elba1","elba2","elba3","insv")

	heatmap_table <- as.matrix(heatmap_table)

	# heatmap plots
	name <- str_sub(file, start=5, end=-7)
	red.pal <- RColorBrewer::brewer.pal(9, "Reds")
	pheatmap(heatmap_table, color= red.pal, cluster_rows= TRUE, cluster_cols=TRUE, cellwidth=70, cellheight=30, fontsize=20, fontsize_row=25, fontsize_col=25, main=paste("Heatmap_PValue_",name, sep=""))	
}
dev.off()

###From mroast results:
pdf(file="Heatmaps(2).pdf", width=29, height=24)

path = "/home/phumbert/raw_data/RNAseq/pathway_analysis/" # to save this path for retrieve files in this directory after

all.files <- dir(path, pattern = "rst_") # to list all "rst_" file names in this directory

for (file in all.files) # loop to load all rst_"".RData files one by one to create for each one a table for heatmap
{
	load(file)
	# create a new data frame with only PValue and direction of elba1,2,3 and insv
	elba1 <- subset(rst.elba1VSyw, select=c("Direction","PValue"))
	elba2 <- subset(rst.elba2VSyw, select=c("Direction","PValue"))
	elba3 <- subset(rst.elba3VSyw, select=c("Direction","PValue"))
	insv <- subset(rst.insvVSyw, select=c("Direction","PValue"))

	all <- data.frame(elba1,elba2,elba3,insv)

	#select just data with PValue<=0.05 for elba1, or elba2 or elba3 or insv
	select_PValue <- subset(all, elba1$PValue<=0.05 | elba2$PValue<=0.05 | elba3$PValue<=0.05 | insv$PValue<=0.05)
	colnames(select_PValue)<-c("elba1_Direction","elba1_PValue","elba2_Direction","elba2_PValue","elba3_Direction","elba3_PValue","insv_Direction","insv_PValue")

	# rank rows according to PValues 
	select_rankPValue <- select_PValue[order(select_PValue$elba1_PValue, select_PValue$elba2_PValue, select_PValue$elba3_PValue, select_PValue$insv_PValue), , drop = FALSE] 
			# drop = FALSE for keeping row.names
	
	#select top 50 genes for heatmaps
	ifelse(nrow(select_rankPValue)<=50, top_50_genes <- select_rankPValue, top_50_genes <- select_rankPValue[1:50,])

	# Apply -log10() or log10() on each PValue in elba 1, 2, 3 and insv according to genes are up or down-regulated in elba 1, 2, 3 and insv
	# => distinct genes up and down regulated in heatmap
	top_50_genes$elba1_PValue <- ifelse(top_50_genes$elba1_Direction == "Up", -log10(top_50_genes$elba1_PValue), log10(top_50_genes$elba1_PValue))
	top_50_genes$elba2_PValue <- ifelse(top_50_genes$elba2_Direction == "Up", -log10(top_50_genes$elba2_PValue), log10(top_50_genes$elba2_PValue))
	top_50_genes$elba3_PValue <- ifelse(top_50_genes$elba3_Direction == "Up", -log10(top_50_genes$elba3_PValue), log10(top_50_genes$elba3_PValue))
	top_50_genes$insv_PValue <- ifelse(top_50_genes$insv_Direction == "Up", -log10(top_50_genes$insv_PValue), log10(top_50_genes$insv_PValue))

	# create table for heatmaps without columns Directions
	heatmap_table <- subset(top_50_genes, select=c("elba1_PValue","elba2_PValue","elba3_PValue","insv_PValue")) 
	colnames(heatmap_table)<- c("elba1","elba2","elba3","insv")

	heatmap_table <- as.matrix(heatmap_table)

	# heatmap plots
	name <- str_sub(file, start=5, end=-7)
	red.pal <- RColorBrewer::brewer.pal(9, "Reds")
	pheatmap(heatmap_table, color= red.pal, cluster_rows= TRUE, cluster_cols=TRUE, cellwidth=70, cellheight=30, fontsize=20, fontsize_row=25, fontsize_col=25, main=paste("Heatmap_PValue_",name, sep=""))	
}

dev.off()