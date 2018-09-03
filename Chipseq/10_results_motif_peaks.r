#!/usr/bin/env Rscript                                                                                                                                      2

library(dplyr)
library(ggplot2)
library(Rmisc)

setwd("/home/phumbert/raw_data/RNAseq/results_DEG/") 
load("matrix_DEG.RData") # load data (expression genes level) in this file in the previous directory


setwd("/home/phumbert/raw_data/Chipseq/Motif_Analysis/") # change directory

dir.name <- list.dirs(path = "/home/phumbert/raw_data/Chipseq/Motif_Analysis", full.names = FALSE, recursive = FALSE)

pdf(file="barplots.pdf")
#loop to do all these commands for Elba1,2,3 and insv
for (dir in dir.name)
{
	in_dir <- paste("/home/phumbert/raw_data/Chipseq/Motif_Analysis/",dir,"/1_PWM_motif/",sep="")
	setwd(in_dir)

	if (dir == "wtElba1_elba1Elba1") {
		n <- 1 # number for the dataframe in the list (expression genes level) from "matrix_DEG.RData"
		TF="elba1"
	} else if (dir == "wtElba2_elba2Elba2"){
		n <- 2
		TF="elba2"
	} else if (dir == "wtElba3_elba3Elba3"){
		n <- 3
		TF="elba3"
	} else if (dir == "wtInsv_insvInsv"){
		n <- 4			
		TF="insv"
	}		

	listDEGenes[[n]] <- subset(listDEGenes[[n]], select=c("GeneID", "logFC")) #keep only in this list the two columns GeneID and logFC

	#Read peakannofile to access on genes
	genes <- paste(dir,"_peakannofile.txt",sep="")
	genes <- read.delim(genes, col.names=c("PeakID","Chr","Start","End","Strand","Peak_score","Focus_Ratio/Region_size","Annotation","Detailed_Annotation",
	"Distance_to_TSS", "Nearest_PromoterId", "EntreID", "Nearest_Unigene","Nearest_Refseq", "Nearest_Enembl", "Gene_Name", "Gene_Alias", 
	"Gene_Description","Gene_Type"),fill=TRUE)
	genes[7] <- NULL # supress columns without values
	genes[14] <- NULL

	#Read peaks file
	peaks <- paste(dir,"_peaks.narrowPeak",sep="")
	peaks <- read.table(peaks, col.names=c("chrom","chromStart","chromEnd","name","score","strand","signalValue","pValue","qValue","peak"), fill=TRUE)

	#Read motif file
	in_dir_fimo <- paste("/home/phumbert/raw_data/Chipseq/Motif_Analysis/",dir,"/1_PWM_motif/fimo/",sep="")
	setwd(in_dir_fimo)
	motif <- read.table("fimo.txt", col.names=c("motif_id","sequence_name","start","stop","strand","score","p-value","q-value","matched_sequence"), fill=TRUE, skip=1)
	

	#Genes containing all peaks without motif
	genes_peaks_NOmotif <- left_join(genes, peaks, by=c("PeakID"="name"))
	genes_peaks_NOmotif <- left_join(genes_peaks_NOmotif, motif, by=c("PeakID"="sequence_name")) # merge these different dataframes acccording to  arguments in "by"
	genes_peaks_NOmotif <- subset(genes_peaks_NOmotif, !(genes_peaks_NOmotif$PeakID %in% motif$sequence_name), select=colnames(genes_peaks_NOmotif))
	# keep genes with peaks which aren't in motif


	#Genes containing all peaks with motif
	genes_peaks_motif <- left_join(genes, peaks, by=c("PeakID"="name"))
	genes_peaks_motif <- left_join(genes_peaks_motif, motif, by=c("PeakID"="sequence_name")) # merge these different dataframes acccording to  arguments in "by"
	genes_peaks_motif <- subset(genes_peaks_motif, genes_peaks_motif$PeakID %in% motif$sequence_name, select=colnames(genes_peaks_motif))
	#keep genes with peacks which are in motif

	#Genes containing top 200 peaks with motif
	genes_top200peaks_motif <- genes_peaks_motif[order(genes_peaks_motif$Peak_score, decreasing=TRUE),] # order genes_peaks_motif with Peak_score 
	genes_top200peaks_motif <- genes_top200peaks_motif[1:200,] #select the 200 first rows

	#remove genes_top200peaks_motif from genes_peaks_motif
	genes_peaks_motif <- genes_peaks_motif[order(genes_peaks_motif$Peak_score, decreasing=TRUE),]
	genes_peaks_motif <- genes_peaks_motif[-(1:200),]

	#Genes containing any peaks
	genes_NOpeaks <- subset(listDEGenes[[n]], !(listDEGenes[[n]]$GeneID %in% genes_peaks_NOmotif$Gene_Name), select=colnames(listDEGenes[[n]])) # keep all genes in the list (expression genes level) which aren't in genes
	genes_NOpeaks <- subset(genes_NOpeaks, select=c("GeneID")) #keep only GeneID column
	colnames(genes_NOpeaks)[1] <- "Gene_Name" #rename the column

	#Homogenize tables and add logFC(expression)

	## just select Gene_Name column 
	genes_peaks_NOmotif <- subset(genes_peaks_NOmotif, select=c("Gene_Name"))
	genes_peaks_motif <- subset(genes_peaks_motif, select=c("Gene_Name"))
	genes_top200peaks_motif <- subset(genes_top200peaks_motif, select=c("Gene_Name"))

	##add column "set" (group)
	genes_NOpeaks <- genes_NOpeaks%>%mutate(set = "Genes_NOpeaks")
	genes_peaks_NOmotif <- genes_peaks_NOmotif%>%mutate(set = "Genes_allpeaks_Nomotif")
	genes_peaks_motif <- genes_peaks_motif%>%mutate(set = "Genes_allpeaks_allmotif")
	genes_top200peaks_motif <- genes_top200peaks_motif%>%mutate(set = "Genes_top200peaks_motif")

	#put 4 dataframes in 1 and merge it with list (expression gene level) to have access to logFC for each gene
	All <- rbind(genes_NOpeaks,genes_peaks_NOmotif,genes_peaks_motif,genes_top200peaks_motif)
	All <- left_join(All, listDEGenes[[n]], by=c("Gene_Name"="GeneID"))	
	All[,3] <- as.numeric(All[,3])
	
	# Summarize data: summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
	All_sum <- summarySE(All, measurevar="logFC", groupvar=c("set"), na.rm=TRUE)
	
	#save the last dataframe in a All_sum following elba1 elba2 elba3 and insv
	assign(paste("All_sum_",TF, sep=""), All_sum)
}

#Bar plots
ggplot(All_sum_elba1, aes(x=set, y=logFC))+ geom_bar(position=position_dodge(), stat="identity",size=.3, width=.4)+
	geom_errorbar(aes(ymin=logFC-se, ymax=logFC+se),size=.3, width=.1, position=position_dodge(.9))+
	xlab("Set")+ ylab("logFC")+ggtitle(paste("Genes Expression_peaks_motifs_elba1",sep=""))

ggplot(All_sum_elba2, aes(x=set, y=logFC))+ geom_bar(position=position_dodge(), stat="identity", colour="black",size=.3, width=.4)+
	geom_errorbar(aes(ymin=logFC-se, ymax=logFC+se),size=.3, width=.1, position=position_dodge(.9))+
	xlab("Set")+ ylab("logFC")+ggtitle(paste("Genes Expression_peaks_motifs_elba2",sep=""))

ggplot(All_sum_elba3, aes(x=set, y=logFC))+ geom_bar(position=position_dodge(), stat="identity", colour="black",size=.3, width=.4)+
	geom_errorbar(aes(ymin=logFC-se, ymax=logFC+se),size=.3, width=.1, position=position_dodge(.9))+
	xlab("Set")+ ylab("logFC")+ggtitle(paste("Genes Expression_peaks_motifs_elba3",sep=""))

ggplot(All_sum_insv, aes(x=set, y=logFC))+ geom_bar(position=position_dodge(), stat="identity", colour="black",size=.3, width=.4)+
	geom_errorbar(aes(ymin=logFC-se, ymax=logFC+se),size=.3, width=.1, position=position_dodge(.9))+
	xlab("Set")+ ylab("logFC")+ggtitle(paste("Genes Expression_peaks_motifs_insv",sep=""))

dev.off()
