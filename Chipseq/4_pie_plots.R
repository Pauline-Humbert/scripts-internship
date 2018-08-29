#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)

### PIE PLOTS

setwd("/home/phumbert/raw_data/Chipseq/GenomeOntology_dir/") # to be in the directory where the file of reference is to select and plot the basic part of genome
dir_ontology <- list.files() # list all files in the directory
#read in R the next file and transform it in dataframe 
file_1 <- read.table("basic.genomeOntology.txt", col.names=c("Name","PeakFile/Annotation","#features[ref=8701]","Coverage(bp)[ref=3847448]","AvgFeatureSize[ref=442]","Overlap(#peaks)","Overlap(bp)","Expected_Overlap(bp,gsize=2.00e+09)","Log_Ratio_Enrichment","Log_P-value(+ underrepresented)","P-value"), fill=TRUE, skip=1)
file_1 <- data.frame(file_1)

# select just column Name of the previous dataframe
ref <- file_1$Name
# rename these two value to be the same Name used in the results "_stat.txt" files created previously in 3_peak_calling_annotation via annoPeak
levels(ref)[levels(ref) == "utr5"]<- "5UTR"
levels(ref)[levels(ref) == "utr3"]<-"3UTR"

# save the directory where are all results files of 3_peak_calling_annotation
path = "/home/phumbert/raw_data/Chipseq/MACS2results/"
dir.name <- dir(path)

pdf(file="/home/phumbert/raw_data/Chipseq/pie_plots.pdf")
# loop to list all directories (results for elba1, elba2, elba3 and insv) and for each select just "_stat.txt" file created previously in 3_peak_calling_annotation via annoPeak
for (dir in dir.name)
{
	in_dir <- paste("/home/phumbert/raw_data/Chipseq/MACS2results/",dir,"/",sep="")
	setwd(in_dir) # te be in the good direcory to select and read after "_stat.txt" file

	#list all files and just select "_stat.txt" file
	files <- list.files()
	file_2 <- grep("_stat.txt$", files, perl=TRUE, value=TRUE)

	#read in R the file and transform it in dataframe
	anno_data <- read.table(file_2, col.names=c("Annotation","Number_of_peaks","Total_size_(bp)","Log2_Enrichment"), fill=TRUE, skip=11)
	anno_data <- data.frame(anno_data)

	#rename these values in column "Annotation" to have the same name of rows (which will be retained) like in "basic.genomeOntology.txt" file => ref
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "Promoter"] <- "promoters"
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "Intron"] <- "introns"
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "Exon"] <- "exons"
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "TTS"] <- "tts"
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "CpG-Island"] <- "cpgIsland"
	levels(anno_data$Annotation)[levels(anno_data$Annotation) == "Intergenic"] <- "intergenic"

	#select in anno_data (current dataframe) just the rows which are the same name in Annotion and in ref and save the result in new dataframe
	new_anno_data <- subset(anno_data, Annotation %in% ref)

	#add two columns with the percent of peaks for each rows (in each part of the genome) and with a label which contains Annotation and percent
	new_anno_data <- new_anno_data%>%mutate(percent = paste(round((Number_of_peaks/ sum(Number_of_peaks))*100,1), "%", sep=""))
	new_anno_data <- new_anno_data%>%mutate(labels = paste(Annotation, " (",percent, ")", sep=""))

	#save the last dataframe in a new name to distinc the 4 dataframes created at the end of the loop (1 elba1, 1 elba2, 1 elba3 and 1 insv)
	assign(paste("new_anno_data_",dir, sep=""), new_anno_data)
}	

#code and command to plot the previous 4 dataframes with ggplot 
ggplot(new_anno_data_wtElba1_elba1Elba1, aes(x=factor(1), y=Number_of_peaks, fill=labels))+ 
	geom_bar(width=1, stat="identity")+ 
	scale_fill_manual(values = c("Navy", "Yellow", "Green", "Orange", "Purple","Red", "Pink","Cyan","Blue","Brown","Magenta"))+
	labs(x="", y="", title = "Pie Plot Proportion of Peaks Elba1 \n", fill = "Annotation")+
	coord_polar("y",start=0)

ggplot(new_anno_data_wtElba2_elba2Elba2, aes(x=factor(1), y=Number_of_peaks, fill=labels))+ 
	geom_bar(width=1, stat="identity")+ 
	scale_fill_manual(values = c("Navy", "Yellow", "Green", "Orange", "Purple","Red", "Pink","Cyan","Blue","Brown","Magenta"))+
	labs(x="", y="", title = "Pie Plot Proportion of Peaks Elba2 \n", fill = "Annotation")+
	coord_polar("y",start=0)

ggplot(new_anno_data_wtElba3_elba3Elba3, aes(x=factor(1), y=Number_of_peaks, fill=labels))+ 
	geom_bar(width=1, stat="identity")+ 
	scale_fill_manual(values = c("Navy", "Yellow", "Green", "Orange", "Purple","Red", "Pink","Cyan","Blue","Brown","Magenta"))+
	labs(x="", y="", title = "Pie Plot Proportion of Peaks Elba3 \n", fill = "Annotation")+
	coord_polar("y",start=0)


ggplot(new_anno_data_wtInsv_insvInsv, aes(x=factor(1), y=Number_of_peaks, fill=labels))+ 
	geom_bar(width=1, stat="identity")+ 
	scale_fill_manual(values = c("Navy", "Yellow", "Green", "Orange", "Purple","Red", "Pink","Cyan","Blue","Brown","Magenta"))+
	labs(x="", y="", title = "Pie Plot Proportion of Peaks Insv \n", fill = "Annotation")+
	coord_polar("y",start=0)
																																																
dev.off()