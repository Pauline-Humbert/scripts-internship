#!/usr/bin/env Rscript


library(Biostrings)

setwd("/home/phumbert/raw_data/Chipseq/Motif_Analysis/")

dir.name <- list.dirs(path = "/home/phumbert/raw_data/Chipseq/Motif_Analysis", full.names = FALSE, recursive = FALSE)

for (dir in dir.name)
{
	in_dir <- paste("/home/phumbert/raw_data/Chipseq/Motif_Analysis/",dir,"/1_PWM_motif/",sep="")
    setwd(in_dir)

	fasta_file <- paste(dir,"_peaks_narrowPeak.fa",sep="")
	#  rewrite the fasta file from 7_fasta.sh with the name of peak id instead of the localization on chromosome
    dseq = readDNAStringSet(fasta_file)
	#dseq
	dseq = DNAStringSet(toupper(dseq))
	names(dseq) = gsub("::(.*)", "", names(dseq), perl=T)
	writeXStringSet(dseq, file=fasta_file, format="fasta", width=20000) 
}

#cf 9_motif_peaks.sh to run fimo