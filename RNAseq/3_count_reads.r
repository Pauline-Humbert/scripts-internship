#!/usr/bin/env Rscript

# Count #reads per gene for all genes => raw data
setwd("/home/phumbert/raw_data/RNAseq/")   

library(Rsubread) # load the library

files <- list.files() # list all files in the directory 
sorted_bam_files <- grep(".sorted.bam$", files, perl=TRUE, value=TRUE) 
# save in the variable "sorted_bam_files" all files in "files" for each gene which are .sorted.bam
print (sorted_bam_files)

fc_PE <- featureCounts(sorted_bam_files, annot.ext = "/home/phumbert/genome_index/dm6flybase/dmel-all-r6.12_chr.gtf", 
GTF.attrType="gene_symbol", nthreads=5, minMQS=10, countMultiMappingReads = F, strandSpecific = 1, isGTFAnnotationFile=T,
GTF.featureType = "exon", useMetaFeatures=TRUE, primaryOnly=TRUE, isPairedEnd=TRUE)

# annot.ext => needs gene annotation file (.gtf) (with the path if it isn't in the same directory)
# GTF.attrType="gene_id" => the label in the gtf describing a gene
# nthreads=5 => use multiple cores for speed
# minMQS=10 => exclude low quality
# countMultiMappingReads = F => If specified, multi-mapping reads/fragments will be counted.
# strandSpecific = 1 => this is strand-specific data so take strand into account
# isGTFAnnotationFile=T => Specify the format of the annotation file.  Acceptable formats include ‘GTF’ and ‘SAF’.
# GTF.featureType = "exon" => Specify the feature type. Only rows which have the matched feature type in the provided GTF annotation file will be included for read counting.
# useMetaFeatures=TRUE => This indicates we want expression over genes; for exon-level expression set this to FALSE 
# primaryOnly=TRUE => only include primary mappings
# isPairedEnd=TRUE => this is pair-ended data

save(fc_PE, file = "Count_reads_raw_data.RData")
