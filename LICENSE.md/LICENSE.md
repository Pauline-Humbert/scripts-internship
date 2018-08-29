#!/bin/bash

cd /home/phumbert/raw_data/RNAseq/ #put directly in the correct directory 

output=/home/phumbert/raw_data/RNAseq/FASTQCresults/ #variable saves the directory where files created will be put
for file in '/home/phumbert/raw_data/RNAseq/*.R?.fastq.gz'
# loop to analyze the quality of the fastq files
do
fastqc -f fastq -o ${output} ${file} 
# write the report to the new directory (variable output)  
done

for FName in `ls *.R1.fastq.gz` 
#loop to run different commands gene by gene  
do
f=`basename $FName .R1.fastq.gz` #retrieve only the name of gene 
echo $f

trimmomatic PE -threads 3 -phred33 $f.R1.fastq.gz $f.R2.fastq.gz $f.pair1.fastq.gz $f.unpair1.fastq.gz $f.pair2.fastq.gz $f.unpair2.fastq.gz ILLUMINACLIP:TrueSeq_adaptor.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30
# trim adaptor on .R1.fastq.gz and .R2.fastq.gz files depending of gene, and create 4 files: 2 pair.fastq.gz and 2 unpair.fastq.gz

hisat2 --fr --threads 10  --rna-strandness RF -x /home/phumbert/genome_index/dm6flybase/dm6flybase -1 $f.pair1.fastq.gz -2  $f.pair2.fastq.gz -S $f.sam
# it's to align the diffeent sequences with the whole genome

samtools view -bSh -@ 3 $f.sam -o $f.bam #compress more datas via .bam
samtools sort -@ 3 $f.bam -o $f.sorted.bam #sort the files 
samtools index $f.sorted.bam #index(organize) the files sorted

done
