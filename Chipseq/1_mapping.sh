#!/bin/bash


cd /home/phumbert/raw_data/Chipseq/ #put directly in the correct directory

output=/home/phumbert/raw_data/Chipseq/FASTQCresults/ #variable saves the directory where files created will be put
for file in '/home/phumbert/raw_data/Chipseq/*.fastq.gz'
# loop to analyze the quality of the fastq files
do
fastqc -f fastq -o ${output} ${file}
#write the report to the new directory (variable output)
done 
 
for FName in `ls *[^.forward.].fastq.gz`
#loop to run different commands gene by gene
do
f=`basename $FName .fastq.gz` #retrieve only the name of gene
echo $f 

trimmomatic SE -threads 3 -phred33 $f.fastq.gz $f.forward.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:5:1:true LEADING:3 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:30

bowtie2 -p 5 -x /home/phumbert/genome_index/bowtie2_index_dm6/genome $f.forward.fastq.gz -S $f.sam

samtools view -bSh -@ 3 $f.sam -o $f.bam #compress more datas via .bam
samtools sort -@ 3 $f.bam -o $f.sorted.bam #sort the files
samtools index $f.sorted.bam #index(organize) the files sorted

done
