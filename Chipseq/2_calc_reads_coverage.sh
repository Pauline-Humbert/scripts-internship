#!/bin/bash

cd /home/phumbert/raw_data/Chipseq/ #put directly in the correct directory

for FName in `ls *[^.forward.].fastq.gz` # loop to run different commands gene by gene
do
f=`basename $FName .fastq.gz` # retrieve only the gene's name
echo $f

# Calculate number of reads in .fastq files before trimming adaptor
r1=`gzip -dc $f.fastq.gz | wc -l` # gzip -dc to decompress .gz files compressed and wc -l gives the number of lines in the files
r1=`echo "scale=2; $r1 / 4" | bc -l`
#echo $r1
#r1=number of reads in this .fastq file / scale: number after "," and bc -l used for float number
echo "Number reads before trimming adaptor $r1" > /home/phumbert/raw_data/Chipseq/libsize/$f.libsize.txt
#save this value in a new file name_gene.libzise.txt (creation file ">")

# Calculate number of reads in .fastq files after trimming adaptor
r2=`gzip -dc $f.forward.fastq.gz | wc -l` # gzip -dc to decompress .gz files compressed and wc -l gives the number of lines in the files
r2=`echo "scale=2; $r2 / 4" | bc -l`
#echo $r2
#r2=number of reads in this .fastq file / scale: number after "," and bc -l used for float number
echo "Number reads after trimming adaptor $r2" >> /home/phumbert/raw_data/Chipseq/libsize/$f.libsize.txt
#save this value in he previous file (add information in a file ">>")

# Calculate number of mapped reads in .sorted.bam files
r3=`samtools view -@ 8 -F 0x105 -c $f.sorted.bam`
#echo $r3
echo "Number mapped reads: $r3" >> /home/phumbert/raw_data/Chipseq/libsize/$f.libsize.txt
#save this value in the previous file (add information in a file ">>")

# Statistical calculation
r=`echo "scale=1; ($r3 / $r2) *100" | bc -l`
#echo $r
# r=percent of mapped reads/reads after triming // scale: number after "," and bc -l used for float number
echo "Mappability in percent: $r %" >> /home/phumbert/raw_data/Chipseq/libsize/$f.libsize.txt
#save this value in the previous file (add information in a file ">>")

#Produce .bedGraph files and normalized reads, using bamCoverage
bamCoverage --bam $f.sorted.bam -o $f.bedGraph --binSize 1 --scaleFactor 0.001 --normalizeUsingRPKM --outFileFormat bedgraph -p 3

#Produce .bw files more compact
wigToBigWig $f.bedGraph ~/genome_index/dm6flybase/dm6.chrom.sizes $f.bedGraph.bw

done 