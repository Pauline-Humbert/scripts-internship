#!/bin/bash

cd /home/phumbert/raw_data/Chipseq/ #put directly in the correct directory
 
for FName in `ls [^wt]*.sorted.bam` # loop to run different commands gene by gene
do
f=`basename $FName .sorted.bam` # retrieve only the gene's name

# Identify areas in the genome that have been enriched with aligned reads:

if [ "$f" = "elba1Elba1" ]
then
	w=wtElba1
elif [ "$f" = "elba2Elba2" ]
then
	w=wtElba2
elif [ "$f" = "elba3Elba3" ]
then
	w=wtElba3
elif [ "$f" = "insvInsv" ]
then
	w=wtInsv
fi
# statement to define $w following $f

name=$w"_"$f 
#define the variable name => name of input

macs2 callpeak -t $w.sorted.bam  -c $f.sorted.bam --format=BAM  --gsize=dm  --qvalue=0.1 --nomodel --extsize 200 --cutoff-analysis --bdg --call-summits --outdir=/home/phumbert/raw_data/Chipseq/MACS2results/$name --name $name
#for identifying transcript factor binding sites/MACS captures the influence of genome complexity to evaluate the significance of enriched ChIP regions
#main MACS2: Call peaks from alignment results

annotatePeaks.pl /home/phumbert/raw_data/Chipseq/MACS2results/$name/$name"_"peaks.narrowPeak dm6  -go GO_dir -genomeOntology GenomeOntology_dir -annStats /home/phumbert/raw_data/Chipseq/MACS2results/$name/$name"_"stat.txt > /home/phumbert/raw_data/Chipseq/MACS2results/$name/$name"_"peakannofile.txt
#genomic feature association analysis (Genome Ontology), associate peaks with gene expression data, calculate ChIP-Seq Tag densities from different experiments, and find motif occurrences in peaks
done