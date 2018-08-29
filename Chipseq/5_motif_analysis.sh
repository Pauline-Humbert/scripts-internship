#!/bin/bash                                                                                                          

# Before running this script: create a "a meme environnement" to run "meme-chip" command in "Motif Analysis"
# Create this environnement with these following commands (in bash)
#conda create -n env_meme  
#source activate env_meme

# Loop to select .bed files and modify them for the motif analysis
di= cd /home/phumbert/raw_data/Chipseq/MACS2results/ #to be directly in the correct directory

for DName in `ls $di`
do
for FName in $DName # loop to run different commands gene by gene
do
d=$DName

f=`basename $FName .bed` 
# retrieve only the gene's name

# sort summits.bed file by peak score (5th column) and select only the 500 first rows
sort -r -k5 -n /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"summits.bed | head -500 > /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"summits"_sorted_top500".bed

# increase the size of each sequence in summit_sorted_top500.bed file: increase of 500 base pairs in each diection
slopBed -i /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"summits"_sorted_top500".bed \
-g /home/phumbert/genome_index/dm6flybase/dm6.chrom.sizes -b 500 > /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"slopBed"_"top500.bed

# create a fasta file for the motif analysis
bedtools getfasta -fi /home/phumbert/genome_index/dm6.fa \
-bed /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"slopBed"_"top500.bed \
-fo /home/phumbert/raw_data/Chipseq/MACS2results/$f/$f"_"sorted"_"slopBed"_"top500.fa

done
done

# Motif analysis 

# This is this part which need "meme environnement"

cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/ #to be in the directory where are "meme" results and where the "sorted_slopBed_top500.fa" was been copied

meme-chip -norand -ccut 100 -dna -meme-mod anr -nmeme 1000 -meme-minw 5 -meme-maxw 15 \
-meme-nmotifs 30 -dreme-e 0.05 -meme-maxsize 2000000 -oc wtElba1_elba1Elba1 wtElba1_elba1Elba1_sorted_slopBed_top500.fa

meme-chip -norand -ccut 100 -dna -meme-mod anr -nmeme 1000 -meme-minw 5 -meme-maxw 15 \
-meme-nmotifs 30 -dreme-e 0.05 -meme-maxsize 2000000 -oc wtElba2_elba2Elba2 wtElba2_elba2Elba2_sorted_slopBed_top500.fa

meme-chip -norand -ccut 100 -dna -meme-mod anr -nmeme 1000 -meme-minw 5 -meme-maxw 15 \
-meme-nmotifs 30 -dreme-e 0.05 -meme-maxsize 2000000 -oc wtElba3_elba3Elba3 wtElba3_elba3Elba3_sorted_slopBed_top500.fa

meme-chip -norand -ccut 100 -dna -meme-mod anr -nmeme 1000 -meme-minw 5 -meme-maxw 15 \
-meme-nmotifs 30 -dreme-e 0.05 -meme-maxsize 2000000 -oc wtInsv_insvInsv wtInsv_insvInsv_sorted_slopBed_top500.fa