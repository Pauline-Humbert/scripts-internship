#!/bin/bash


# Before running this script: create a "a meme environnement" to run "fimo" command in "Motif Analysis"                                                  
# Create this environnement with these following commands (in bash)
#conda create -n env_meme
#source activate env_meme

cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/

di=cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/

for DName in `ls $di`
do
d=$DName

cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/${d}/1_PWM_motif/

fimo --parse-genomic-coord --oc fimo matrix_2motifs.txt ${d}"_"peaks"_"narrowPeak.fa
#create a file to know what peaks contain a motif

done