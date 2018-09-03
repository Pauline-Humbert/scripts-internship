#!/bin/bash

cd /home/phumbert/raw_data/Chipseq/MACS2results/

di=cd /home/phumbert/raw_data/Chipseq/MACS2results/

for DName in `ls $di` 
do 
d=$DName

cd /home/phumbert/raw_data/Chipseq/MACS2results/${d}/

# create a fasta file with all peaks (.narrowPeaks file) for fimo
bedtools getfasta -name -fi /home/phumbert/genome_index/dm6.fa \
-bed /home/phumbert/raw_data/Chipseq/MACS2results/${d}/${d}"_"peaks.narrowPeak \
-fo /home/phumbert/raw_data/Chipseq/MACS2results/${d}/${d}"_"peaks"_"narrowPeak.fa

# I duplicated peaks_narrowPeak.fa in /Motif_Analysis/$d/1_PWM_motif/
# cf 8_fasta_peak_id.r to obtain this final fasta file

done



