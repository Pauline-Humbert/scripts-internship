#!/bin/bash

cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/

di=cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/

for DName in `ls $di` # direcory in motif analysis loop to create for each TF a txt file with the 2 motifs (elba and insv) like the model 
do 
d=$DName
cd /home/phumbert/raw_data/Chipseq/Motif_Analysis/$d/meme_out/ #go in this directory to access meme.txt

f=/home/phumbert/raw_data/Chipseq/Motif_Analysis/$d/1_PWM_motif/matrix_2motifs.txt
echo ${d}

#copy the 4th line and remove the end of the line 
sed -n '4p' meme.txt | sed -e 's/(Release date: Tue Jun 27 16:22:50 2017 -0700)$//' > $f 
echo `\t` >> $f 

#copy the 32th and 305 line: alphabet and strand
sed -n '32p' meme.txt >> $f
echo `\t` >> $f

sed -n '305p' meme.txt >> $f
echo `\t` >> $f

#copy informatio about background letter frequencies
echo "Background letter frequencies" >> $f
sed -n '310p' meme.txt >> $f 
echo `\t` >> $f

if [ "$d" = "wtElba1_elba1Elba1" ]
then
	### matrix motif insv
	#copy the line without space at the beginning of the line and without the end of the line
	motif1=`sed -n '1274p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-4 position-specific probability matrix$//'`
	#copy the matrix
	echo "MOTIF $motif1" >> $f
	sed -n '1276,1288p' meme.txt >> $f
	echo `\t` >> $f
	echo `\t` >> $f

	### matrix motif elba
	motif2=`sed -n '1760p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-6 position-specific probability matrix$//'` 
	echo "MOTIF $motif2" >> $f
	sed -n '1762,1771p' meme.txt >> $f

elif [ "$d" = "wtElba2_elba2Elba2" ]
then	
	### matrix motif insv
	motif1=`sed -n '533p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-1 position-specific probability matrix$//'`
	echo "MOTIF $motif1" >> $f 
	sed -n '535,547p' meme.txt >> $f
	echo `\t` >> $f
	echo `\t` >> $f
	
	### matrix motif elba
	motif2=`sed -n '1531p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-5 position-specific probability matrix$//'`
	echo "MOTIF $motif2" >> $f 
	sed -n '1533,1542p' meme.txt >> $f

elif [ "$d" = "wtElba3_elba3Elba3" ]
then	
	### matrix motif insv
	motif1=`sed -n '1262p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-4 position-specific probability matrix$//'`
	echo "MOTIF $motif1" >> $f 
	sed -n '1264,1276p' meme.txt >> $f
	echo `\t` >> $f
	echo `\t` >> $f
	
	### matrix motif elba
	motif1=`sed -n '2212p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-8 position-specific probability matrix$//'`
	echo "MOTIF $motif2" >> $f 
	sed -n '2214,2223p' meme.txt >> $f

elif [ "$d" = "wtInsv_insvInsv" ]
then	
	### matrix motif insv
	motif1=`sed -n '1039p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-3 position-specific probability matrix$//'`
	echo "MOTIF $motif1" >> $f 
	sed -n '1041,1053p' meme.txt >> $f
	echo `\t` >> $f
	echo `\t` >> $f

	### matrix motif elba
	motif2=`sed -n '1534p' meme.txt | sed -e 's/^[ \t Motif]*//' | sed -e 's/MEME-5 position-specific probability matrix$//'`
	echo "MOTIF $motif2" >> $f 
	sed -n '1536,1548p' meme.txt >> $f

else
	echo ""
fi
done 