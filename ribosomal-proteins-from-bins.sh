#!/bin/bash
#SBATCH -c 12                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=50gb                         # mem in gb
#SBATCH -t 14-0:0:0                         # How long will your job run for? If omitted, the default is 3 hours.
#SBATCH -J ribosomal_p                      	# Name of job
#SBATCH -e /home/gcolb087/scratch/ribosomal_p-muscle.err
#SBATCH -o /home/gcolb087/scratch/ribosomal_p-muscle.out

BASEDIR=/home/gcolb087/scratch
OUTPUT=$BASEDIR/output_anvio5_1_Sept6
TREE=$OUTPUT/phylogenetic_tree
BINS=$OUTPUT/SAMPLES-SUMMARY-concoct-checkm
NCBI=$BASEDIR/ncbi_genomes
TRANS=/home/gcolb087/bin

module load muscle/3.8.31
module load mafft/7.310

####MUSCLE
#fetch the ribosomal proteins for the bins, align them with translatorX and concatenate alignments
mkdir -p $BINS/ribosomal_proteins/
mkdir $BINS/ribosomal_proteins/muscle

for bin in $(awk '{print $0}' $BINS/goodbins.names)
do for annotation in `awk '{print $1}' $TREE/RP_annotations.txt`
do grep $annotation $BINS/bin_by_bin/${bin}/${bin}-gene_calls.txt | awk -v OFS="" -v x=${annotation} -v y=${bin} '{print ">", x, "_", y, "\n", $NF}' > $BINS/ribosomal_proteins/${bin}_${annotation}.fa &
done
wait
done

#concatenate ribosomal fasta and bins fastsa together 
for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do
(cat $NCBI/ribosomal_proteins/${annotation}.fa $BINS/ribosomal_proteins/*_${annotation}.fa > $BINS/ribosomal_proteins/translatorx/muscle/${annotation}.fa
perl $TRANS/translatorx_vLocal.pl -i $BINS/ribosomal_proteins/translatorx/muscle/${annotation}.fa -o $BINS/ribosomal_proteins/translatorx/muscle/${annotation}_muscle_translatorx -c 11 ) &
done
wait

for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do perl $TRANS/translatorx_vLocal.pl -i $BINS/ribosomal_proteins/translatorx/muscle/${annotation}.fa -o $BINS/ribosomal_proteins/translatorx/muscle/${annotation}_muscle_translatorx -p M -c 11 -t T
done
wait


		##############
		## muscle-reading-frame
		#############
for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do cat $NCBI/ribosomal_proteins/${annotation}.fa $BINS/ribosomal_proteins/*_${annotation}.fa > $BINS/ribosomal_proteins/translatorx/muscle-reading-frame/${annotation}.fa
done
wait

for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do perl $TRANS/translatorx_vLocal.pl -i $BINS/ribosomal_proteins/translatorx/muscle-reading-frame/${annotation}.fa -o $BINS/ribosomal_proteins/translatorx/muscle-reading-frame/${annotation}_muscle_translatorx -p M -c 11 -t T
done
wait








###### MAFFT
#fetch the ribosomal proteins for the bins, align them with translatorX and concatenate alignments
#concatenate ribosomal fasta and bins fastsa together 
for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do
(cat $NCBI/ribosomal_proteins/${annotation}.fa $BINS/ribosomal_proteins/*_${annotation}.fa > $BINS/ribosomal_proteins/mafft/${annotation}.fa
perl $TRANS/translatorx_vLocal.pl -p F -i $BINS/ribosomal_proteins/mafft/${annotation}.fa -o $BINS/ribosomal_proteins/mafft/${annotation}_mafft_translatorx -c 11 ) &
done
wait


###### prank 
mkdir -p $BINS/ribosomal_proteins/prank
#concatenate ribosomal fasta and bins fastsa together 
for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do
(cat $NCBI/ribosomal_proteins/${annotation}.fa $BINS/ribosomal_proteins/*_${annotation}.fa > $BINS/ribosomal_proteins/prank/${annotation}.fa
perl $TRANS/translatorx_vLocal.pl -p P -i $BINS/ribosomal_proteins/prank/${annotation}.fa -o $BINS/ribosomal_proteins/prank/${annotation}_prank_translatorx -c 11 ) &
done
wait


###########
#this might just be another way to do the above combine ncbi with ribosomal proteins from the bins
# for annontation in RP_annotations.txt copy fasta file over, then append Bin_#_#_annotation in to annotation file 
for annotation in $(awk '{print $1}' $TREE/RP_annotations.txt)
do 
cp $NCBI/ribosomal_proteins/${annotation}.fa $BINS/ribosomal_proteins/${annotation}.fa
cat $BINS/ribosomal_proteins/*_${annotation}.fa >> $BINS/ribosomal_proteins/${annotation}.fa
done
###########

