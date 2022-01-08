#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''


# arguments
sam=$1
outfolder=$(basename $sam .sorted.bam)-quick-profile

# cd to mapping results folder

cd /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/mappingResults

# inStrain quick profile command

inStrain quick_profile -p 2 -s ../binningResults/relative_abundance/abigail_new_bins_scaffolds_to_bins.tsv -o ../inStrain/quick_profiles/$outfolder $sam ../binningResults/relative_abundance/all-new-bins.fasta
