#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments

genome=$1
mapping=$2

samplename=$(basename $mapping -spRep.sorted.bam)
genomeName=$(basename $genome .fa)
resultName=$genomeName-vs-$samplename.IS

fasta=/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/ref_genomes/$genomeName.fa
genes=/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/ref_genomes/$genomeName.genes.fna

# cd to mapping folder

cd /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/mappingResults

# profile command

inStrain profile $mapping $fasta -o /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/inStrain/$resultName -p 8 -g $genes -s /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/binningResults/relative_abundance/Abigail-scaffolds-to-bins.tsv