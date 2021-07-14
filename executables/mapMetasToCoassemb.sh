#! /bin/bash

# Map metagenomic reads to refs/assemblies

export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

BBPATH=/opt/bifxapps/bbmap-38.32

# cd to directory of metagenomes
cd /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes

ref=$1
meta=$2
out=$3

$BBPATH/bbmap.sh ref=$ref in=$meta outm=$out idtag minid=0.95 nodisk -Xmx50g
samtools sort $out -o ${out%.bam}.sorted.bam
