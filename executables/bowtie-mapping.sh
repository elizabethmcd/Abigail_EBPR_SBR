#! /bin/bash

# setup environment for samtools dependencies to work 

export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

# arguments
reads=$1
samplename=$(basename $reads .fastq)

cd /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/mappingResults

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/binningResults/relative_abundance/bt2/Abigail-final-bins.fasta --interleaved $reads > $samplename-spRep.sam


# BAM, sort, index
samtools view -S -b  $samplename-spRep.sam >  $samplename-spRep.bam
samtools sort  $samplename-spRep.bam -o  $samplename-spRep.sorted.bam
samtools index $samplename-spRep.sorted.bam $samplename-spRep.sorted.bam.bai