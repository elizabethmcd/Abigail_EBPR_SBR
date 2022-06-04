#! /bin/bash

# Map metagenomic reads to assemblies

# setup environment for samtools dependencies to work 
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

# program paths
BBPATH=/opt/bifxapps/bbmap-38.32
METABATDIR=/opt/bifxapps/metabat-2.12.1/

# Arguments and paths
ref=$1
meta=$2
out=$3
metagenome_path="/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/cleaned_fastqs"

read1=$(basename $meta)
read2=$(basename $meta _R1.qced.fastq)_R2.qced.fastq
samplename=$(basename $out .bam)
assembly_name=$(basename $ref .fasta)

# cd to directory
cd /home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/mappingResults

# mapping job
$BBPATH/bbmap.sh ref=$ref in=${metagenome_path}/${read1} in2=${metagenome_path}/${read2} outm=$out idtag minid=0.95 nodisk -Xmx100g

# samtools commands 
samtools sort  $samplename.bam -o $samplename.sorted.bam
samtools index $samplename.sorted.bam $samplename.sorted.bam.bai

# get coverage table with jgi_summarize_contigs 
$METABATDIR/jgi_summarize_bam_contig_depths --outputDepth ${assembly_name}-depth.txt ${assembly_name}*.sorted.bam


