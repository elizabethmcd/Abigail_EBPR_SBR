#! /bin/bash

######################
# MetaBAT Binning on GLBRC
# Elizabeth McDaniel
######################

# Queue assembly files
ASSEMB=/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/assemblies/Abigail_coassembly_contigs.fasta
MAPDIR=/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/mappingResults
METABATDIR=/opt/bifxapps/metabat-2.12.1/
OUTBIN=/home/GLBRCORG/emcdaniel/EBPR/Abigail/metagenomes/binningResults

# Get depth 
$METABATDIR/jgi_summarize_bam_contig_depths --outputDepth $OUTBIN/Abigail-coassembly-depth.txt $MAPDIR/Abigail_coassembly_contigs*.sorted.bam

# Run metabat
$METABATDIR/metabat2 -i $ASSEMB -a $OUTBIN/Abigail-coassembly-depth.txt -o $OUTBIN/coassembly-bins/coassembly-bin
