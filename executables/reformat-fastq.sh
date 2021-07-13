#! /bin/bash

###################
# Reformattingmetagenomic reads with BBTools to take paired end files to interleaved
# For use on WEI GLBRC servers running HT Condor
# Elizabeth McDaniel 
##################

# set path where fastp is installed in local home directory bin
BBPATH=/opt/bifxapps/bbmap-38.32/

# queueing r1 r2 metagenomic reads and output folder/file names
file1=$1
file2=$2
out=$3

$BBPATH/reformat.sh in=$file1 in2=$file2 out=$out