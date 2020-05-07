#! /bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=20G
#$ -m bea

########### gatk find covered regions
# do this after mapping
# going to do this for each file and then can use it in variant calling.
# this will maximize called sites
# then can try to find neutral regions among these. 

# modules

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
module load samtools
module load bedtools
# header:
header=$1 # input header into file
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
paleomixOutput=$wd/paleomix/${header} # specific to this header
outdir=$wd/coveredIntervals
mkdir -p $outdir
# programs
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
# ferret reference:
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel

# find covered intervals:
# starting out lax on coverage, then will filter out later.
# note. for modern use cov = 5 to get intervals; for ancient use cov = 1. 

java -jar $GATK \
	-T FindCoveredIntervals \
	-R $REFERENCE \
	-I $paleomixOutput/${header}.${REFPREFIX}.bam \
	-cov 5 \
	-minBQ 20 \
	-minMQ 30 \
	-o $outdir/${header}.coveredIntervals.list
########## Notice: this already gives you non-overlapping intervals so you should NOT do bedtools merge

# you'll then use this as -L when you call variants.
# minimum map quality : 30 (JAR pipeline)
# minimum base quality: 20 (JAR pipeline)


# Also want to make a bed file version
# how to: the interval list has format: 
# GL896898.1:37826-37995 (scaffold:start-stop) and is 1-based.
# So for .bed file, I need to make it 0-based on the start coord and non-inclusive for end coord, so keep the same
# so want scaffold\tstart-1\tstop for bed version (to use for other things)
########## also make a bed version: 
awk -F [:-] '{OFS="\t"; print $1,$2-1,$3}' $outdir/${header}.coveredIntervals.list > $outdir/${header}.coveredIntervals.int.bed
# want to check that nothing got made negative 1 (if it started at 0)
sed -i'' 's/-1/0/g' $outdir/${header}.coveredIntervals.int.bed
# update: 20180730: DO NOT MERGE. Nothing gets merged and it is just slow and sometimes errors out.
# and want to add dots for empty six fields 
awk '{OFS="\t"; print $1,$2,$3,".",".","."}' $outdir/${header}.coveredIntervals.int.bed > $outdir/${header}.coveredIntervals.bed
rm $outdir/${header}.coveredIntervals.int.bed
sleep 10m
