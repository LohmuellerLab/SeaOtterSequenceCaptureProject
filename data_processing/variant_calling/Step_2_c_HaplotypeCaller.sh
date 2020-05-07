#! /bin/bash
#$ -cwd
#$ -l h_rt=150:00:00,h_data=21G,highp,arch=intel*
#$ -m abe
# use a submission script

# carry out haplotype caller

source /u/local/Modules/default/init/modules.sh
module load java/1.8.0_111
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel
# header:
header=$1


# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
paleomixOutput=$wd/paleomix/$header
outdir=$wd/gvcfs
mkdir -p $outdir

# these are covered intervals (result of previous step)
intervals=$wd/coveredIntervals
intervalFile=$intervals/${header}.coveredIntervals.list

java -jar $GATK \
-T HaplotypeCaller \
-R $REFERENCE \
-ERC BP_RESOLUTION \
-mbq 20 \
-L $intervalFile \
--interval_padding 100 \
-out_mode EMIT_ALL_SITES \
-I $paleomixOutput/${header}.${REFPREFIX}.bam \
-o $outdir/${header}.g.vcf.gz
# BP_RESOLUTION keeps all sites
# this just calls across covered targets (+- 100bp); much faster
sleep 10m
