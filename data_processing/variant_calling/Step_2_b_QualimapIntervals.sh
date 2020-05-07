#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=3G,highp
#$ -pe shared 5
#$ -N qualimap
#$ -m abe
############ carry out qualimap

# file locations:

header=$1


SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
fastqs=$wd/fastqs
paleomixOutput=$wd/paleomix/
outdir=$wd/qualimap/$header
mkdir -p $outdir/

bed=$wd/coveredIntervals/${header}.coveredIntervals.bed
qmap=/u/home/a/ab08028/klohmueldata/annabel_data/bin/qualimap_v2.2.1/qualimap
REFPREFIX=Mustela_putorius_furo.MusPutFur1.0.dna.toplevel
#modules
source /u/local/Modules/default/init/modules.sh
module load R
module load java

$qmap bamqc --java-mem-size=10G -nt 5 -bam $paleomixOutput/$header/${header}.${REFPREFIX}.bam -gff $bed -outdir $outdir -outfile ${header}.qualimap.out.html -outformat HTML
 
sleep 10m
