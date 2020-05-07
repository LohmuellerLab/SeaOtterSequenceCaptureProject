#! /bin/bash
#$ -cwd
#$ -l h_rt=72:00:00,h_data=2G,highp,h_vmem=36G
#$ -m abe
#$ -M ab08028
#$ -pe shared 16
#$ -N angsdStep1aiii
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd
################# 1e-06 snps only -- transversions only ##############
######### Call GLs and GPs, counts and mafs using ANGSD based on full coverage modern + aDNA mapped to elut/mfur **using allele freqs as prior for GPS** ######

#### run specific settings ####
trimValue=7 # set value you want to trim from either end of read (looking at mapdamage plots)
posterior=1 # setting for angsd -doPost : 1 for using allele frequencies as prior, 2 for using a uniform prior 
snpCutoff=1e-06
#todaysdate=`date +%Y%m%d`'-highcov-AFprior-MajorMinor4'
#todaysdate=`date +%Y%m%d`'-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL'
todaysdate="20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL-RedoneToReplaceDeletedFiles"
#### ANGSD v 0.923 ####
source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env

######### dirs and files ###########
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling_aDNA
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/aDNA-ModernComparison
bamdir=$wd/bams/
GLdir=$wd/angsd-GLs
mkdir -p $GLdir
mkdir -p $GLdir/$todaysdate
outdir=$GLdir/$todaysdate

### list of bam files to include: high coverage modern + aDNA: includes COM AL KUR now
elutBamList=$scriptDir/bamLists/angsd.bamList.mappedtoElutfullpaths.HighCovPlusADNA.PlusCOM.KUR.AL.txt # list of bam files mapped to sea otter, including downsampled AND non-downsampled
mfurBamList=$scriptDir/bamLists/angsd.bamList.mappedtoMfurfullpaths.HighCovPlusADNA.PlusCOM.KUR.AL.txt  # list of bam files mapped to ferret, including downsampled AND non-downsampled

# reference genomes:
elutRef=/u/home/a/ab08028/klohmueldata/annabel_data/sea_otter_genome/dedup_99_indexed_USETHIS/sea_otter_23May2016_bS9RH.deduped.99.fasta
mfurRef=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta

echo -e "THIS USES HIGH COVERAGE MODERN + ANCIENT ONLY\nBamLists used:\n$elutBamList\n$mfurBamList \ntrimvalue = $trimValue\ndoPost posterior setting = $posterior (1 = use allele freq as prior; 2 = use uniform prior)" > $GLdir/$todaysdate/HIGHCOVERAGEONLY.txt
echo -e "NOTE: lost my original elut runs due to SCRATCH backup failing, so I'm redoing it now for the record. But all the PCA and admixture stuff is done and fine based on 20191212 runs. This is just a redo."

######### ANGSD settings:##############

####### Mfur mapped bams ############
spp="mfur"
ref=$mfurRef
bamList=$mfurBamList

angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 4 -doMaf 1 \
-beagleProb 1 -doPost $posterior \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
-out $outdir/angsdOut.mappedTo${spp}.${snpCutoff}.snpsOnly.TransvOnly \
-doGlf 2 \
-doCounts 1 -dumpCounts 2 -doDepth 1 \
-SNP_pval $snpCutoff \
-rmTrans 1



####################### elut ###############
####### Elut mapped bams ############
spp="elut"
ref=$elutRef
bamList=$elutBamList

angsd -nThreads 16 \
-ref $ref \
-bam $bamList \
-GL 2 \
-doMajorMinor 4 -doMaf 1 \
-beagleProb 1 -doPost $posterior \
-remove_bads 1 -uniqueOnly 1 \
-C 50 -baq 1 -trim $trimValue -minQ 20 -minMapQ 25 \
-out $outdir/angsdOut.mappedTo${spp}.${snpCutoff}.snpsOnly.TransvOnly \
-doGlf 2 \
-doCounts 1 -dumpCounts 2 -doDepth 1 \
-SNP_pval $snpCutoff \
-rmTrans 1

