#! /bin/bash
#$ -cwd
#$ -l h_rt=100:00:00,h_data=16G,highp,arch=intel*
#$ -N vcf1b_filtering
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028

# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load python/2.7
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

#### parameters:
rundate=20181119 # date genotypes were called (vcf_20180806 includes capture 02)
# only want to get rid of egregious sites where hardly anyone is called, instead of being stringent across populations
noCallFrac=1.0 # instead of maxnocallfrac, have min DP requirements that were set in 1a. So doing away with this filter for now. (setting to 1.0) : maximum fraction of genotypes that can be "no call" (./.) ; previously was 0.9 for 20180806 gt filtering
snpNoCallFrac=0.2 # max frac of no-call genotypes allowed in snp file that you're going to use for PCA
# saving the stringent filtering for later!

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
indir=$SCRATCH/captures/vcfs/vcf_${rundate}
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
scriptdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/data_processing/variant_filtering

bespokeFilterScript=$scriptdir/filtering_bespokeFiltersAndChecks.py

outdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $outdir


######################## Remove bad individuals (low coverage) ###############
# These individuals have missingness that is 1SD above the mean value (/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/QC_Reports/missingnessPerInd/20180806/noCall_per_Ind_all_5_passingFilters_raw_variants.txt)
# calculated using this script: /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/generate_plots_and_tables/plotMissingGTsPerIndividual/plotMissingGenotypeCountsPerIndividual.R
# low coverage:


# However, I may want to use this part of the script to remove the RWAB individuals that pull away in the PCA
ind1=112_Elut_AL_AM_AM92022
ind2=128_Elut_AK_AF3630
ind3=131_Elut_BER_34
ind4=143_Elut_KUR_5original
ind5=144_Elut_KUR_1
ind6=145_Elut_KUR_18
ind7=150_Elut_AK_AF3631
ind8=151_Elut_AK_AF3714
ind9=152_Elut_AK_AF3720
ind10=153_Elut_AK_AF3181
ind11=154_Elut_KUR_11
ind12=155_Elut_KUR_9
ind13=161_Elut_MED_19
ind14=76_Elut_BER_31
ind15=86_Elut_MED_25
ind16=RWAB003_15_ELUT_CA_159
ind17=RWAB003_22_ELUT_CA_410


### also decided to remove all RWAB individuals (capture 2 that was highly multiplexed)
# because they are outliers at low frequencies

# keeping in relatives/admixed, because they aren't found til later in the pipeline; remove them in step 1f

badInds="$ind1 $ind2 $ind3 $ind4 $ind5 $ind6 $ind7 $ind8 $ind9 $ind10 $ind11 $ind12 $ind13 $ind14 $ind15 $ind16 $ind17 and all RWAB samples"

echo "starting step 6: remove bad individuals"
echo "These are the bad individuals that are getting removed: $badInds because they have > mean + 1sd missing genotypes"
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${outdir}/'all_5_passingFilters_'${infile} \
-o ${outdir}/'all_6_rmBadIndividuals_passingFilters_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
-xl_sn ${ind3} \
-xl_sn ${ind4} \
-xl_sn ${ind5} \
-xl_sn ${ind6} \
-xl_sn ${ind7} \
-xl_sn ${ind8} \
-xl_sn ${ind9} \
-xl_sn ${ind10} \
-xl_sn ${ind11} \
-xl_sn ${ind12} \
-xl_sn ${ind13} \
-xl_sn ${ind14} \
-xl_sn ${ind15} \
-xl_sn ${ind16} \
-xl_sn ${ind17} \
-xl_se 'RWAB003_.+_ELUT'

### make sure there are no spaces after the \ ! ###
# note that 4 extremely low samples were already excluded prior to genotype calling from the medgenome captures
# note I need to add the other bad individuals from 
echo "done with step 6: remove bad individuals"
# double check that it worked:
echo "checking that bad individuals are gone"
> ${outdir}/badIndsThatSurvivedFiltering.shouldBeEmpty.txt
for i in $badInds RWAB
do
zcat ${outdir}/'all_6_rmBadIndividuals_passingFilters_'${infile} | grep $i >> ${outdir}/badIndsThatSurvivedFiltering.shouldBeEmpty.txt
done

#################################################################################
############################ RUN BESPOKE FILTERS and UPDATE AN/AC ##########################
#################################################################################
# These filters will check for:
# usage: python [script] [full path to invcf] [full path to out vcf] [full path to error file] [max no call fraction]
# 1. ref or alt alleles that aren't a single letter (AGCT) or . (alt)
# 2. genotypes that aren't in 0/0, 0/1, 1/1 or ./. (maybe it's phased, etc)
# 3. must have qual score
# 4. must be PASS for site
# 5. make sure not missing DP, AN, GT, AD, DP or GQ/RGQ
# X. (does not do this anymore) make sure no called genotype is missing any info from the genotype info field [20180831: this now happens per-population]
# X. (does not do this anymore) gets rid of sites where all calls are 0/1 (all hets)  [20180831: this now happens per-population]
# 8. updates AN and AC based on final sets of calls (these aren't updated when GATK does genotype filtering)
# 9. filters sites that exceed the maximum no-call fraction (supplied by user)
# this script does NOT: change any genotypes; do any genotype filtering; change any FT fields for genotypes (./. gts will still be PASS if they started as ./. -- bit of GATK weirdness that isn't fatal)
echo "starting step 7a: carrying out bespoke filtering (includes filtering out sites exceeding max no-call frac: $noCallFrac)"

python $bespokeFilterScript ${outdir}/'all_6_rmBadIndividuals_passingFilters_'${infile} \
${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
${outdir}/'fail_all_7_FAILINGBespoke_passingFilters_'${infile%.vcf.gz}.txt \
$noCallFrac
# bgzip the result: (note: must use bgzip not gzip)
$bgzip -f  ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz}
$tabix -f -p vcf ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} # index the vcf

echo "done with step 7a: carrying out bespoke filtering"
# you should gzip your bad sites eventually.
#################################################################################
############################ Get sets of variant / invariant sites for initial PCA ###############
############################ Note that you're only going to use these for relatedness/admixture detection; ###############
############################ Eventually going to remake them without your bad individuals ###############
############################ have a missingness filter of 0.2 ###############
#################################################################################
echo "starting step 7b: select final passing snps from merged file"
# Some sites that may have started as variant may have become INVARIANT by the end.
# Want these to end up in the INVARIANT category. 
# select the biallelic snps: 
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
-V ${outdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${outdir}/'snp_7_maxNoCallFrac_'${snpNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${snpNoCallFrac}
echo "done step 7b: select final passing snps from merged file"

