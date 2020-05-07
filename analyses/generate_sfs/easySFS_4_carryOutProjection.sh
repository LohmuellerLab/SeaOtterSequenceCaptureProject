#! /bin/bash
#$ -cwd
#$ -l h_rt=10:00:00,h_data=20G,highp
#$ -N easySFSProjection
#$ -o /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -e /u/flashscratch/a/ab08028/captures/reports/SFS
#$ -m abe
#$ -M ab08028

####### Easy SFS
# https://github.com/isaacovercast/easySFS
# install:
# git clone git@github.com:isaacovercast/easySFS.git
# cd easySFS
# chmod +x *.py
# easySFS.py
source /u/local/Modules/default/init/modules.sh
module load python/2.7
module load R

bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
todaysdate=`date +%Y%m%d`
genotypeDate=20181119
noCallFrac=1.0
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt
allSamplesHetFilter=0.75 # het filtering done across all samples
perPopHetFilter=0.75 # the vcf file has already had some degree of maxHetFiltering. now per-population, easy sfs will do per-population filtering at this level  
# after troubleshooting and looking at HWE exact test, 0.75 was chosen.

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/


easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
# this version of script filters het sites that are excessively heterozygous 

## choose your projections: choosing for now: 
# CA,AK,AL,COM,KUR 
# choosing projections to maximize snps:
#Choices based on 20181220 projection preview: 
#AK Ð 14 (max)
#AL Ð 20 (max)
#CA Ð 12 (max would be at 8)
#COM Ð 34 (max)
#KUR Ð 12 (max would be at 10)
populations="CA,AK,AL,COM,KUR"
projections="12,14,20,34,12" # haploids; updated these values on 20181220

### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
# write out projection values:

############################### NEUTRAL SITE PROJECTIONS #############################

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/neutral/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# write projection choices into a readme
# this gets written out as part of monomorphic script now; don't need to do here.
# make sure vcf isn't zipped

allVCF=neutral.all_9_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
snpVCF=neutral.snp_9b_forEasySFS_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf


### NOTE: projection values must be in same order as populations are in your popFile (this isn't ideal -- at some point I am going to modify the easySFS script)
# note that order is CA,AK,AL,COM,KUR 
$easySFS -i $vcfdir/neutralVCFs/${snpVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter
# test run only had --proj 20 (maybe just projects 1 population? will have to see)
# -f forces overwrite of outdir
########## get counts of monomorphic sites to add to the SFSes ############
# will output two files, one for individual population's counst and one file with the pairs of populations 
# for the pairs of populations, it is a count of sites where both populations were monomorphic at whatever their 
# individual projection levels were. Note that easy sfs takes care of cases where one pop is monomorphic and the other is polymorphic
# this is just about sites where all individuals in all pops is monomorphic but where there might be missing data
# for example: a monomorphic site that is 0/0 for all individuals. Alaska has > Projection value gts at that site, so does KUR, so it would be included in the AK-KUR 2D sfs
# if at another site, AK has > proj value but KUR has < its proj value, that site wouldn't be included in the 2D sfs (though would be included in AK's 1D sfs)
### updated this script to also write out the projection values for use in fsc wrapper script.

python $scriptdir/getMonomorphicProjectionCounts.1D.2DSFS.py --vcf $vcfdir/neutralVCFs/${allVCF} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir

############## adding monomorphic sites to fsc and dadi SFSes and writing out totals #####################
 # this script add monomorphic sites to 0 bin of fsc sfses. (overkill in dadi because it gets masked, but does it for completeness ) 
Rscript $scriptdir/easySFS_addInMonomorphicSites.R --dataDir $outdir --popFile $popFile --class neutral # will write them out in your data dir in new directories
# need to modify this R script so that it also adds to 2D sfses.
################################################################################
##################################### coding sites #############################
################################################################################

cdsVCF=cds_all_9_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz
synVCF=syn_vep_cds_snp_9b_forEasySFS_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf
misVCF=missense_vep_cds_snp_9b_forEasySFS_maxHetFilter_${allSamplesHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/synonymous/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
$easySFS -i $vcfdir/cdsVCFs/${synVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter

outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/missense/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
$easySFS -i $vcfdir/cdsVCFs/${misVCF} -p $popFile -a -v --proj $projections -f -o $outdir -maxHetFilter $perPopHetFilter


########## get counts of monomorphic cds sites -- still have to think about how to scale these for syn/mis ############
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/cds/monomorphicCDSSites/projection-${todaysdate}-hetFilter-${perPopHetFilter}
mkdir -p $outdir
#python $scriptdir/getMonomorphicProjectionCounts.py --vcf $vcfdir/cdsVCFs/${cdsVCF} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir
python $scriptdir/getMonomorphicProjectionCounts.1D.2DSFS.py --vcf $vcfdir/cdsVCFs/${cdsVCF} --popMap $popFile --proj $projections --popIDs CA,AK,AL,COM,KUR --outdir $outdir



