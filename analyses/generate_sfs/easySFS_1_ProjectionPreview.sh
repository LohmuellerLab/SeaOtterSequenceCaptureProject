#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=32G,highp
#$ -N easySFSPreview
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
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/neutral_and_cds_VCFs/neutralVCFs
popFile=/u/flashscratch/a/ab08028/captures/samples/samplesPop.Headers.forEasySFS.3.20181119.txt # this doesn't have baja on it; doesn't have any admixed/bad inds on it. 
# this has admixed in it , but they aren't in pop file

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=${gitdir}/analyses/generate_sfs/

easySFS=$scriptdir/easySFS.abModified.3.noInteract.Exclude01Sites.HetFiltering.20181121.py  # this is my modification
outdir=/u/flashscratch/a/ab08028/captures/analyses/SFS/$genotypeDate/easySFS/projection_preview
mkdir -p $outdir
# had to modify easySFS so that it wouldn't prompt a "yes/no" response about samples that are missing from VCF file
# want it to just output that info and continue , not prompt yes/no.
# this vcf has all snps across all categories (cds, neutral, etc.) with 0.9 max no call frac (v. liberal)
# and has had all individuals removed that won't go into the SFS
# going to do the actual projection for each category of site
vcf=neutral.snp_9b_forEasySFS_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf

$easySFS -i $vcfdir/${vcf} -p $popFile --preview -a -v > $outdir/neutral.snp_9b.easySFS.projPreview.txt


### now plot projections in R and decide on your levels. Actually DO the projections on 

