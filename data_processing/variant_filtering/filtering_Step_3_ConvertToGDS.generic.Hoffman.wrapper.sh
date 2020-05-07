############## wrapper for GDS conversion script #############

source /u/local/Modules/default/init/modules.sh
module load R
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/

scriptdir=${gitdir}/scripts/scripts_20180521/data_processing/variant_filtering/
script=filtering_Step_3_ConvertToGDS.generic.Hoffman.R

genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered
outdir=$vcfdir/gdsFormat/
mkdir -p $outdir
# list of gzipped vcfs you want to convert (can be any VCFs)
# also want to make it with downsampled commanders vcf:
# this script gets run multiple times; first time is for snp_7 because that lets you do an initial PCA
# to figure out which individuals to exclude; then can do it on snp_9a after all filtering for a final gds file
# based on the final set of individuals passing
# May also want to downsample populations, etc. which you can actually do in R as part of running your PCA --
# so don't have to downsample the COMs here, but will want to in R.
vcfs2gds="snp_9a_forPCAetc_maxHetFilter_0.75_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz" # round 2

for vcf in $vcfs2gds
do
Rscript $scriptdir/${script} --vcf $vcfdir/$vcf --outdir $outdir
done
# output file name will be input file name with .gds extension instead of .vcf.gz extension
