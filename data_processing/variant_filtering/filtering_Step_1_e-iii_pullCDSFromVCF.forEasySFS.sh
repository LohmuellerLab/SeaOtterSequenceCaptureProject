#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=16G,highp
#$ -N exonVEPSFS
#$ -o /u/flashscratch/a/ab08028/captures/reports/VEP
#$ -e /u/flashscratch/a/ab08028/captures/reports/VEP
#$ -m abe
#$ -M ab08028
#$ -pe shared 3

# 20181018 decided this is an important step (Tanya's scripts for pi/total snps dont use chr designations so don't work for my data; remaking my scripts)
# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load bedtools
module load perl/5.10.1
module load htslib

GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
vepdir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/ensembl-vep/
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

rundate=20181119
noCallFrac=1.0 # no filter
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
cdsBed=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed
# based on MusPutFuro1.0.91.cdsOnly.gff
vcfdir=$wd/${rundate}_filtered # date you called genotypes
outdir=$vcfdir/neutral_and_cds_VCFs/cdsVCFs
mkdir -p $outdir
mkdir -p ${vcfdir}/bedCoords/cdsCallableSites/
allVCF=all_9_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters_raw_variants.vcf.gz

snpVCF=snp_9b_forEasySFS_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_passingBespoke_passingAllFilters_postMerge_raw_variants.vcf.gz

# intersect cds bed with population vcf files that have all sites called in all individuals

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/${snpVCF} \
-o $outdir/cds_${snpVCF} \
-L $cdsBed


# all vcf: doing this to get the total number of cds sites (will need to do something easy-sfsy to get total number for use in projections
# is this a big file?
# have to think this through a little more... What will my L be? Number of sites that are at the projection level for cds per population (with a multiplier for miss/syn)
# I think so. So will use my monomorphic python script on them. 
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/${allVCF} \
-o $outdir/cds_${allVCF} \
-L $cdsBed


### note: Jazlyn says vep works even with gzipped vcf (throws and error but actually works; confirm for myself)
# only run VEP on snp vcf file because those are the only things that will be annotated anyway!

# run VEP :
$vepdir/vep -v -i $outdir/cds_${snpVCF}  \
--cache --force_overwrite --species mustela_putorius_furo \
--variant_class --vcf --canonical \
-o $outdir/vep_cds_${snpVCF%.gz} \
--pick
# note vep output is not gzipped.
# Using PICK to pick one variant per site:
# https://uswest.ensembl.org/info/docs/tools/vep/script/vep_other.html#pick_options
# "This is the option we anticipate will be of use to most users. VEP chooses one block of annotation per variant, using an ordered set of criteria. This order may be customised using --pick_order. "
   # canonical status of transcript
   # APPRIS isoform annotation
   # transcript support level
   # biotype of transcript ("protein_coding" preferred)
   # CCDS status of transcript
   # consequence rank according to this table
   # translated, transcript or feature length (longer preferred)

$vepdir/filter_vep --filter "Consequence is synonymous_variant and CANONICAL is YES" \
--input_file $outdir/vep_cds_${snpVCF%.gz} \
--output_file $outdir/syn_vep_cds_${snpVCF%.gz} \
--force_overwrite

$vepdir/filter_vep --filter "Consequence is missense_variant and CANONICAL is YES" \
--input_file $outdir/vep_cds_${snpVCF%.gz} \
--output_file $outdir/missense_vep_cds_${snpVCF%.gz}  \
--force_overwrite 

# also do stopgained (added on 20190806)
$vepdir/filter_vep --filter "Consequence is stop_gained and CANONICAL is YES" \
--input_file $outdir/vep_cds_${snpVCF%.gz} \
--output_file $outdir/stopgained_vep_cds_${snpVCF%.gz}  \
--force_overwrite

# bgzip and tabix the VEP vcf 

# $bgzip -f $outdir/${pop}_cds_${suffix}.vcf
# $tabix -f -p vcf $outdir/${pop}_cds_${suffix}.vcf.gz
# 
# $bgzip -f $outdir/${pop}_VEP_cds_${suffix}.vcf
# $tabix -f -p vcf $outdir/${pop}_VEP_cds_${suffix}.vcf.gz
# 
# $bgzip -f  $outdir/${pop}_VEP_syn_${suffix}.vcf
# $tabix -f -p vcf $outdir/${pop}_VEP_syn_${suffix}.vcf.gz
# 
# $bgzip -f $outdir/${pop}_VEP_missense_${suffix}.vcf
# $tabix -f -p vcf $outdir/${pop}_VEP_missense_${suffix}.vcf.gz
# bgzip and tabix the VEP vcf
# okaky this works now. 

# get totals:
echo -e 'totalCalledcdsSites' >  ${vcfdir}/filteringStats/summary.cdsCallableSites.txt

######### For records, get a bed file of the coords that overlap between pop vcfs and neutBed, and the total amnt of neut sequence per pop: 

# for reference, want to get bed file of callable sites in cds regions : 
bedtools intersect -a $outdir/cds_${allVCF%.gz} -b $cdsBed | awk '{OFS="\t"; print $1,$2-1,$2}' | sort -k1,1 -k2,2n | bedtools merge -i stdin > ${vcfdir}/bedCoords/cdsCallableSites/${allVCF}.cdsOnly.callableSites.0based.bed
# and the total amount of coding sequence (cds) sequence:
totalSeq=`awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  ${vcfdir}/bedCoords/cdsCallableSites/${allVCF}.cdsOnly.callableSites.0based.bed`
echo -e ${pop}'\t'${totalNeut} >> ${vcfdir}/filteringStats/summary.cdsCallableSites.txt
# at some point get syn and mis totals.
# For CA, it's only 12 Mb called. -- 20180806

