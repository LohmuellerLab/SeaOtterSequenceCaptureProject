#! /bin/bash
#$ -cwd
#$ -l h_rt=30:00:00,h_data=16G,highp
#$ -N vcf1d_removeRelatives
#$ -o /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -e /u/flashscratch/a/ab08028/captures/reports/GATK
#$ -m abe
#$ -M ab08028
# Realized that maxnocallfrac doesn't act on the individuals you've selected, but rather the whole original VCF!
# so was being far too stringent! 
# So what I do now is I make a population vcf without a nocallfraction set - it will just have the baseline maxnocallfrac 0.9
# that was in the original file it is being drawn from. This might be useful for projecting down in dadi.
# THEN I am going to filter the allCalled
######## These indivduals were identified in PLINK as having a kinship coeff > 0.2 with another individual in the same population
### you need to manually enter which individuals to remove; this is not an automated script; so be careful!
# I am removing one from each pair 
# modules
source /u/local/Modules/default/init/modules.sh
module load java
module load python
GATK=/u/home/a/ab08028/klohmueldata/annabel_data/bin/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar
tabix=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/tabix
bgzip=/u/home/a/ab08028/klohmueldata/annabel_data/bin/tabix-0.2.6/bgzip

#### parameters:
rundate=20181119 # date genotypes were called (vcf_20180806 includes capture 02)
noCallFrac=1.0 # maximum fraction of genotypes that can be "no call" (./.) that was used in previous steps in previous "all sites" files (currently no cutoff)
snpNoCallFrac=0.2 # max missingness allowed in snp file (stricter cutoff)
perPopNoCallFrac=0 # max missingness allowed in final file for each pop for sfs (super strict cutoff)
maxHetFilter=0.75 # (maximum heterozygous genotypes per site filter)
#### file locations
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering
mkdir -p $wd
infile=raw_variants.vcf.gz ### make sure this doesn't have a path as part of its name! just infile names
REFERENCE=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
vcfdir=$wd/${rundate}_filtered # date you called genotypes
mkdir -p $vcfdir/populationVCFs
mkdir -p $vcfdir/populationVCFs/admixedVCFs
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/data_processing/variant_filtering
hetFilterScript=filtering_perPopulation.noCall.maxHetFilter.py

## Relatives to remove
ind1="RWAB003_19_ELUT_CA_352"
ind2="106_Elut_AL_AD_GE91109"
#ind3="101_Elut_KUR_3"
#ind4="102_Elut_KUR_4"
# switched ind3 and ind4 to remove other member of relative pair (removing RWAB) -- 20181204
#ind3="RWAB003_28_ELUT_KUR_12" # already removed RWAB samples
#ind4="RWAB003_26_ELUT_KUR_15"  # already removed RWAB samples
ind5="104_Elut_KUR_7"
ind6="77_Elut_KUR_14"

### Heavily admixed invididuals / PCA outliers to keep out of population vcfs (can stay in all inds vcf):
ind7="90_Elut_KUR_20"
ind8="80_Elut_KUR_19"
ind9="91_Elut_KUR_21"
ind10="92_Elut_KUR_22"
ind11="93_Elut_KUR_23"
ind12="68_Elut_AK_AF3370"
ind13="56_Elut_AK_AF24903"
ind14="57_Elut_AK_AF24906"
ind15="163_Elut_AK_AF24915"
## mildly admixed invididuals / PCA outliers to remove:
ind16="125_Elut_AK_AF3369" # (adding)
ind17="65_Elut_AK_GE91060" # (adding)
ind18="165_Elut_AK_AL4661" # (doesn't appear admixed but is PCA outlier; so am keeping out of pop file and out of admixed file)
ind19="141_Elut_CA_419" # doesn't appear admixed but is a PCA outlier 
###### ** UPDATED THIS LIST *** for 20181119 ## 

## 20180910 update: I am not going to remove admixed individuals from the all_8 vcf files
# Because they are valid sequences for PCA, etc. Just removing them from pop-specific files

############################# remove admixed and relatives and other outliers from the all file ###############
# update in 20181206: not making an all file that removes relatives but keeps admixed; that is excessive file size
# and you can get that from all_7 if you need it (just exclude relatives)
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/'all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-xl_sn ${ind1} \
-xl_sn ${ind2} \
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
-xl_sn ${ind18} \
-xl_sn ${ind19} \
-trimAlternates
# adding trim alternates? 20181214


######## removing all outliers, relatives and admixed here ########

#######################################################################################
###################################### imposing excess heterozygosity filter ################################
#######################################################################################

##### going to filter out sites that exceed a number of het-calls across all samples (then I do more specific population filtering)
# of heterozygosity within easysfs. But the goal of this stage is to make sure that my files that go into PCA, etc. aren't too affected by these
# weird sites.
# don't use no maxnocallfrac filter at this pont (though the script is capable)(comes at next stage)

python $scriptdir/$hetFilterScript --vcf ${vcfdir}/'all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
--outfile ${vcfdir}/'all_9_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} \
--errorfile ${vcfdir}/'sitesFailingMaxHetFilter_'${maxHetFilter}'.txt' \
--maxNoCallFrac $noCallFrac \
--maxHetFilter $maxHetFilter

$bgzip -f ${vcfdir}/'all_9_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile%.gz} 
$tabix -f -p vcf ${vcfdir}/'all_9_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} 


#######################################################################################
######################### after pop-vcfs made, make a final version  #########################
######################### of the all sites vcfs that uses a 20% cutoff #######################
#######################################################################################
######## snps file: this one has a 0.2 missingness cutoff for use in PCA, etc.:
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
--variant ${vcfdir}/'all_9_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile}  \
-trimAlternates \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${vcfdir}/'snp_9a_forPCAetc_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${snpNoCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--maxNOCALLfraction ${snpNoCallFrac}  \
--excludeNonVariants
# this call to maxnocallfrac is okay because it's for the snp file with a 20% cutoff for use in pca, etc.
########### this also keeps PCA outliers.

######### snp file: this one has no missingness cutoff for use in EASY SFS -- will become neutral regions.
java -jar -Xmx4G ${GATK} \
-T SelectVariants \
-R ${REFERENCE} \
--variant ${vcfdir}/'all_9_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile}  \
-trimAlternates \
--restrictAllelesTo BIALLELIC \
--selectTypeToInclude SNP \
-o ${vcfdir}/'snp_9b_forEasySFS_maxHetFilter_'${maxHetFilter}_'rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_passingBespoke_passingAllFilters_postMerge_'${infile} \
--excludeNonVariants
# no maxnocall frac.
# added trimAlternates in case removal of some individuals made some sites not variable anymore.
# added excludeNonVariants as extra precaution since I was seeing some sites that used to be variable not getting removed by trim alt
# need to do trim alt upstream (all 7 --> all 8 )
#### CHECK final output to make sure there are no sneaky sites that are just made up of 0/0s and ./.s and so aren't really variant anymore 

#######################################################################################
############## put admixed individuals (but not relatives) into their own VCFs for later ##############
######### note these have not been through het filtering #################
#######################################################################################
## KURIL ADMIXED (you take out of the all_7 file that still has admixed inds. in there and you positively select them with sn):
# 20181206: don't need to make all_8 that keeps admixed; can just pull them out of all_7:
java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
-trimAlternates \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/'admixIndOnly_KUR_all_8_notHetFiltered_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-sn ${ind7} \
-sn ${ind8} \
-sn ${ind9} \
-sn ${ind10} \
-sn ${ind11} 

## Alaska ADMIXED:

java -jar $GATK \
-R $REFERENCE \
-T SelectVariants \
-trimAlternates \
--variant ${vcfdir}/'all_7_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
-o ${vcfdir}/populationVCFs/admixedVCFs/'admixIndOnly_AK_all_8_notHetFiltered_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
-sn ${ind12} \
-sn ${ind13} \
-sn ${ind14} \
-sn ${ind15} \
-sn ${ind16} \
-sn ${ind17} 
# skipping ind 18 because it odesn't appear admixed in FASTRUCTURE; just appears like a PCA outlier


#######################################################################################
############### not making population VCFs anymore (because projecting with EasySFS from one file) ####################
#######################################################################################
# Bering/Medny --> Commanders (COM)
## java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/COM_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_BER_.+' \
# -se '.+_Elut_MED_.+'
# note: no more maxnocall frac at this stage; file will be big -- but in all_9 my script filters it anyway

# note I'm keeping the old final name so that it flows into other scripts easily.

## California --> CA (include the RWAB hiseq4000 samples)
#java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/CA_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_CA_.+' \
# -se 'RWAB003_.+_ELUT_CA_.+' 
# relative was removed in main file above, so don't need to remove it here. 
## Alaska --> AK : remove admixed and make pop specific vcf
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/AK_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_AK_.+' \
# -xl_sn ${ind12} \
# -xl_sn ${ind13} \
# -xl_sn ${ind14} \
# -xl_sn ${ind15} \
# -xl_sn ${ind16} \
# -xl_sn ${ind17} \
# -xl_sn ${ind18} \
# -xl_sn ${ind19}

## Aleutian --> AL
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/AL_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_AL_.+'



## Kuril --> KUR (include the RWAB hiseq4000 samples)
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/KUR_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_KUR_.+' \
# -se 'RWAB003_.+_ELUT_KUR_.+' \
# -xl_sn ${ind7} \
# -xl_sn ${ind8} \
# -xl_sn ${ind9} \
# -xl_sn ${ind10} \
# -xl_sn ${ind11} 

## Baja --> BAJ
# java -jar $GATK \
# -R $REFERENCE \
# -T SelectVariants \
# --variant ${vcfdir}/'all_8_rmRelatives_keepAdmixed_passingBespoke_maxNoCallFrac_'${noCallFrac}'_rmBadIndividuals_passingFilters_'${infile} \
# -o ${vcfdir}/populationVCFs/BAJ_'all_8_rmRelativesAdmixed_passingAllFilters_maxNoCallFrac_'${noCallFrac}'.vcf.gz' \
# -se '.+_Elut_BAJ_.+' 




