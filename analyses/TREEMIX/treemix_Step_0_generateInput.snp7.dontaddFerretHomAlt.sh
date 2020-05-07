################ TREEMIX ############# 
# this should be run *after* you've converted your vcf to plink format
###### Involves some manual steps!!!!!!! ############ 

## It will make 3 versions of treemix input: 
module load treemix
module load plink

# converter script:
# wget https://bitbucket.org/nygcresearch/treemix/downloads/plink2treemix.py
# from treemix 3/12/12:
# Added a small script to convert stratified allele frequencies output from plink into TreeMix format. 
# This will be incorporated into the next release, but for the moment must be downloaded 
# separately. To run this, let's say you have data in plink format (e.g., data.bed, data.bim, 
# data.fam) and a plink cluster file matching each individual to a population (data.clust).
# data was formatted for PLINK to run faststructure (see those scripts)
## change third column of .fam to be the pop identifier, as save as population.clusters
# only want first three columns
#awk '{OFS="\t"; print $1,$2,$3}' population.clusters.temp > population.clusters
############## get plink format files from fastStructure_Step_1_vcf2plinkBed.20181119.sh in the FASTSTRUCTURE analysis section ##########

gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/
scriptdir=$gitdir/analyses/TREEMIX/
genotypeDate=20181119
vcfdir=/u/flashscratch/a/ab08028/captures/vcf_filtering/${genotypeDate}_filtered/
plinkFileDir=$vcfdir/plinkFormat/ 
treeFileDir=$vcfdir/treemixFormat/
mkdir -p $treeFileDir

#### DO FOR BOTH HEADERS (separate scripts) ##########

header=snp_7_maxNoCallFrac_0.2_passingBespoke_passingAllFilters_postMerge_raw_variants


# must add a snp name to my .bim file:
# make a backup
/bin/cp $plinkFileDir/$header.bim $plinkFileDir/${header}.original.bim
awk '{$2 = $1":"$4; print}' $plinkFileDir/${header}.original.bim >  $plinkFileDir/$header.bim
awk '{print $1,$2}' $plinkFileDir/$header.fam > $plinkFileDir/$header.samples
/bin/cp $plinkFileDir/$header.samples $plinkFileDir/$header.population.clusters 
clusters=$plinkFileDir/$header.population.clusters # 3 columns: 0 sampleID popID
# make a back up once you've added populations
/bin/cp $clusters $clusters.backup


##################################### With separate CA/BAJA and removing 3 relatives   ##############################
############# Excluding 3 Relatives: can use plink --filter ######################

# make a filter file:
# plink --file data --filter myfile.raw 1 --freq

#implies a file myfile.raw exists which has a similar format to phenotype and cluster files: that is, the first two columns are family and individual IDs; the third column is expected to be a numeric value (although the file can have more than 3 columns), and only individuals who have a value of 1 for this would be included in any subsequent analysis or file generation procedure. e.g. if myfile.raw were 
# copy cols 1, 2 and then a "1"

awk '{print $1,$2,1}' $plinkFileDir/$header.fam > $plinkFileDir/$header.exclList.allIndsIncluded
echo " YOU MUST MANUALLY EXCLUDE THE RELATIVES HERE" #### 
# excluding relatives from my excel sheet:
# 104_Elut_KUR_7 (keep 79_Elut_KUR_17)
# 77_Elut_KUR_14 (keep 81_Elut_KUR_2)
# 106_Elut_AL_AD_GE91109 (keep 118_Elut_AL_AD_GE91101)

###### MANUALLY edit and rename as : $plinkFileDir/$header.exclList.rmRelatives
# manually edit the file to set individuals you want to exclude to "0"
clusters=$plinkFileDir/$header.population.clusters.sepCA-BAJ # 3 columns: 0 sampleID popID
marker="sepCA-BAJ.exclRelatives"
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order \
--filter $plinkFileDir/$header.exclList.rmRelatives 1

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz

######################## EXCLUDE BAJA ########################################
############# Excluding 3 Relatives + BAJA : can use plink --filter ######################

# make a filter file:
# plink --file data --filter myfile.raw 1 --freq

#implies a file myfile.raw exists which has a similar format to phenotype and cluster files: that is, the first two columns are family and individual IDs; the third column is expected to be a numeric value (although the file can have more than 3 columns), and only individuals who have a value of 1 for this would be included in any subsequent analysis or file generation procedure. e.g. if myfile.raw were 
# copy cols 1, 2 and then a "1"

awk '{print $1,$2,1}' $plinkFileDir/$header.fam > $plinkFileDir/$header.exclList.allIndsIncluded
# excluding relatives from my excel sheet:
# 104_Elut_KUR_7 (keep 79_Elut_KUR_17)
# 77_Elut_KUR_14 (keep 81_Elut_KUR_2)
# 106_Elut_AL_AD_GE91109 (keep 118_Elut_AL_AD_GE91101)
# 168_Elut_BAJ_TS2
# 169_Elut_BAJ_R1
echo " YOU MUST MANUALLY EXCLUDE THE RELATIVES + BAJA  HERE" #### 

###### MANUALLY edit and rename as : $plinkFileDir/$header.exclList.rmRelatives.rmBAJA
# then manually edit the file to set individuals you want to exclude to "0"
clusters=$plinkFileDir/$header.population.clusters.sepCA-BAJ # 3 columns: 0 sampleID popID
# but it's confusing. the distinction between minor and non-ref is not well documented
marker=noBAJA.exclRelatives
plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order \
--filter $plinkFileDir/$header.exclList.rmRelatives.rmBAJA 1

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz



########### 20191007: new addition: merge Aleutian Islands ""****CombineAl****"" ####################
# manually edit clusters and want to exclude BAJ still ()
marker=noBAJA.exclRelatives.CombineAL
clusters=$plinkFileDir/$header.population.clusters.sepCA-BAJ.CombineAL # 3 columns: 0 sampleID popID

plink --bfile $plinkFileDir/$header \
--freq \
--missing \
--within $clusters \
--allow-extra-chr \
--out $plinkFileDir/$header.${marker} \
--nonfounders \
--keep-allele-order \
--filter $plinkFileDir/$header.exclList.rmRelatives.rmBAJA 1

gzip -f $plinkFileDir/${header}.${marker}.frq.strat

python $scriptdir/plink2treemix.py $plinkFileDir/${header}.${marker}.frq.strat.gz $treeFileDir/${header}.${marker}.frq.strat.treemixFormat.gz

