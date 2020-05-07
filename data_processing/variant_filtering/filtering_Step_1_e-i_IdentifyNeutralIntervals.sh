################################# Set up ###########################
# modules
source /u/local/Modules/default/init/modules.sh
module load bedtools
module load blast

rundate=20181119
maxHetFilter=0.75 # het filter used across all samples (per population het filter occurs during easy sfs)

SCRATCH=/u/flashscratch/a/ab08028
filterDir=$SCRATCH/captures/vcf_filtering/${rundate}_filtered/
bedDir=$SCRATCH/captures/vcf_filtering/${rundate}_filtered/bedCoords
wd=$SCRATCH/captures/vcf_filtering/${rundate}_filtered/checkingNeutralSites
mkdir -p $wd
mfurDir=/u/home/a/ab08028/klohmueldata/annabel_data/ferret_genome
drerDir=/u/home/a/ab08028/klohmueldata/annabel_data/zebra_fish_genome
REFERENCE=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fasta
gff=$mfurDir/Mustela_putorius_furo.MusPutFur1.0.91.gff3
cdsRegions=$mfurDir/MusPutFuro1.0.91.cdsOnly.0based.sorted.merged.bed

exonicRegions=$mfurDir/MusPutFur1.0.91.exonsOnly.0based.sorted.bed # make this below if haven't already

mfurCpG=$mfurDir/CpG_Islands/cpgIslandExtUnmasked.already0based.bed 
mfurRepeat=$mfurDir/repeatMasker_UCSC/repeatRegions.0based.bed

# script to get gc content
getGC=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/scripts/scripts_20180521/scripts_from_others/GCContent/get_gc_content.pl

######## Filtering steps:
# 1. >10kb from exons
# 2. not inside CpG Island
# 3. not inside repeat region
# 4. normal GC content
# 5. doesn't blast to zebra fish
####### get exonic regions (once) and make sure is sorted: ############### 
### (do once) grep exon $gff | awk '{OFS="\t";print $1,$4-1,$5,$9}' | sort -k1,1 -k2,2n > $exonicRegions
# should be ~200,000 lines

############ HQ site coords ####################
# results of filtering snps (all populations; all nv and snps)
noCallFrac=1.0 # no filter at all. no limits on how many individuals must be called. but there is min dp 500 at site level
#prefix=all_8_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_${noCallFrac}_rmBadIndividuals_passingFilters
prefix=all_9_maxHetFilter_${maxHetFilter}_rmRelatives_rmAdmixed_passingBespoke_maxNoCallFrac_1.0_rmBadIndividuals_passingFilters

hqSites=$bedDir/${prefix}.sorted.merged.coords.bed # bed coords (sorted, merged) of sites from step 7 of filtering (comes from filtering_Step_2.sh)

############# Get distance of every set of sites from exonic regions in ferret genome ################
mkdir -p $wd/distanceFromExons # this dir will have info on distance of sites from exons
mkdir -p $wd/CpG_Islands
mkdir -p $wd/repeatRegions
mkdir -p $wd/get_fasta
mkdir -p $wd/GC_Content
mkdir -p $wd/zebra_fish
mkdir -p $wd/passing_sites


 # this dir will have neutral regions going through 3 checks: CpG Islands, GC content, and blast to fish
bedtools closest -d -a ${hqSites} -b ${exonicRegions} > $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt
#### NOTE: the output of this will be
# [Hq site info] [closest exon info] [distance between]; so I want the HQ sites that are >10000bp away from the closest exon.
# don't want it the other way around (getting info on each exon). want info on hq sites.

# last column (8) is the distance; want it to be at least 10,000, and want to keep
# track of the distance. Collect all that are >10,000 away. 
# pick the ones with high distance (awk) (get site totals below)
awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>100000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min100kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>1000000)print $1,$2,$3}' $wd/distanceFromExons/${prefix}.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/${prefix}.min1Mb.fromExon.0based.sorted.merged.bed

### Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)
######### get total amounts of sequence in each file: ########

> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
for i in `ls $wd/distanceFromExons/*bed`
do
echo $i >> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> $filterDir/filteringStats/totalSequenceByDistanceFromExons.txt
done


###### for now use this 10kb bed file for making the SFS
###################### Check neutral regions (10kb)###################

######### Check to see if any regions intersect with CpG Islands ##############
# got CpG islands from UCSC browser;
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/CpG_Islands)
# trying to intersect 
# get regions that DO NOT intersect with CpG islands
# *** -v **** this will output regions in "A" that ***DO NOT*** intersect with "B" (CpG Islands)
bedtools intersect -v -a $wd/distanceFromExons/${prefix}.min10kb.fromExon.0based.sorted.merged.bed -b $mfurCpG > $wd/CpG_Islands/${prefix}.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' 

########## 20181119 value:
#9,212,233 sites >10kb from genes (no max call frac) w no max no call frac filter
#9,062,287 after CpG filter (no max call frac filter)

########## 20180806 value
# total sequence before CpG filter: 9,030,072 w 0.9 max no call frac (7,818,882 when using 0.2 max no call frac)
# total sequence after CpG filter: 8,885,106 w 0.9 max no call frac  (6,816,729 when using 0.2 max no call frac)
# so lost ~100kb of sequence
# get CpG content
# overly conservative to filter >10kb from CpG Islands.
#bedtools closest -d -a $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed -b ${mfurCpG} > $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.distanceFromCpGIslands.0based.txt
#awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.distanceFromCpGIslands.0based.txt > $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.min10kb.fromCpG.0based.bed
 

######### Check to see if any regions intersect with RepeatMasker Regions ##############
# got repeat mask from UCSC browser
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/repeatMasker_UCSC)
# trying to intersect 
# get regions that DO NOT intersect with repeat regions

bedtools intersect -v -a $wd/CpG_Islands/${prefix}.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed -b $mfurRepeat > $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# check amount  of sequence lost:
# awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# 2018119 values: (no max call frac)
# 6,853,475 (so lost 2,208,812 sites (same as before))
# 20180806 values:
# 6,701,482 remain with 0.9 max no call frac (5,953,300 with 0.2 max no call frac) (lost 2Mb). Lost more with 0.9 setting at this stage: maybe because things that map poorly are in repeat regions more

##### just for the fasta file (not for the eventual list of HQ sites) re-merge the bed file with -d 10 setting (if things are 10bp apart you can still merge them)
# use this to get fasta (but not for SFS)
# doing this because if there is a single isolated base,
# it becomes its own entry in the fasta, which is going to make BLASTING a pain
# for now looking at overall region (allowing gaps of up to 10bp), not just called sites.
bedtools merge -d 10 -i $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed > $wd/get_fasta/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.mergedMaxDistance10.forFasta.notForSFS.bed
# * you are here *

###### get fasta sequence
bedtools getfasta -fi $REFERENCE -bed $wd/get_fasta/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.mergedMaxDistance10.forFasta.notForSFS.bed -fo $wd/get_fasta/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta

############# Get GC content of each part of Fasta (exclude if >50%?) ##############
### for now I'm not filtering on this; just generating it for interest.
perl $getGC $wd/get_fasta/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta > $wd/GC_Content/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.GC_Content.txt
# found that most regions were <=~70% GC, not going to filter further since I already got rid of CpG islands.
############# Blast against zebra fish genome to look for conservation #############
# retrieved 20180620
# do this once:
# cd $drerDir
# wget ftp://ftp.ensembl.org/pub/release-86/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
# gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
# makeblastdb -in Danio_rerio.GRCz10.dna.toplevel.fa -out Drer_blastdb -dbtype nucl
blastn -query $wd/get_fasta/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta -db  $drerDir/Drer_blastdb -outfmt 7 > $wd/zebra_fish/neutralBlast_ZebraFish_blastn.out
# based on output, get regions with e-value < 1e-10 to exclude. You are getting their coordinates from their fasta name, not from teh blast output
# so it is still 0-based even though blast output is 1-based.
# only lose 11kb of sequence
grep -v "#"  $wd/zebra_fish/neutralBlast_ZebraFish_blastn.out | awk '{if($11<1e-10)print $1}' | awk -F"[:-]" '{OFS="\t"; print $1,$2,$3}' | sort | uniq > $wd/zebra_fish/fish.matches.eval.1e-10.0based.bed
# then want to exclude those
bedtools intersect -v -a $wd/repeatRegions/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed -b $wd/zebra_fish/fish.matches.eval.1e-10.0based.bed > $wd/zebra_fish/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.bed
# awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $wd/zebra_fish/${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.bed
# new amount of sequence: 6,687,614 with 0.9 no call frax; (5,942,506 w 0.2 no call frac) # lost ~13kb of sequence from 6,701,48 to 6,687,614
# can then use the final bed file to make the SFS using Tanya's script.

# you can choose which of the sets of filters you want and update this accordingly. For now (20180820), it is the file that:
# is 10kb from exons
# is not in CpG island
# is not in repeat region
# does not blast to zebra fish
# no other GC content filter
finalBedDir=$wd/zebra_fish
finalBed=${prefix}.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.bed
cp $finalBedDir/$finalBed $wd/passing_sites/${finalBed%.bed}.useThis.bed
# also copy it to bedCoords
cp $finalBedDir/$finalBed $filterDir/bedCoords/${finalBed%.bed}.useThis.bed

# get final amount of sequence:
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $wd/passing_sites/${finalBed%.bed}.useThis.bed > $wd/passing_sites/totalPassingSequence.txt
# and copy this info to filteringStats as well
cp $wd/passing_sites/totalPassingSequence.txt $SCRATCH/captures/vcf_filtering/${rundate}_filtered/filteringStats/totalPassingNeutralSequence.txt

# 20181119 value: 6839537 ~6.8 with no max call frac
# 20180806 valu: 6687614 ~6.7 Mb with 0.9 max no call frac. (up from 5942506 under this is with 0.2 max no call frac) so I gained 745kb of sequence. Not very much.
# That is okay, shows it doesn't matter much which way you do it. So can keep neutral SFSs as-is, or remake based on these new coords
# going to remake to be thorough.
 
