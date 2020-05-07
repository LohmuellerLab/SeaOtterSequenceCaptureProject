################################# Set up ###########################
# modules
source /u/local/Modules/default/init/modules.sh
module load bedtools
module load blast

rundate=20180806

SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures/vcf_filtering/${rundate}_filtered/checkingNeutralSites
mkdir -p $wd/
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
hqSites=$wd/bedCoords/all_7_passingBespoke.sorted.merged.coords.bed # bed coords (sorted, merged) of sites from step 7 of filtering (comes from filtering_Step_2.sh)

############# Get distance of every set of sites from exonic regions in ferret genome ################
mkdir -p $wd/distanceFromExons # this dir will have info on distance of sites from exons
mkdir -p $wd/CpG_Islands
mkdir -p $wd/repeatRegions
mkdir -p $wd/get_fasta
mkdir -p $wd/GC_Content
mkdir -p $wd/zebra_fish
mkdir -p $wd/passing_sites


# this dir will have neutral regions going through 3 checks: CpG Islands, GC content, and blast to fish
bedtools closest -d -a ${hqSites} -b ${exonicRegions} > $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt
#### NOTE: the output of this will be
# [Hq site info] [closest exon info] [distance between]; so I want the HQ sites that are >10000bp away from the closest exon.
# don't want it the other way around (getting info on each exon). want info on hq sites.

# last column (8) is the distance; want it to be at least 10,000, and want to keep
# track of the distance. Collect all that are >10,000 away. 
# pick the ones with high distance (awk)
awk -F'\t' '{OFS="\t";if($8>10000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>100000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min100kb.fromExon.0based.sorted.merged.bed

awk -F'\t' '{OFS="\t";if($8>1000000)print $1,$2,$3}' $wd/distanceFromExons/all_7_passingBespoke.distanceFromExons.0based.txt |  sort -k1,1 -k2,2n | bedtools merge -i stdin > $wd/distanceFromExons/all_7_passingBespoke.min1Mb.fromExon.0based.sorted.merged.bed

### Note: 1,2,3 columns are the HQ SITES position, NOT the position of the exon. (If you mess up what is a and b in bedtools closest this would be messed up)
######### get total amounts of sequence in each file: ########

> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
for i in `ls $wd/distanceFromExons/*bed`
do
echo $i >> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $i >> $wd/filteringStats/totalSequenceByDistanceFromExons.txt
done


###### for now use this 10kb bed file for making the SFS
###################### Check neutral regions (10kb)###################

######### Check to see if any regions intersect with CpG Islands ##############
# got CpG islands from UCSC browser;
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/CpG_Islands)
# trying to intersect 
# get regions that DO NOT intersect with CpG islands
# *** -v **** this will output regions in "A" that ***DO NOT*** intersect with "B" (CpG Islands)
bedtools intersect -v -a $wd/distanceFromExons/all_7_passingBespoke.min10kb.fromExon.0based.sorted.merged.bed -b $mfurCpG > $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed
# total sequence before CpG filter: 7,818,882
# total sequence after CpG filter: 6,816,729
# so lost ~100kb of sequence


######### Check to see if any regions intersect with RepeatMasker Regions ##############
# got repeat mask from UCSC browser;
# added .1 to each scaffold name so it matches my reference (see script in mfurDir/repeatMasker_UCSC)
# trying to intersect 
# get regions that DO NOT intersect with repeat regions
bedtools intersect -v -a $wd/CpG_Islands/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.0based.sorted.merged.bed -b $mfurRepeat > $wd/repeat_regions/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# check amount  of sequence lost:
# awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'  $wd/repeat_regions/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed
# 5,953,300 remaining (lost 800kb)

##### re-merge the bed file with -d 10 setting (if things are 10bp apart you can still merge them)
# use this to get fasta
# doing this because if there is a single isolated base,
# it becomes its own entry in the fasta, which is going to make BLASTING a pain
# for now looking at overall region (allowing gaps of up to 10bp), not just called sites.
bedtools merge -d 10 -i $wd/repeat_regions/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed > $wd/get_fasta/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.mergedMaxDistance10.forFasta.notForSFS.bed

###### get fasta sequence
bedtools getfasta -fi $REFERENCE -bed $wd/get_fasta/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.mergedMaxDistance10.forFasta.notForSFS.bed -fo $wd/get_fasta/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta

############# Get GC content of each part of Fasta (exclude if >50%?) ##############
### for now I'm not filtering on this; just generating it for interest.
perl $getGC $wd/get_fasta/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta > $wd/GC_Content/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.GC_Content.txt
# found that most regions were <=~70% GC, not going to filter further since I already got rid of CpG islands.
############# Blast against zebra fish genome to look for conservation #############
# retrieved 20180620
# do this once:
# cd $drerDir
# wget ftp://ftp.ensembl.org/pub/release-86/fasta/danio_rerio/dna/Danio_rerio.GRCz10.dna.toplevel.fa.gz
# gunzip Danio_rerio.GRCz10.dna.toplevel.fa.gz
# makeblastdb -in Danio_rerio.GRCz10.dna.toplevel.fa -out Drer_blastdb -dbtype nucl
blastn -query $wd/get_fasta/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.fasta -db  $drerDir/Drer_blastdb -outfmt 7 > $wd/zebra_fish/neutralBlast_ZebraFish_blastn.out
# based on output, get regions with e-value < 1e-10 to exclude. You are getting their coordinates from their fasta name, not from teh blast output
# so it is still 0-based even though blast output is 1-based.
# only lose 11kb of sequence
grep -v "#"  $wd/zebra_fish/neutralBlast_ZebraFish_blastn.out | awk '{if($11<1e-10)print $1}' | awk -F"[:-]" '{OFS="\t"; print $1,$2,$3}' | sort | uniq > $wd/zebra_fish/fish.matches.eval.1e-10.0based.bed
# then want to exclude those
bedtools intersect -v -a $wd/repeat_regions/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.0based.sorted.merged.bed -b $wd/zebra_fish/fish.matches.eval.1e-10.0based.bed > $wd/zebra_fish/all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.bed
# new amount of sequence: 5,942,506 # lost ~10kb of sequence
# want to find a way to exclude any regions that match (there shouldn't be any, since I already did this with elut sequences)
# can then use the final bed file to make the SFS using Tanya's script.
finalBedDir=$wd/zebra_fish
finalBed=all_7_passingBespoke.min10kb.fromExon.noCpGIsland.noRepeat.noFish.0based.sorted.merged.bed
cp $finalBedDir/$finalBed $wd/passing_sites/${finalBed%.bed}.useThis.bed
# also copy it to bedCoords
cp $finalBedDir/$finalBed ${wd%checkingNeutralSites}/bedCoords/${finalBed%.bed}.useThis.bed
# you can choose which of the sets of filters you want and update this accordingly. For now (20180820), it is the file that:
# is 10kb from exons
# is not in CpG island
# is not in repeat region
# does not blast to zebra fish
# no other GC content filter

# get final amount of sequence:
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' $wd/passing_sites/${finalBed%.bed}.useThis.bed > $wd/passing_sites/totalPassingSequence.txt
# 5942506 ~ 6Mb