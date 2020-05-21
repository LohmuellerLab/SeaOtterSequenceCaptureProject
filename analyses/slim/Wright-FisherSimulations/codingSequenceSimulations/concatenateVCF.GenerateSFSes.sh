########### this script will gather up simulation replicates within a single rundate (run again for different rundates)
# for the populations and models you specify
# you can also modify for differnent run dates by supplying the popsModelsRundates list  manually 
# it will gather up summaries of the simulations (one line per mutation) and will concanenate VCF files and make SFSes from them
##### SOMETHING TO LOOK INTO: Jacqueline found that the VCF files from SLIM were kind of weird, which is why she does summaries
# Make sure they are concordant!

# gather up replicates:
numReps=25 # total number of reps you ran 
DesiredReps=20 # how many you'll take out of the 25 (some randomly fail, so picking 20 / 25 for all)
numChunks=20 #
numStates=4 #
checkNumber=$((numChunks*numStates))

# process output of slim 
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/SLIM
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
# com is 3epoch (differnt model) 

rundate=20200310 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs="0.5 0" # set of hs you're interested in
popMods="CA/1D.3Epoch.LongerRecovery" # population and corresponding models you're interested in

# you can have multiple models per population just list as: AK/1D.2Epoch.1.5Mb.cds AK/OtherModel in the popMods variable
popsModelsRundates=""
for i in $popMods
do
for h in $hs
do
pm=$i/$rundate/h_$h
popsModelsRundates=`echo $popsModelsRundates $pm`
done
done # this sets up your list of pops, models, rundates and Hs in a list that looks like:
# AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/

# or you can set it up manually if you are using an odd mixture of models/dates --> 
# popsModelsRundates='AK/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ AL/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ CA/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/ KUR/1D.2Epoch.1.5Mb.cds/20190423/h_0.5/' # maybe? -- this is kind of awkward, maybe have to deal with diff populations differently?
# not ready yet: COM/1D.3Epoch.1.5Mb.cds/20190423/h_0.5/
# loop through models, populations and 25 replicates
scriptdir=$gitdir/slim_scripts
wd=$SCRATCH/Paper_data/SLIM 

################### first get lists of which runs worked #################

for popsModelsRundate in $popsModelsRundates
do
echo $popsModelsRundate
#pull out the population name for when you make sfses (note that if you change how you label popmodelsrundates you'll have to change how to do this)
pop=${popsModelsRundate%/*/*/*} #

############### get top passing runs ############################
# use a little script I made to detect the ones that passed and make lists
# usage script.sh [popModel] [total reps run] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]
chmod +x $scriptdir/FigureOutWhichRunsPassed.vcf.sh # make sure it's executable
$scriptdir/FigureOutWhichRunsPassed.vcf.sh $popsModelsRundate $numReps $DesiredReps $numChunks $numStates # this will make files of the ones that passed 
passingFile=passingReps.FIRST.${DesiredReps}.usethis.txt 
# make sure the file exists and is $DesiredReps long 
if [ ! -f $wd/$popsModelsRundate/$passingFile ]
then
break
fi

##################################################################
# make the slim script from the maker script:
vcfOutDir=$wd/concattedVCFs/$popsModelsRundate/
sfsDir=$wd/SFSes/$popsModelsRundate/
sfsScriptDir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/scripts/sfs_scripts

mkdir -p $vcfOutDir
mkdir -p $sfsDir

states="PreContraction PostContraction PostRecovery EndSimulation"


for state in $states # different states sampled from 
do
cat $wd/$popsModelsRundate/$passingFile | while read i # instead of a for loop, going to read through my passing file. i is replicate_x not just a number 
do 
echo "$i ${state}"
repdir=$SCRATCH/slim/$popsModelsRundate/${i}


# vcf with all mutations together:
outVCF=$vcfOutDir/${i}.slim.output.${state}.allConcatted.vcf

# mutation subtypes vcfs:
mut1VCF=$vcfOutDir/${i}.slim.output.${state}allConcatted.mutType1.vcf
# mutation subtypes vcfs:
mut2VCF=$vcfOutDir/${i}.slim.output.${state}.allConcatted.mutType2.vcf

#vcf headers:
grep "#" $repdir/slim.output.${state}.1.vcf > $outVCF
grep "#" $repdir/slim.output.${state}.1.vcf > $mut1VCF
grep "#" $repdir/slim.output.${state}.1.vcf > $mut2VCF


# then loop over all chunks
for j in $(seq 1 $numChunks)
do

# concat vcfs and gzip
echo "concatenating vcf chunks"
# # want to select everything but first column (setting it to "") and replace with a chromosome 'identifier' that is the chunk nunber
# exclude header 
# make sure to be tab separated! otherwise python script won't know what to do
grep -v "#" $repdir/slim.output.${state}.${j}.vcf | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >>  $outVCF
# # separate mutation types
# mutation type 1: 
grep -v "#" $repdir/slim.output.${state}.${j}.vcf | grep "MT=1"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >> $mut1VCF
# mutation type 2:
grep -v "#" $repdir/slim.output.${state}.${j}.vcf | grep "MT=2"  | awk -v chr=$j '{OFS="\t";$1=""; print chr,$0}' >> $mut2VCF

done
# gzip the outputs:
echo "Gzipping results"
gzip -f $outVCF
gzip -f $mut1VCF
gzip -f $mut2VCF
# make SFSes: mut type 1
# to just get population name, not whole model date etc do
echo "making SFSes"
python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut1VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}.slim.mutType1

# make SFSes: mut type 2

python $sfsScriptDir/generate1DSFS.py \
--vcf ${mut2VCF}.gz \
--pop ${pop}_sim \
--outdir $sfsDir \
--outPREFIX $i.${state}.slim.mutType2
done
done
done
#done
