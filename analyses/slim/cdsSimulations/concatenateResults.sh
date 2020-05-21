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
numStates=1 #
checkNumber=$((numChunks*numStates))

# process output of slim 
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/final_paper_scripts
# choose the specific combination of populations/models/rundates you want? this is awkward... what is best way to do it?
# com is 3epoch (differnt model) 

rundate=20200311 # set of simulations you're interested in (if is across different rundates you can list popsModelsRundates explicitly)
hs="0 0.5" # set of hs you're interested in
#popMods="CA_AK/2D.3Epoch.NoTranslocation CA_AK/2D.3Epoch.Translocation.25for2Gen CA_AK/2D.3Epoch.Translocation.1perGen CA_AK/2D.3Epoch.Translocation.5perGen CA_AK/2D.3Epoch.Translocation.10perGen CA_AK/2D.3Epoch.Translocation.25perGen" # population and corresponding models you're interested in
popMods=AK/1D.5Epoch
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
scriptdir=$gitdir/SLIM_final
wd=$SCRATCH/Paper_Data/SLIM/raw_data/ 

################### first get lists of which runs worked #################

for popsModelsRundate in $popsModelsRundates
do
echo $popsModelsRundate
#pull out the population name for when you make sfses (note that if you change how you label popmodelsrundates you'll have to change how to do this)
pop=${popsModelsRundate%/*/*/*} #

############### get top passing runs ############################
# use a little script I made to detect the ones that passed and make lists
# usage script.sh [popModel] [total reps run] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]
chmod +x $scriptdir/FigureOutWhichRunsPassed.summary.sh # make sure it's executable
$scriptdir/FigureOutWhichRunsPassed.summary.sh $popsModelsRundate $numReps $DesiredReps $numChunks $numStates # this will make files of the ones that passed 
passingFile=passingReps.FIRST.${DesiredReps}.usethis.txt 
# make sure the file exists and is $DesiredReps long 
if [ ! -f $wd/$popsModelsRundate/$passingFile ]
then
break
fi

##################################################################
# make the slim script from the maker script:
summaryOutDir=$wd/concattedSummaries/$popsModelsRundate/

mkdir -p $summaryOutDir


cat $wd/$popsModelsRundate/$passingFile | while read i # instead of a for loop, going to read through my passing file. i is replicate_x not just a number 
do 
echo "$i ${state}"
repdir=$SCRATCH/Paper_Data/SLIM/raw_data/$popsModelsRundate/${i}

outSummary=$summaryOutDir/${i}.slim.output.allConcatted.summary.txt

# summary header:

grep "replicate" $repdir/slim.output.1.summary.txt > $outSummary

# then loop over all chunks
for j in $(seq 1 $numChunks)
do
# concat summaries and gzip
echo "concatenating summaries"
#  grep -v "^$" removes extra blank lines

grep -v "replicate" $repdir/slim.output.${j}.summary.txt | grep -v "^$" >> $outSummary

done
# gzip the outputs:
echo "Gzipping results"
gzip -f $outSummary

done
done

