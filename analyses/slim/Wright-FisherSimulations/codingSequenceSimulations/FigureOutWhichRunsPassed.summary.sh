### Come up with a way to detect / deal with failures

# Did 25 replicates
# pick the first 20 complete replicates, leave the rest behind
# to be complete need to have 2x the number of chunks

# usage script.sh [popModel]  [number of reps carried out total] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]


popModel=$1 # in format pop/model/h_$h
numReps=$2 # total number of reps run 
DesiredReps=$3 # how many you'll take
numChunks=$4 #
numStates=$5 #
checkNumber=$((numChunks*numStates))


wd=/u/scratch/p/pkalhori/Paper_Data/SLIM/raw_data/$popModel/

> $wd/passingReps.summaries.txt
> $wd/failingReps.summaries.txt
for i in $(seq 1 $numReps)
do
summaryCount=0 # restart it each time so that it's never null
summaryCount=`ls $wd/replicate_${i}/*summary.txt | wc -l`
echo $summaryCount
# eq is equal
if [ $summaryCount -eq $checkNumber ]
then
echo replicate_$i >> $wd/passingReps.summaries.txt
passing=$((passing+1))
# ne is not equal
elif [ $summaryCount -ne $checkNumber ]
then
echo $replicate_$i >> $wd/failingReps.summaries.txt
fi
done

# then get the top DesiredReps passing reps 
head -n$DesiredReps $wd/passingReps.summaries.txt > $wd/passingReps.FIRST.$DesiredReps.usethis.txt

# make sure 20 made it through:

lineCount=`wc -l $wd/passingReps.FIRST.$DesiredReps.usethis.txt | awk '{print $1}'`
if [ $lineCount -ne $DesiredReps ]
then
echo "$wd: FEWER THAN $DesiredReps FINISHED --- REDO SIMULATION"
mv $wd/passingReps.FIRST.$DesiredReps.usethis.txt $wd/NOT.ENOUGH.PASSED.txt
fi

