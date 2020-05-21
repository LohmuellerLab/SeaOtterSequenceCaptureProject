### Come up with a way to detect / deal with failures

# Did 11 replicates
# pick the first 10 complete replicates, leave the rest behind
# to be complete need to have 2x the number of chunks

# usage script.sh [popModel]  [number of reps carried out total] [Desired replicates to pull out] [numChunks per replicate to check for completeness] [number of states per replicate (e.g. pre and post contraction)]

# put this in more context 
model=CA.1D.2Epoch.35Gen.200Inds/20200224 # in format pop/model/h_$h

numReps=11 # total number of reps run 
DesiredReps=10 # how many you'll take
numChunks=60 #
numStates=2 #
checkNumber=$((numChunks*numStates))


wd=/u/flashscratch/a/ab08028/captures/analyses/slim/neutralSimulations/CA.1D.2Epoch.35Gen.200Inds/20200224/
> $wd/passingReps.txt
> $wd/failingReps.txt
for i in $(seq 1 $numReps)
do
vcfCount=0 # restart it each time so that it's never null
vcfCount=`ls $wd/replicate_${i}/*vcf | wc -l`
echo $vcfCount
# eq is equal
if [ $vcfCount -eq $checkNumber ]
then
echo replicate_$i >> $wd/passingReps.txt
passing=$((passing+1))
# ne is not equal
elif [ $vcfCount -ne $checkNumber ]
then
echo $replicate_$i >> $wd/failingReps.txt
fi
done

# then get the top DesiredReps passing reps 
head -n$DesiredReps $wd/passingReps.txt > $wd/passingReps.FIRST.$DesiredReps.usethis.txt

# make sure 10 made it through:

lineCount=`wc -l $wd/passingReps.FIRST.$DesiredReps.usethis.txt | awk '{print $1}'`
if [ $lineCount -ne $DesiredReps ]
then
echo "$wd: FEWER THAN $DesiredReps FINISHED --- REDO SIMULATION"
mv $wd/passingReps.FIRST.$DesiredReps.usethis.txt $wd/NOT.ENOUGH.PASSED.txt
fi

