pop=AK
model=AK.1D.2Epoch.35Gen.250Inds
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/neutralSimulations/populationSpecificModels/${model}
#todaysdate=`date +%Y%m%d` # don't want to use todays date because then different time starting arrays could get messed up
todaysdate=20200129
# 100 replicates; starting with 10
for i in {1..11}
do
# name job slimRepX
qsub -N slimRep${pop}_${i} $scriptdir/slim.${model}.array.sh $i $todaysdate
done

# this will submit 100 jobs, each of which is an array of 60 portions.

