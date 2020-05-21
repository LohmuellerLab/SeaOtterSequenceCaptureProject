#! /bin/bash
#$ -cwd
#$ -l h_rt=05:00:00,h_data=2G
#$ -o /u/flashscratch/a/ab08028/captures/reports/slim
#$ -e /u/flashscratch/a/ab08028/captures/reports/slim
#$ -m abe
#$ -M ab08028
#$ -t 1-60
####### this is the simplified model I used for all my AK cds simulations ######

######### 2 Epoch script generates 100 x 1kb independent blocks ###########
# want a total of 6000 blocks. so 60 instances of this script for one replicate.
pop=AK
model=AK.1D.2Epoch.35Gen.250Inds
rep=$1 # doing one replicate, then will set from command line from submission script
rundate=$2 # date arrays are submitted; set in submitter so as not to have jobs on different days
outdir=$SCRATCH/captures/analyses/slim/neutralSimulations/$model/$rundate/replicate_${rep} # set this in submission script 
mkdir -p $outdir
######### programs #########
# load the proper gcc 
source /u/local/Modules/default/init/modules.sh
module load gcc/6.3.0
slim=/u/home/a/ab08028/klohmueldata/annabel_data/bin/SLiM/slim_build/slim # location of slim

############## files and dirs ############
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # project github location
scriptdir=$gitdir/scripts/scripts_20180521/analyses/slim/neutralSimulations/populationSpecificModels/$model # location of slim scripts
slimscript=AK.1D.2Epoch.35Gen.250Inds.slim # specific slim script


######## parameters #############
todaysdate=`date +%Y%m%d`

# this seed will account for runs run on differet days, and if they start during the same second there will also be adjustments for the run number/task id
# it is not perfectly random, but in general few runs will set the seed at the exact time and the exactness is pretty good (598928/600000 seeds were unique when tested with a for loop)
# in real life, the array jobs will be starting at lots of different times, so mostly it isn't an issue. just a small percentage might start at exact same time
seed=$(($todaysdate+$RANDOM+(($RANDOM*$rep*10))+$SGE_TASK_ID)) # uses date, plus random , plus the replicate and SGE task id. so no task should be the same, and no jobs run on different days should be the same, even if they have same task id
# so if two tasks with the same task id across different reps start in the same second (get same random), then they will still be different because of 10*rep 
mu=8.64e-9
r=1e-8
ss=6 # sample size in individuals
nanc=4500 # ancestral size from dadi for CA, simplified
nu=250 #
# keeping track of which sript was run to make sure no mixups are made --> 
echo "script run on $todaysdate: $slimscript" > $outdir/parameters.txt
echo "pop: $pop; ss: $ss ; nu=$nu ; nanc= $nanc ; mu: $mu ; r: $r ; seed: $seed" >> $outdir/parameters.txt

$slim \
-long \
-seed $seed \
-d v_MU=$mu \
-d v_R=$r \
-d v_SS=$ss \
-d v_NANC=$nanc \
-d v_NU=$nu \
-d v_OUTFILE="'$outdir'" \
-d chunk=${SGE_TASK_ID} \
$scriptdir/$slimscript
# script comes at end
# long prints out what vars are set as 
# this will output vcf and full population state.

sleep 10m

