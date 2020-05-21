#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=6G
#$ -m abe
#$ -M pkalhori
#$ -t 1-20
# eventually do 1-20 (1.5Mb * 20 = 30Mb, ~ an exome)
# note: this array can be used for any slim model
# going to simulate 1.5MB X 14 (=21MB)
# each 1.5Mb 'chromosome' will have 1000 genes
# each gene is 1500bp long and has an internal 1e-08 recomb rate
# with 1e-03 recomb rate between genes
# genes can have NS or S mutations with ratio 2.31:1
# S mutations are neutral
# NS mutations are drawn from bernard gamma dist DFE, and can have any dom coefficient (set as a variable in this script)
# this script produces a VCF file which you can parse 
# not sure if I want fixed mutations to be included or not. I kind of think I do want them because I count them in my empirical data

######### run parameters -- change these across models ###########
# want a total of 6000 blocks. so 60 instances of this script for one replicate.
pop=$1
model=$2 #1D.2Epoch.1.5Mb.cds 
rep=$3 # doing one replicate, then will set from command line from submission script
rundate=$4 # date arrays are submitted; set in submitter so as not to have jobs on different days
h=$5
# submitter usage: qsub -N name -o outdir -e errordir $script $pop $model $rep $rundate

#wd=$SCRATCH/captures/analyses/slim/cdsSimulations/$pop/$model/h_$h/$rundate/
######### programs #########
# load the proper gcc 
source /u/local/Modules/default/init/modules.sh
module load gcc/6.3.0
slim=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/software/slim_build/slim # location of slim

############## files and dirs ############
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/SLIM/ # project github location
scriptdir=$gitdir/slim_scripts/$pop/$model # location of slim scripts
# this seed will account for runs run on differet days, and if they start during the same second there will also be adjustments for the run number/task id
# it is not perfectly random, but in general few runs will set the seed at the exact time and the exactness is pretty good (598928/600000 seeds were unique when tested with a for loop)
# in real life, the array jobs will be starting at lots of different times, so mostly it isn't an issue. just a small percentage might start at exact same time
todaysdate=`date +%Y%m%d`
seed=$(($todaysdate+$RANDOM+(($RANDOM*$rep*10))+$SGE_TASK_ID)) # uses date, plus random , plus the replicate and SGE task id. so no task should be the same, and no jobs run on different days should be the same, even if they have same task id
# so if two tasks with the same task id across different reps start in the same second (get same random), then they will still be different because of 10*rep 

#for h in 0 0.5
#do
wd=$SCRATCH/slim/$pop/$model/$rundate/h_$h
outdir=$wd/replicate_${rep} # set this in submission script 
mkdir -p $wd
mkdir -p $outdir
slimscript=slim_elut_${model}_${pop}_h${h}.job # specific slim script 
cp $scriptdir/$slimscript $wd/$slimscript.AsRunOn.$todaysdate # make a record of the script as it was run; this is inefficient, copies it for each task in the array
######## parameters #############
$slim \
-long \
-seed $seed \
-d outdir="'$outdir'" \
-d v_CHUNK=${SGE_TASK_ID} \
-d v_REP=${rep} \
$scriptdir/$slimscript
# script comes at end
# long prints out what vars are set as 
# this will output vcf and full population state.

#done
sleep 2m


