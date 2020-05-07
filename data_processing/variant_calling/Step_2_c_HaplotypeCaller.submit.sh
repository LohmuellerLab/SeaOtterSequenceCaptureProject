#! /bin/bash
#$ -cwd
#$ -l h_rt=3:00:00,h_data=3G,highp
#$ -o /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -e /u/flashscratch/a/ab08028/captures/reports/submissions/
#$ -m bea
#$ -M ab08028
#$ -N hapCaller_submit
# submit step 2c
user=ab08028 # where emails are sent

QSUB=/u/systems/UGE8.0.1vm/bin/lx-amd64/qsub

# location of github:  which may be on remote server or local laptop
gitDir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/
# scripts:
scriptDir=$gitDir/scripts/scripts_20180521/data_processing/variant_calling
# script to run: 
scriptname=Step_2_c_HaplotypeCaller.sh # change this to final script name!! 
# 
# file locations:
SCRATCH=/u/flashscratch/a/ab08028
wd=$SCRATCH/captures
headers=$wd/samples/capture2.txt # all Elut samples, modern, ancient and blank
reports=/u/flashscratch/a/ab08028/captures/reports/GATK/
mkdir -p $reports

cat $headers | while read header
do
errorLocation=${reports}/${header} # report location
mkdir -p $errorLocation
$QSUB -e $errorLocation -o $errorLocation -M $user -N hapCaller.${header} $scriptDir/$scriptname $header
sleep 1m
done
