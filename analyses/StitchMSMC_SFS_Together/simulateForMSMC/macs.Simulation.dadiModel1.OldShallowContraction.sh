#!/bin/bash
#$ -cwd
#$ -l h_rt=5:00:00,h_data=28G,highp
#$ -N simDadiModel1
#$ -m abe
#$ -M ab08028
#$ -o /u/flashscratch/a/ab08028/captures/reports/MaCS
#$ -e /u/flashscratch/a/ab08028/captures/reports/MaCS
#$ -t 1-10
source /u/local/Modules/default/init/modules.sh
module load anaconda
# conda create -n MaCsSimulations python=3.6 # only once
source activate MaCsSimulations
wd=/u/flashscratch/a/ab08028/captures/analyses/simulateForMSMC
cd $wd
macsFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/macs
msformatterFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msformatter
ms2multiFile=/u/home/a/ab08028/klohmueldata/annabel_data/bin/msmc-tools/ms2multihetsep.py
rundate=`date +%Y%m%d`
replicate=$SGE_TASK_ID
model=dadiModel1.OldShallowContraction
mkdir -p ${model}
for j in {1..6}
do
outdir=$wd/${model}/rep_${replicate}/group_$j.${model}
mkdir -p $outdir
cd $outdir
cp -n $macsFile $outdir
cp -n $msformatterFile $outdir
cp -n $ms2multiFile $outdir
for i in {1..10}
do
# dadi model 1 for msmc
mu=8.64e-09
r=1e-08
Na=4500.0
rho=0.00018
theta=0.00015552
date=`date +%Y%m%d`
SEED=$((date+$RANDOM+((j-1)*10)+i))
# this is a new addition! need to have a different random seed for each simulation; if they start within a second of each other, they will have the same seed. not an issue for big simulations of 30Mb because those are slow, but 100kb can start within a second of each other!
./macs 2 30000000 -t 0.00015552 -r 0.00018 -s $SEED -eN 0.0 0.4 -eN 0.06 1  > $outdir/group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt
#convert to ms format
./msformatter < $outdir/group_${j}_block_${i}.${model}.macsFormat.OutputFile.${rundate}.txt > $outdir/group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt
#convert to msmc input format
python ./ms2multihetsep.py $i 30000000 < $outdir/group_${j}_block_${i}.${model}.msFormat.OutputFile.${rundate}.txt > $outdir/group_${j}_block_${i}.${model}.MSMCFormat.OutputFile.${rundate}.txt
done
cd $wd
done
sleep 5m
source deactivate # deactivate conda env
