#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=8G
#$ -N dadi_inference
#$ -o /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -e /u/flashscratch/a/ab08028/captures/reports/dadi
#$ -m abe
#$ -M ab08028


######################### set up virtual environment #####################
# if using on Hoffman, have to use virtual environment:
# you can set up the virtual env by doing:
# module load python/2.7.13_shared
# virtualenv $HOME/env_python2.7.13 # once
# then activate it every future time with  source $HOME/env_python2.7.13/bin/activate
# make sure you pip install numpy scipy matplotlib 
# and then you can deactivate

source /u/local/Modules/default/init/modules.sh
module load python/2.7.13_shared
source /u/home/a/ab08028/env_python2.7.13/bin/activate # activate virtual environment 

SCRATCH=/u/flashscratch/a/ab08028/
gitdir=/u/home/a/ab08028/klohmueldata/annabel_data/OtterExomeProject/ # hoffman
scriptdir=$gitdir/scripts/scripts_20180521/analyses/dadi_inference
#script=$1 # generic -- give it any of your dadi scripts -- make a version of this for Hoffman too.
#model=${script%.dadi.py}
mu=8.64411385098638e-09
genotypeDate=20181119 # newer gts
sfsDate=20190301 # projection with 0.75 het filter and these projection values
hetFilter=0.75
todaysdate=`date +%Y%m%d`-MSMCMODEL
captures=$SCRATCH/captures/
sfsdir=$captures/analyses/SFS/$genotypeDate/easySFS/neutral/projection-${sfsDate}-hetFilter-${hetFilter}/dadi-plusMonomorphic/
dadidir=$captures/analyses/dadi_inference/
sfssuffix=plusMonomorphic.sfs
### Make sure this is the correct file #####

scripts="1D.CA.PSMC.Trim22.dadi.py 1D.CA.PSMC.Trim22.Simplified.dadi.py" # trim22 models for inference 
for pop in CA # only AL for now
do
# get total sites from total sites file that was written out as part of my easySFS scripts
L=`grep $pop $sfsdir/$pop-[0-9]*.totalSiteCount.L.withMonomorphic.txt | awk '{print $2}'`

for script in $scripts
do
model=${script%.dadi.py}
echo "starting inference for $pop for model $model"
outdir=$dadidir/$genotypeDate/$pop/inference_$todaysdate/$model/
mkdir -p $outdir
# carry out inference with 50 replicates that start with different p0 perturbed params:
for i in {1..50}
do
echo "carrying out inference $i for model $model for pop $pop" 
# [0-9] indicates that it's a number, but not specific about proj value
python $scriptdir/$script --runNum $i --pop $pop --mu $mu --L $L --sfs ${sfsdir}/${pop}-[0-9]*.${sfssuffix} --outdir $outdir
done


echo "concatenating results"
grep rundate -m1 $outdir/${pop}.dadi.inference.${model}.runNum.1.output > $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
for i in {1..50}
do
grep rundate -A1 $outdir/${pop}.dadi.inference.${model}.runNum.${i}.output | tail -n1 >> $outdir/${pop}.dadi.inference.${model}.all.output.concatted.txt
done

done
done
############################ deactivate virtualenv ###############
deactivate # deactivate virtualenv
