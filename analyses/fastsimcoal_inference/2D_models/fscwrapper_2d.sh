\#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=10G
#$ -N Fsc2D_Models
#$ -o /u/scratch/p/pkalhori/fastsimcoal/reports
#$ -e /u/scratch/p/pkalhori/fastsimcoal/reports
#$ -m abe
#$ -M pkalhori
#$ -t 1:50

deme0=AK
deme1=CA
pops="neutral.${deme1}.${deme0}"
models="2D.2Epoch 2D.2Epoch.Mig.Symmetric 2D.3Epoch.FixedContraction.Time.Sizes 2D.3Epoch.Mig.Symmetric.FixedContraction.Time.Sizes"
rundate=`date +%Y%m%d`
#this is the string of populations to loop through
# wd stands for "working directory"
wd=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/fastsimcoal
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/fastsimcoal
fsc=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/software/fsc26_linux64/fsc26
for model in $models
#iterate through each model
do
for pop in $pops
#iterate through each population
do
cd $wd
header=${model}_${pop}_$rundate
mkdir $header
cd $header
cp $gitdir/new_neutral.CA.AK/${model}_${pop}/${model}_${pop}.tpl ./
cp $gitdir/new_neutral.CA.AK/${model}_${pop}/${model}_${pop}.est ./
cp $gitdir/SFSes/${pop}_jointMAFpop1_0.obs ${model}_${pop}_jointMAFpop1_0.obs
mkdir $wd/$header/run_${SGE_TASK_ID}
cd $wd/$header/run_${SGE_TASK_ID}
cp $wd/$header/${model}_${pop}.tpl $wd/$header/${model}_${pop}.est $wd/$header/${model}_${pop}_jointMAFpop1_0.obs ./
#ss0=`grep -w $deme0 $wd/projectionValues.txt|awk '{print$2}'`
#sed -i "s/SAMPLE_SIZE_0/$ss0/g" ${model}_${pop}.tpl
#ss1=`grep -w $deme1 $wd/projectionValues.txt|awk '{print$2}'`
#sed -i "s/SAMPLE_SIZE_1/$ss1/g" ${model}_${pop}.tpl
$fsc -t ${model}_${pop}.tpl -n100000 -m -e ${model}_${pop}.est -M -L 50 
done
done
cd $wd
sleep 10m 


