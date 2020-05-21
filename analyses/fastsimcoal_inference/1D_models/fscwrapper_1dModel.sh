#! /bin/bash
#$ -cwd
#$ -l h_rt=20:00:00,h_data=28G
#$ -N {name of job}
#$ -o {path to output files}
#$ -e {path to error files}
#$ -m abe
#$ -M {user}
#$ -t 1-50:1

pops='AK AL CA COM KUR'
models='1D.1Epoch 1D.2Epoch 1D.3Epoch'
rundate=`date +%Y%m%d`
#this is the string of populations to loop through
# wd stands for "working directory"
wd={set working directory}
md={path to .est and .tpl files}
fsc={path to fsc26 program}
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
cp $md/$model.tpl ${model}_${pop}.tpl
cp $md/$model.est ${model}_${pop}.est
cp $wd/SFSdir/${pop}_MAFpop0.obs ${model}_${pop}_MAFpop0.obs
mkdir $wd/$header/run_${SGE_TASK_ID}
cd $wd/$header/run_${SGE_TASK_ID}
cp $wd/$header/${model}_${pop}.tpl $wd/$header/${model}_${pop}.est $wd/$header/${model}_${pop}_MAFpop0.obs ./
ss=`grep -w $pop $wd/projectionValues.txt|awk '{print$2}'`
sed -i "s/SAMPLE_SIZE/$ss/g" ${model}_${pop}.tpl
$fsc -t ${model}_${pop}.tpl -n100000 -m -e ${model}_${pop}.est -M -L 50 -q
done
done
cd $wd
sleep 10m 


