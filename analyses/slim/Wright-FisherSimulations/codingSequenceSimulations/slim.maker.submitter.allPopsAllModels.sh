#gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/SLIM
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/SLIM
#models='2D.3Epoch.NoTranslocation 2D.3Epoch.Translocation.1perGen 2D.3Epoch.Translocation.5perGen 2D.3Epoch.Translocation.10perGen 2D.3Epoch.Translocation.25perGen 2D.3Epoch.Translocation.25for2Gen'
models=1D.5Epoch
#populations='AK AL CA COM KUR'
populations=AK # do COM separately below
# loop through models, populations and 25 replicates
scriptdir=$gitdir/slim_scripts

todaysdate=`date +%Y%m%d`
for pop in $populations
do
for model in $models
do
# make the slim script from the maker script:
wd=$SCRATCH/slim/$pop/$model/$todaysdate/
mkdir -p $wd
# make a dir to put logs in:
logdir=$wd/logs
mkdir -p $logdir
# make the slim script:
for h in 0 0.5
do
sh $scriptdir/$pop/$model/make_slim_elut.${model}.${pop}.sh $h
#done

for i in {1..25}
do
# qsub -N name -o outdir -e errordir $script $pop $model $rep $rundate
qsub -N slimRep${i}.${pop}.${model}.${h} -o $logdir -e $logdir $scriptdir/array_slim_elut.generic.sh $pop $model $i $todaysdate $h
done
done
done
done

# need to do COM 3 Epoch separately #
#for pop in COM
#do
#for model in '1D.3Epoch.1.5Mb.cds'
#do
#for h in 0.5 # 0
#do
# make the slim script from the maker script:
#wd=$SCRATCH/captures/analyses/slim/cdsSimulations/$pop/$model/$todaysdate/
#mkdir -p $wd
# make a dir to put logs in:
#logdir=$wd/logs
#mkdir -p $logdir
# make the slim script:
#sh $scriptdir/$pop/$model/make_slim_elut.${model}.${pop}.sh $h
#done
#for i in {1..1}
#do
# qsub -N name -o outdir -e errordir $script $pop $model $rep $rundate
#qsub -N slimRep${i}.${pop} -o $logdir -e $logdir $scriptdir/array_slim_elut.generic.sh $pop $model $i $todaysdate
#done
#done
#done
