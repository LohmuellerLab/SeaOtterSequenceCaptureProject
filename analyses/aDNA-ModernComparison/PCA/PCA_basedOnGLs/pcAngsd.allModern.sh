#! /bin/bash
#$ -cwd
#$ -l h_rt=24:00:00,h_data=2G
#$ -m abe
#$ -pe shared 10
#$ -M ab08028
#$ -N PCAngsd
#$ -e /u/flashscratch/a/ab08028/captures/reports/angsd
#$ -o /u/flashscratch/a/ab08028/captures/reports/angsd

##### pcangsd -- run PCAngsd to do PCA based on genotype likelihoods
# set a MAF of 0.05
# maybe do some sort of SNP likelihood filtering? 1e-6?

source /u/local/Modules/default/init/modules.sh
module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
pcangsddir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pcangsd
wd=$SCRATCH/captures/aDNA-ModernComparison

dates="20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL"
minMafs="0.12 0.2 0.05 0.025"
for angsdDate in $dates
do

GLdir=$wd/angsd-GLs/$angsdDate 
PCAdir=$wd/pca/covarianceMatrices/$angsdDate

refs="mfur elut"
mkdir -p $PCAdir

for state in 1e-06.snpsOnly.TransvOnly

do
for ref in $refs
do
echo $ref
for minMaf in $minMafs
do
echo $minMaf
python $pcangsddir/pcangsd.py \
-beagle $GLdir/angsdOut.mappedTo${ref}.${state}.beagle.gz \
-o $PCAdir/pcAngsd.${ref}.${state}.minMaf.${minMaf} \
-minMaf $minMaf -threads 10 1> $PCAdir/pcAngsd.${ref}.${state}.minMaf.${minMaf}.log

# this generates a covariance matrix called pcAngsd.ref.state.cov.npy which is a numpy binary file
done
done
done
done

source deactivate

# then have to move to python or R to deal with this
# using a modification of:
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R


# then can
# using a modification of:
# https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/plotPCA.R
# to plot the pca
# using script plot.pcAngsd.R
