##### adding admixture to PC angsd #####
# tutorial http://www.popgen.dk/software/index.php/PCAngsdTutorial


module load anaconda # load anaconda
source activate angsd-conda-env # activate conda env
pcangsddir=/u/home/a/ab08028/klohmueldata/annabel_data/bin/pcangsd
wd=$SCRATCH/captures/aDNA-ModernComparison
todaysdate=`date +%Y%m%d`

################ high coverage ################
dates="20191212-highcov-AFprior-MajorMinor4-plusCOM-KUR-AL" # with more individuals
refs="elut mfur" # for now
minMafs="0.025 0.05 0.12 0.2"

state=1e-06.snpsOnly.TransvOnly
alpha=50 # soft uppper limit on alpha, a small penalty on admixture proportions

for date in $dates
do
echo $date
GLdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/angsd-GLs/$date
outdir=/u/flashscratch/a/ab08028/captures/aDNA-ModernComparison/admixture/$date
mkdir -p $outdir

for ref in $refs
do
echo $ref
input=angsdOut.mappedTo${ref}.${state}.beagle.gz
for minMaf in $minMafs
do
echo $minMaf
python $pcangsddir/pcangsd.py \
-beagle $GLdir/$input \
-o $outdir/pcAngsd.${ref}.${state}.minMaf.${minMaf} \
-minMaf $minMaf -threads 10  \
-admix -admix_auto $alpha 1> $outdir/pcAngsd.${ref}.${state}.minMaf.${minMaf}.autoAlpha.${todaysdate}.log
1> saves the output so you can see that things converge, etc. useful! 
# this gives 2-3 groupings (that make sense along a gradient)
# you're not supposed to set K but I want to try to get more groups; going to adjust "e" (number of eigenvalues) to 5
# for the 5 pops (CA, AK, AL, COM, KUR)
python $pcangsddir/pcangsd.py \
-beagle $GLdir/$input \
-o $outdir/pcAngsd.${ref}.${state}.minMaf.${minMaf}.e5 \
-minMaf $minMaf -threads 10  \
-e 5 \
-admix -admix_auto $alpha 1> $outdir/pcAngsd.${ref}.${state}.minMaf.${minMaf}.autoAlpha.e5.${todaysdate}.log


done
done
done


source deactivate
