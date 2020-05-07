# ab changes: took out parallel; am actually adding the --bootstrap flag to treemix to make it resample the data over k blocks 
###########################################
##### Treemix pipeline with bootstrap #####
#License: GPLv3 or later
#Modification date: 2017-05-10
#Written by: Elia Vajana, Marco Milanesi
#Contact: Elia Vajana <vajana.elia@gmail.com>, Marco Milanesi <marco.milanesi.mm@gmail.com>
#Description: Run Treemix software for a specific migration number previously chosen applying a bootstrap procedure 
###########################################
infile=$1 		## treemix input file
numk=$2 		## the migration number you choose and you want test
ncore=$3 		## max number of cores to use
blockk=$4 		## block size
outgroup=$5 	## name of the selected outgroup population (if you want to do an unrooted ML tree put here 'NoOutgroup' (without quotes))
nboot=$6		## number of bootstrap replicates of the TREE
pathP=$7		## path to Phylip consense program. Example: /biosoftware/phylip/phylip-3.696/exe/consense
outdir=$8
outname=$9		## name for output file

##########################################################################################
#########################
##### Settings file #####
#########################
echo "Input file name = $1" > $outdir/$outname"_Settings.txt"
echo "Output file name = $8/$9" >> $outdir/$outname"_Settings.txt"
echo "Number of migration modeled = $2" >> $outdir/$outname"_Settings.txt"
echo "Number of SNPS per block = $4" >> $outdir/$outname"_Settings.txt"
if [ $outgroup = "NoOutgroup" ]; then
	echo "Unrooted ML tree" >> $outdir/$outname"_Settings.txt"
else
	echo "Outgroup = $5" >> $outdir/$outname"_Settings.txt"
fi
echo "Number of replications = $6" >> $outdir/$outname"_Settings.txt"

mkdir $outdir/bootstrap



#############################

#### BOOTSTRAP PROCEDURE ####
#############################
#AB: adding the -bootstrap flag to treemix to resample the data over k size blocks
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### N.B. If you need to modify the treemix parameters please do it here ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### Seed number settings ###
echo "**** RUNNING TREEMIX ****"
if [ $outgroup = "NoOutgroup" ]; then
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -global -bootstrap -i $2 -k $3 -se -m $4 -seed $c -o $outdir"/bootstrap/"$5"_treemix_bootrep_"$1
	}
	export -f dowork
	for i in `seq 1 $nboot`
	do
	echo "boot: " $i
	dowork $i $infile $blockk $numk $outname > $outdir/$outname"_logfile_treemix_bootrep.log"
	done
else
echo "outgroup is: " $outgroup
	dowork() { 
		a=$RANDOM
		b=1`date +%N`
		let "c=$a+$b"
		treemix -global -bootstrap -i $2 -k $3 -se -m $4 -root $6 -seed $c -o $outdir/"bootstrap/"$5"_treemix_bootrep_"$1
	}
	export -f dowork
	for i in `seq 1 $nboot`
	do
	echo "boot:" $i
	echo "dowork $i $infile $blockk $numk $outname $outgroup"
	dowork $i $infile $blockk $numk $outname $outgroup > $outdir/$outname"_logfile_treemix_bootrep.log"
	done
fi
echo "**** RUNNING TREEMIX: DONE ****"


### Create a file with all the bootstrapped trees
for a in `seq 1 $nboot`
do
	bootf=$outdir/"bootstrap/"$outname"_treemix_bootrep_"$a".treeout.gz"
	gunzip -c $bootf | head -1 >> $outdir/$outname"_boottree.tre"
done

echo "   ***** Boostrap procedure: DONE *****"



#########################################################################
#### Run Phylip on the bootstrapped trees to obtain a consensus tree ####
#########################################################################
echo "   ***** Phylip - consensus tree construction: START *****"
### Clean the environment
cd $outdir # will this work?
rm -rf outfile outtree screanout

# Create parameters file
if [ $outgroup = "NoOutgroup" ]; then
	echo $outdir/$outname"_boottree.tre" > $outdir/$outname"_PhylipInputFile"
	echo "Y" >> $outdir/$outname"_PhylipInputFile"
else
	# Find the position of Outgroup population
	posOutgroup=`head -1 $outdir/$outname"_boottree.tre" | tr "," "\n" | grep $outgroup -n | cut -d":" -f1`
	# echo $posOutgroup
	echo $outdir/$outname"_boottree.tre" > $outdir/$outname"_PhylipInputFile"
	echo "O" >> $outdir/$outname"_PhylipInputFile"
	echo $posOutgroup >> $outdir/$outname"_PhylipInputFile"
	echo "Y" >> $outdir/$outname"_PhylipInputFile"
fi

# Run Phylip
$pathP < $outdir/$outname"_PhylipInputFile" > $outdir/screanout

### The output from Phylip will be modified because:
### 1) Treemix accept only one line tree 
### 2) Treemix accept newick format file 
#sed ':a;N;$!ba;s/\n//g' outtree > $outdir/$outname"_outtree.newick"
cat $outdir/outtree | tr -d "\n" > $outdir/$outname"_outtree.newick"
echo >> $outdir/$outname"_outtree.newick"
echo "   ***** Phylip - consensus tree construction: DONE *****"


######################################################################################
### Run TreeMix with the chosen number of migrations by loading the consensus tree ###
######################################################################################
###AB: don't run with -bootstrap setting here (don't want to resample data):
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### N.B. If you need to modify the treemix parameters please do it here ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

echo "**** RUNNING TREEMIX with Consensus tree ****"
if [ $outgroup = "NoOutgroup" ]; then
	treemix -global -i $infile -m $numk -k $blockk -se -tf $outdir/$outname"_outtree.newick" -o $outdir/${outname}.tfConsensus > $outdir/$outname"_logfile_treemix_boot.log"
else
	treemix -global -i $infile -m $numk -k $blockk -root $outgroup -se -tf $outdir/$outname"_outtree.newick" -o $outdir/${outname}.tfConsensus > $outdir/$outname"_logfile_treemix_boot.log"
fi
echo "**** RUNNING TREEMIX: DONE ****"

echo "TreeMix - Bootstrap Analysis --> DONE"

