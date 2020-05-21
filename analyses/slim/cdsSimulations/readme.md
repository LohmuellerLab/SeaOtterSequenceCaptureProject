This directory contains the models simulated using SLiM.

Models:
AK/1D.5Epoch: 5 Epoch model for Alaskan Population Parameters, with two consecutive contractions
CA/1D.3Epoch:3 Epoch model for the California population parameters, with one contraction followed by a recovery period

CA_AK/2D.3Epoch.{Model}: Two population migration model, with a single ancestral population splitting into two populations, following inferred Alaska and Calfornia Parameters. After a period of recovery, various migration rates are introduced for 50 generations:

No Translocation= No migration is introduced
1perGen= Migration rate of approximately 1 individual per generation for both populations
5perGen= Approximately 5 individuals per generation
10perGen= Approximately 10 individuals per generation
25perGen= Approximately 25 individuals per generation
25for2Gen= Approximately 2 individuals migrate from each population for two generations, then migration is stopped.

The make scripts for these models are found in the corresponding directories, along with the most recent SLiM job scripts containing the the most recently used parameters

This directory also contains the following files:

array_slim_elut.generic.sh: This script creates an array to submit an array of 20 simulations in parallel
 
slim.maker.submitter.allPopsAllModels.sh: This submits 25 replicates of the array script, for each corresponding model and dominance coefficient

concatenateResults.sh: This script concatenates the summary files for the 20 chunks per replicate

concatenateVCF.GenerateSFSes.sh: This script concatenates the VCFs from the 20 chunks, and generates SFSes with all mutations, as well as SFSes that are separated by neutral (mut1) and missense (mut2) mutations

FigureOutWhichRunsPassed.summary.sh: This script is run from the summary file concatenating script, and ensures that at least 20 viable replicates were completed

FigureOutWhichRunsPassed.vcf.sh: This is the same as the above script, but used when concatenating VCFs

calculateLoadFromSlimSims.ExcludeFixedBurninVariants.new.R: calculate genetic load from slim simulations, excluding variants that fix prior to the contraction

compareSimualatedEmpiricalSFSes.poonehSImulations.FORMANSUCRIPT.R: compare the simulated SFSes from SLiM to the empirical SFSes
