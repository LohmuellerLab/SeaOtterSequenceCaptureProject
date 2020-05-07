## Analyses
### These scripts were used to carry out analyses for the paper

#### aDNA-ModernComparison
	Scripts to carry out PCA and structure analyses on ancient and modern data using genotype likelihoods

#### calculate_pi
	Scripts to calculate and plot pi and Watterson's Theta per popualtion from the neutral SFS for empirical and simluated data

#### dadi_inference
	Scripts to carry out demographic inference in dadi for the 1-Dimensional 1Epoch, 2Epoch and 3Epoch models and the 2-Dimensinoal Isolation-Migration models described in the paper.
	* 1DModels: single population models. Wrapper script (dadi.Wrapper.Generic.ANYMODEL.Hoffman.sh). Models: 1Epoch (1D.1Epoch.dadi.py), 2Epoch (1D.2Epoch.dadi.py), 3Epoch (1D.1Bottleneck.dadi.py without bottleneck duratin fixed, 1D.1Bottleneck.TB20gen.dadi.py with bottleneck duration fixed at 0.005)
	* 2DModels: joint California-Alaska models, either a split with constant sized populations, with and without migration (2D.ConstantSize.Migration.dadi.py, 2D.ConstantSize.noMigration.dadi.py), or a split+contraction with and without migration (2D.Bottleneck.Migration.dadi.py,2D.Bottleneck.noMigration.dadi.py). Wrapper script (dadi.Wrapper.2D.Models.AllPopPairs.sh)
	* grid.search: carry out grid search along a log10-spaced grid of parameters for 2Epoch model (grid.Search.1D.2Epoch.dadi.py, dadi.grid.search.Wrapper.sh) or 3Epoch (grid.Search.1D.3Epoch.dadi.dadiUnits.py, grid.Search.1D.3Epoch.dadi.dadiUnits.FixTFOnly.py) with or without fixing bottleneck duration. Can extract sets of parameters (extract.RangeofValidParams.GridSearchResults.R) and plot (plot)
	* msmcModelsInDadi: get expected SFS from MSMC models from Beichman et al. (2019)
	* Likelihood.Ratio.Test.ComparingModels.R: carry out LRT to compare models

#### fastsimcoal_inference
	Scripts to carry out fastsimcoal inference
	
#### FASTSTRUCTURE
	Scripts to carry out and plot fastSTRUCTURE analysis
	
#### fst
	Scripts to calculate fst using SNPRelate and plot fst vs distance 
	
#### generate_sfs
	Scripts to generate the projected SFS using a modification of EasySFS 

#### PCA
	Scripts to carry out PCA analysis using SNPRelate

#### RELATEDNESS
	Scripts to calculate kinship using SNPRelate

#### slim
	Scripts to carry out SLiM simulations
	* neutralSimulations: Wright-Fisher simulations of neutral data
	* cdsSimulations: Wright-Fisher simulations of coding sequence
	* nonWFSimulations: non-Wright Fisher simulations

#### StitchMSMC_SFS_Together
	Scripts to combine MSMC models from Beichman et al. (2019) with SFS-based demographic inference

#### TREEMIX
	Scripts to run treemix using the BITE wrapper