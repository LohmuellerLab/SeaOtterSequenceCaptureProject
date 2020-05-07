# SeaOtterSequenceCaptureProject
Scripts related to data processing and analyzing sequence capture of 122 sea otters from across the species range. 
These scripts were used to generate the analyses for Beichman et al. in prep.
Details of the methods are in the supplementary materials of the paper.
Contact: Annabel Beichman <annabel.beichman[at]gmail.com>
Note: these scripts are internal to this project, based on my file systems and directory organization. 
Feel free to adapt and reuse in your own work with the caveat that they are not currently designed for generic use, so results are not guaranteed.

## data_processing
Scripts related to read mapping, genotype calling and filtering for modern and ancient data
	
	* paleomixPipeline: maps reads (ancient and modern)

	* variant_calling: find covered variants and call genotypes (GATK)
	
	* variant_filtering: filter variants on depth and quality (GATK)
	
	* variant_calling_aDNA: call and filter aDNA variants and variants in 15 modern samples (3 from each population) for comparison
	
## analyes
Scripts to carry out structure and demographic analyses

	* aDNA-ModernComparison: carry out PCA and structure analyses on ancient and modern data using genotype likelihoods
	
	* calculate_pi: calculate and plot pi and Watterson's Theta per popualtion from the neutral SFS for empirical and simluated data
	
	* dadi_inference: infer demographic history using dadi 
	
	* fastsimcoal_inference: infer demographic history using fastsimcoal2
	
	* FASTSTRUCTURE: infer structure of populations 
	
	* fst: calculate fst using SNPRelate
	
	* generate_sfs: generate the projected site frequency spectrum using EasySFS
	
	* PCA: carry out PCA analysis using SNPRelate
	
	* RELATEDNESS: calculate kinship using SNPRelate
	
	* slim: forward-in-time simulations using SLiM
	
	* StitchMSMC_SFS_Together: combine MSMC models from Beichman et al. (2019) with SFS-based demographic inference
	
	* TREEMIX: run treemix using the BITE wrapper
