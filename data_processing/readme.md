## Data Processing: Fastqs --> Genotypes
### These scripts will take you through initial QC, the paleomix mapping pipeline, variant calling and filtering


#### 1. paleomixPipeline: maps reads (ancient and modern)
	Makefiles to run paleomix
    * modernMakefiles: makefiles for running modern samples through paleomix mapping pipeline (without indel realignment because IR is done during GATK haplotype caller)

    * modernMakefiles-withIndelR-toCompareWithAncient: Makefiles to run paleomix on a subset of 15 modern samples with IR to compare to ancient samples in ANGSD

    * ancientMakefiles-fullProcessing-usethis: Makefiles to map ancient reads using paleomix

	**Paleomix Notes:** 
	Paleomix Documentation: https://paleomix.readthedocs.io/en/latest/
	Paleomix carries out AdapterRemoval, read mapping, mapping DNA damage and correcting it, and validation.
	It can be used with multiple reference genomes and multiple libraries per sample.
	First, AdapterRemoval trims adapters and low quality bases at the ends of reads. It also collapses overlapping reads.
	For ancient samples, I only use these collapsed reads downstream. Modern samples retain all reads.
	For ancient samples, bwa backtrack with the seed disabled and settings from Kircher's protocol is used to map (-n = 0.01; -o=2; min qual 30)
	For modern samples, bwa mem is used to map (min qual 30)
	All samples have mapDamage plots produced, but the mapDamage is only rescaled for ancient samples.
	Indel realignment is not carried out for modern samples because GATK Haplotype Caller does it internally, but is carried out for samples that are going to be processed in ANGSD (ancient + subset of modern)

#### 2. variant_calling: find covered variants and call genotypes (GATK)

	a. Step_2_a_FindCoveredIntervals.sh (submit with wrapper Step_2_a_FindCoveredIntervals.[modern/ancient].submit.sh): Detect covered intervals using GATK's FindCoveredIntervals to use downstream (min cov. = 1 read; MAPQ min. 30; min base qual. 20)
	
	b. Step_2_b_QualimapIntervals.sh (submit with wrapper Step_2_b_QualimapIntervals.submit.sh): Run Qualimap on covered regions (use multiqc to gather reports)  

	c. Step_2_c_HaplotypeCaller.sh (submit with Step_2_c_HaplotypeCaller.submit.sh): call individual variants using GATK's HaplotypeCaller to generate one g.vcf file per sample
		
	d. Step_3_genotypeGVCFS.modern.sh: call genotypes jointly for all populations together (exclude some extremely low coverage individuals)

#### 3.  variant_filtering: filter variants on depth and quality (GATK)

    To set hard-filtering thresholds, the distribution of QD (quality-by-depth) and DP (depth) were extracted from the genotype VCF file and plotted. Alternate alleles that did not appear in any individuals after joint genotype calling were trimmed. Sites with a minimum of 500 DP across all individuals were selected as an initial filter to remove sites with a high degree of missing data. Biallelic SNPs were selected and filtered using the following hard-filters: FS > 60.0, MQ < 40.0, MQRankSum < -12.5, ReadPosRankSum < -8.0, SOR > 3.0, as recommended in the GATK hard-filtering documentation. Sites with QD < 8.0 (determined by looking at the distribution of QD scores) were also removed. Filters per individual genotype were also applied, with genotypes with a genotype quality < 20 or a per-individual read depth <10 being marked as missing (“./.”).
    Sites that became non-variant after genotype filtering were removed from the SNP VCF file and added to the invariant site VCF file. Clusters of SNPs in which 3 SNPs appeared in a window of 10bp were removed. 
    Invariant sites were filtered to exclude sites with QUAL < 30, and to set genotypes with RGQ < 1 or per-individual DP < 10 as missing. 
    After filtering, the invariant and variant VCF files were combined and the amount of missing data per individual was calculated. Individuals that had more than one standard deviation above the mean missingness were excluded from downstream analysis.
    After removal of 19 individuals with a high degree of missing data, a final custom filtering script was applied to check for any remaining oddities in the data (script available on project Github). The script checked for reference and alternate alleles that were longer than a single letter, unexpected genotypes, sites without QUAL scores, sites that were not labeled as “PASS” after filtering, and sites missing DP, AN, GT, AD, DP or GQ/RGQ annotations. The script then updates the AN and AC annotations for each site based on the final set of genotype calls afterwards. All remaining biallelic SNPs from this filtered dataset were then selected for downstream analysis. A 20% per-site missingness filter was applied and the subsequent SNP set was used for relatedness and population structure analyses (below). 
	Steps:
	* filtering_Step_1_a_filterVariantAndInvariantSites.sh : filter variant and invariant sites
	* filtering_Step_1_b_RemoveBadIndividuals_BespokeFilters.sh : remove individuals that are poor coverage and carry out final filteirng checks
	* filtering_Step_1_c_i_ConvertToGDS.R (and filtering_Step_1_c_instructions.txt) : convert to GDS so you can perform kinship and structure analyses
	* filtering_Step_1_d-i_RemoveRelativesAdmixed.sh : remove close relative or admixed individuals 
	* filtering_Step_1_d-ii_getBedCoordsOfFinalSites.sh : get coordinates of final sites passing all filters
	* filtering_Step_1_e-i_IdentifyNeutralIntervals.sh : identify putatively neutral regions
	* filtering_Step_1_e-ii_pullOutNeutralVariantsFromVCF.forEasySFS.sh : extract putativley neutral regions
	* filtering_Step_1_e-iii_pullCDSFromVCF.forEasySFS.justSG.TESTNOPICK.sh : extract coding regions based on domestic ferret annotation
	* optional: filtering_Step_3_a_downSampleCommandersForStructurePCA.20181119.sh, filtering_Step_3_ConvertToGDS.generic.Hoffman.wrapper.sh, filtering_Step_3_ConvertToGDS.generic.Hoffman.R) (did it for the manucript): downsample populations with high sample sizes for structure analyses
	Helper scripts:
	* extract_QD_DP_QUAL.py: extract QD, DP and QUAL info for choosing filter cutoffs
	* filtering_bespokeFiltersAndChecks.py: modification of a script by Clare Marsden. A final set of filter checks to look for commong GATK bugs
	* filtering_getNoCallPerInd.py and filtering_getNoCallPerInd_wrapper.sh: count up called and missing sites per individual 
	* filtering_perPopulation.noCall.maxHetFilter.py: script to filter on max fraction of heterozygous genotypes per population

#### 4. variant_calling_aDNA: call and filter aDNA variants and variants in 15 modern samples (3 from each population) for comparison
	Used ANGSD with the following settings:
	GL 2 -doMajorMinor 4 -doMaf 1 -beagleProb 1 -doPost $posterior -remove_bads 1 -uniqueOnly 1 -C 50 -baq 1 -trim 7 -minQ 20 -minMapQ 25 -doGlf 2 -doCounts 1 -dumpCounts 2 -SNP_pval 1e-06 -rmTrans 1, which refer to the following: 
	* GL 2: Use GATK style likelihoods
	* doMajorMinor 4: pre-specify major allele based on reference genome (the ‘major’ allele may therefore not necessarily be the highest frequency allele, but will consistently be the reference allele to maintain formatting consistency with the original .vcf file)
	* doMaf 1: calculate allele frequencies with major/minor fixed based on reference genome
	* beagleProb 1: calculate beagle format posterior probabilities
	* doPost 1: do genotype posterior likelihood calculation, using allele frequencies calculated from the data as the prior
	* remove_bads 1: remove bad reads
	* uniqueOnly 1: discard reads that don’t map uniquely
	* baq 1: adjust base quality around indels
	* trim 7: trim the first and last 7 bp from each read (where misincorporation is likely to occur based on MapDamage plots) 
	* minQ 25: min base quality
	* minMapQ: min mapping quality
	* doGLF2: output beagle format genotype likelihoods file
	* doCounts 1 -dumpCounts 2: output coverage information for each individual and site
	* SNP_pval 1e-06: output SNPs with p-value < 1e-06
	* rmTrans 1: remove transitions, leaving only transversions which are less likely to be affected by ancient DNA damage
	Samples used are found in the bamLists/

