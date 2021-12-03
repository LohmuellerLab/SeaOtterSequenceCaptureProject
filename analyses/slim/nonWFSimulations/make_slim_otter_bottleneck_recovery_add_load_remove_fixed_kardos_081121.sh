# Script to make SLIM job script
# USAGE: ./make_slim_otter_bottleneck_recovery_add_load_remove_fixed_DLhmix_010721.sh [Na] [Nb] [nF] [Nr]


# Set Na, the ancestral population size
Na=${1}

# Set Nb, the bottleneck population size
Nb=${2}

# Set nF, the number of founders
nF=${3}

# Set Nr, the recovery population size
Nr=${4}


# Make script
cat > slim_otter_bottleneck_recovery_add_load_remove_fixed_kardos_${Na}Na_${Nb}Nb_${nF}nF_${Nr}Nr_081121.slim << EOM

initialize() {

	initializeSLiMModelType("nonWF");
	defineConstant("K1", ${Na});
	defineConstant("K3", ${Nb});
	defineConstant("K4", ${Nr});
	defineConstant("num_founders", ${nF});
	defineConstant("sampleSize", 60);
	defineConstant("g",20000); //number of genes
	defineConstant("ROHcutoff", 1000000);
	defineConstant("geneLength", 1500);
	defineConstant("seqLength", g*geneLength);
	//cat("Genome length:"+seqLength+"\n");


	initializeMutationRate(8.64e-9);

	defineConstant("h_sublethal", 0.001);
	defineConstant("h_VstrDel", 0.07);
	defineConstant("h_strDel", 0.31);
	defineConstant("h_wkmodDel", 0.48);

	//draw deleterious mutations from human DFE
	//and allow for different dominance coefficients for mutations with different s
	//by creating different mutation types (faster than fitness callbacks)



        //lethals (s=-1)
        initializeMutationType("m1", 0.0, "f", -1.0);

	//sub lethals (s<-0.4)
	initializeMutationType("m2", h_sublethal, "s", "do x=rgamma(1,-0.05,0.5); while (x >= -0.4); return x;");

	//very strongly deleterious mutations (-0.4 <= s < -0.1)
	initializeMutationType("m3", h_VstrDel, "s", "do x=rgamma(1,-0.05,0.5); while (x < -0.4 | x >= -0.1); return x;");

	//strongly deleterious mutations (-0.01 > s >= -0.1)
	initializeMutationType("m4", h_strDel, "s", "do x=rgamma(1,-0.05,0.5); while (x < -0.1 | x >= -0.01); return x;");

	//weakly/moderately deleterious mutations (s >= -0.01)
	initializeMutationType("m5", h_wkmodDel, "s", "do x=rgamma(1,-0.05,0.5); while (x < -0.01); return x;");

	//neutral mutations
	initializeMutationType("m6", 0.5,"f",0.0);


	//ratio of different deleterious mutation types taken from Kardos 2021 DFE (sum to 100 below)
	//assume ratio of deleterious to neutral muts of 2.31:1 (Huber et al 2017)
	//giving 100/2.31=43.3 for neutral mutations below
	initializeGenomicElementType("g1", c(m1,m2,m3,m4,m5,m6), c(5,0.44,14.5,47.2,32.9,43.3));


	// convert fixed nearly neutral muts to substitutions during burn in
	m5.convertToSubstitution = T;
	m6.convertToSubstitution = T;



	// approach for setting up genes on different chromosomes adopted from Jacqueline's wolf scripts

	// vector of # genes on 38 different dog chromosomes (scaled by Jacqueline to 1000 total genes) - need to rescale according to number of desired genes
	gene_nums=c(56,39,42,40,40,35,37,34,28,31,34,33,29,28,29,27,29,25,24,26,23,28,24,21,23,18,21,19,19,18,18,18,14,19,12,14,14,11);
	gene_nums = gene_nums*g/1000; //need to scale to number of desired genes since above array was originally set up for 1000 genes


	for (i in 1:g){
		initializeGenomicElement(g1, ((i-1)*geneLength)+(i-1), (i*geneLength)+(i-2) );
	}


	rates=NULL;

	// Multiple chromosomes:
	for (i in 1:(size(gene_nums)-1)){
		rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[i-1]-1)), 0.5);
	}
	rates=c(rates, 0, rep(c(1e-3, 0), asInteger(gene_nums[size(gene_nums)-1]-1)));

	ends=NULL;
	for (i in 1:g){
		ends=c(ends, (i*geneLength)+(i-2), (i*geneLength)+(i-1));
	}
	ends=ends[0:(size(ends)-2)];

	initializeRecombinationRate(rates, ends);

}





// define function getStats that randomly samples a subpopulation for sampSize # of inds and outputs a string of:
// pop size, mean fitness, heterozygosity, mean Froh, and avg num of variants of different classes per individual (very str del, str del, mod del, wk del)

function (s) getStats(o pop, i sampSize)
{
	i = sample(pop.individuals, sampSize, F);

	m = sortBy(i.genomes.mutations, "position"); //get all mutations in sample
	m_uniq = unique(m); // get rid of redundant muts
	DAF = sapply(m_uniq, "sum(m == applyValue);"); // count number of each mut in pop
	m_uniq_polym = m_uniq[DAF != i.genomes.size()]; //remove fixed mutations??

	//initialize vectors
	ROH_length_sumPerInd = c();
	Num_VstrDel_muts = c();
	Num_strDel_muts = c();
	Num_modDel_muts = c();
	Num_wkDel_muts = c();
	ind_het = c();
	fitness_population = c();
	load_population = c();
	S_population = c();


	for (individual in i) {

		indm = sortBy(individual.genomes.mutations, "position");
		indm = indm[match(indm, m_uniq_polym) >= 0];   // Check that individual mutations are not fixed
		indm_uniq = unique(indm);

		genotype = sapply(indm_uniq, "sum(indm == applyValue);");

		// tally number of mutations for different classes of selection coefficient per individual
		s = individual.genomes.mutations.selectionCoeff;

		Num_VstrDel_muts = c(Num_VstrDel_muts, sum(s<=-0.05));
		Num_strDel_muts = c(Num_strDel_muts, sum(s<=-0.01));
		Num_modDel_muts = c(Num_modDel_muts, sum(s<=-0.001 & s > -0.01));
		Num_wkDel_muts = c(Num_wkDel_muts, sum(s<=-0.00001 & s > -0.001));

		if (isNULL(genotype)) {
			ind_het = c(ind_het, 0); //putting this here to avoid error when trying to sum null vector
			next;
		}

		ind_het = c(ind_het, sum(genotype==1)/(seqLength));

		//code for getting ROHs

		ID_het = (genotype == 1); //outputs T/F for genotypes if they are het or homDer
		ID_homDer = (genotype == 2);
		pos_het = indm_uniq.position[ID_het]; //outputs positions of heterozgoys genotypes

		startpos = c(0, pos_het); //adds 0 to beggining of vector of hets
		endpos = c(pos_het, sim.chromosome.lastPosition); //adds last position in genome to vector of hets
		pos_het_diff = endpos - startpos;
		ROH_startpos = startpos[pos_het_diff > ROHcutoff]; //filter out startpos that dont correspond to ROH > 1Mb
		ROH_endpos = endpos[pos_het_diff > ROHcutoff];
		ROH_length = pos_het_diff[pos_het_diff > ROHcutoff]; //vector of ROHs for each individual
		ROH_length_sum = sum(ROH_length);
		ROH_length_sumPerInd = c(ROH_length_sumPerInd, ROH_length_sum); // add sum of ROHs for each individual to vector of ROHs for all individuals

		// calculate fitness - additive across sites
		allmuts = c(individual.genomes[0].mutationsOfType(m1), individual.genomes[1].mutationsOfType(m1),individual.genomes[0].mutationsOfType(m2), individual.genomes[1].mutationsOfType(m2),individual.genomes[0].mutationsOfType(m3), individual.genomes[1].mutationsOfType(m3),individual.genomes[0].mutationsOfType(m4), individual.genomes[1].mutationsOfType(m4),individual.genomes[0].mutationsOfType(m5), individual.genomes[1].mutationsOfType(m5));
		uniquemuts = c(individual.uniqueMutationsOfType(m1),individual.uniqueMutationsOfType(m2),individual.uniqueMutationsOfType(m3),individual.uniqueMutationsOfType(m4),individual.uniqueMutationsOfType(m5));
		load_individual = c();
		S_individual = c();

		if (size(uniquemuts) > 0){
			for (u in uniquemuts){
				places = (allmuts.id == u.id);
				uu = allmuts[places];

				if (size(uu) == 2) {
					S = sum(uu.selectionCoeff);
				}
				 else if (size(uu) == 1) {
					if (u.mutationType == m1) {
						S = uu.selectionCoeff*0.0;
					}
					if (u.mutationType == m2) {
						S = uu.selectionCoeff*h_sublethal;
					}
                                        if (u.mutationType == m3) {
                                                S = uu.selectionCoeff*h_VstrDel;
                                        }
                                        if (u.mutationType == m4) {
                                                S = uu.selectionCoeff*h_strDel;
					}
                                        if (u.mutationType == m5) {
                                                S = uu.selectionCoeff*h_wkmodDel;
                                        }

				}
				S_individual = c(S_individual, S);
			}
			S_individual = sum(S_individual);
			S_population = c(S_population, S_individual);
		} else {
			S_population = c(S_population, 0);
		}
	S_total = sum(S_population)/(sampSize);
	mean_fitness = exp(S_total);
	load = 1-mean_fitness;
	}

	return(pop.individuals.size() + "," + load + "," + mean(ind_het) + "," + mean(ROH_length_sumPerInd)/seqLength + "," + mean(Num_VstrDel_muts) + "," + mean(Num_strDel_muts)+ "," + mean(Num_modDel_muts) + "," + mean(Num_wkDel_muts));
}



reproduction() {
	subpop.addCrossed(individual, subpop.sampleIndividuals(1));
}



1 early() {
	cat("gen,popSize,meanFitness,meanLoad,FROH,avgVstrDel,avgStrDel,avgModDel,avgWkDel" + "\n");
	sim.addSubpop("p1", 10);
}

1:$((${Na}*11)) early() {
	p1.fitnessScaling = K1 / p1.individualCount;
}

//track statistics pre-bottleneck every 1000 generations
1:$((${Na}*11)) late() {
		if (sim.generation % 1000 == 0) {
			stats = getStats(p1, sampleSize);
			cat(sim.generation + "," + stats + "\n");
		}
}


// stop removing fixed muts
$((${Na}*10)) early(){
	m5.convertToSubstitution = F;
	m6.convertToSubstitution = F;
}




// bottleneck to p3
$((${Na}*11+1)) early(){
	sim.addSubpop("p3",0);
	migrants = sample(p1.individuals, num_founders);
	p3.takeMigrants(migrants);
	cat("gen,K3,p_death,popSizeP3,meanLoad,meanHet,FROH,avgVStrDel,avgStrDel,avgModDel,avgWkDel" + "\n");

	sim.tag = K3; // use sim.tag to keep track of K3 from one generation to the next
}



// fitness scaling for p3 with carrying capacity K3

$((${Na}*11+1)):$((${Na}*11+51)) early() {
	p1.fitnessScaling = 0; // kill off p1

	// kill off individuals at random - not sure if I should then adjust the individualCount
	inds = p3.individuals;

	//simulate beta distribution
	alpha = 0.5;
	beta = 8;
	x1 = rgamma(1, mean = alpha, shape=alpha);
	x2 = rgamma(1, mean = beta, shape=beta);
	beta = x1/(x1+x2); //probability of stochastic mortality this generation


	//set probability of death for each generation equal to outcome of beta
	for(i in inds){
		kill = rbinom(1,1,beta);
		if(kill==1){
			i.fitnessScaling = 0.0;
		}
	}

	//OU model with phi=0.9
	sim.tag = asInteger(exp((1-0.9)*log(K3)+0.9*log(sim.tag)+rnorm(n = 1, mean = 0, sd = log10(1.3))));
	p3.fitnessScaling = sim.tag / p3.individualCount;

	cat(sim.generation + "," + sim.tag + "," + beta + ",");

}


// set sim.tag to K4
$((${Na}*10+52)) early(){
        sim.tag = K4; // use sim.tag to keep track of K3 from one generation to the next
}


// recovery
$((${Na}*11+52)):$((${Na}*11+400)) early() {
        p1.fitnessScaling = 0; // kill off p1

        // kill off individuals at random - not sure if I should then adjust the individualCount
        inds = p3.individuals;

        //simulate beta distribution
        alpha = 0.5;
        beta = 8;
        x1 = rgamma(1, mean = alpha, shape=alpha);
        x2 = rgamma(1, mean = beta, shape=beta);
        beta = x1/(x1+x2); //probability of stochastic mortality this generation


        //set probability of death for each generation equal to outcome of beta
        for(i in inds){
                kill = rbinom(1,1,beta);
                if(kill==1){
                        i.fitnessScaling = 0.0;
                }
        }

        //OU model with phi=0.9
        sim.tag = asInteger(exp((1-0.9)*log(K4)+0.9*log(sim.tag)+rnorm(n = 1, mean = 0, sd = log10(1.3))));
        p3.fitnessScaling = sim.tag / p3.individualCount;

        cat(sim.generation + "," + sim.tag + "," + beta + ",");

}



// track statistics for P3 every generation and terminate when the population goes to 1 individual or after 1000 generations
$((${Na}*11+1)):$((${Na}*11+400)) late() {
	if(p3.individuals.size() < 2){
		stats_P3 = c("NA,NA,NA,NA,NA,NA,NA,NA"); //cant get stats from just one individual
	}
	if(p3.individuals.size() < sampleSize & p3.individuals.size() > 1){	// case when p3 size is less than sample size but greater than 1
		stats_P3 = getStats(p3, p3.individuals.size());
	}
	if(p3.individuals.size() >= sampleSize){ //case when p3 size is greater than or equal to sample size
		stats_P3 = getStats(p3, sampleSize);
	}

	cat(stats_P3 + "\n");

	if(p3.individuals.size() < 2){
			sim.simulationFinished();
			cat("The population has gone extinct");
	}
}
EOM
