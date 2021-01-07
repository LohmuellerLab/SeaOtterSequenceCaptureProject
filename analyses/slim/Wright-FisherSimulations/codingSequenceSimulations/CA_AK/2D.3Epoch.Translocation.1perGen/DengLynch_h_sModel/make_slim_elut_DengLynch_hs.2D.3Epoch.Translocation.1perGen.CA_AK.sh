pop=CA_AK
model=2D.3Epoch.Translocation.1perGen
#model=1D.2Epoch.1.5Mb.cds.LongerContract
gitdir=/u/home/p/pkalhori/project-klohmueldata/pooneh_data/github_repos/otter_exome/SLIM
scriptdir=$gitdir/slim_scripts/$pop/$model
mkdir -p $scriptdir # set script dir
todaysdate=`date +%Y%m%d`
########## population specific parameters ########
# set sample size for vcf (should match empirical for AK)
# can set manually or pull from a table
ss_CA=6
ss_AK=7 # in diploids # will depend on population, note you can find these in projectionValues.txt file
#population	ProjectionValueHaploids Diploids
#KUR	12	6
#CA	12	6
#COM	34	17
#AK	14	7
#AL	20	10
# variables:
nanc_combo=4500
nanc_AK=4500 # a ball park Nanc for AK
nu_AK=250
nanc_CA=3500
nu_CA=200
 # contraction size  -- a medium contraction size
tcontract=35 # contraction duration before you sample -- a longer contraction time than dadi ifnerence (more like FSC)
trecovery=17 #1911-2019 is around 17 generations (6y/gen)
tfuture=50 #allow simulation to continue running at recovered pop size  
nrec_AK=2500 #not sure what population size to use; best I could find for WPWSx
nrec_CA=1000
# this nu / T combo is within the MLE estimate for Alaska from my grid search. also want to try the California version (100 for 30 gens with smaller Nanc of ~3500)
######### general parameters ; can set here or in command line ##########


migAtoC=1e-3
migCtoA=4e-4

# Set g, number of genes (exons)
g=1000
# Set gLen, the length of exons
gLen=1500
# Set t, number of burn-in generations
t=50000
tdiv=4000
# set mutation rate
mu=8.64e-9 # mutation rate
# Set h, dominance coefficient
h=$1  # loop through hs
# Set j, the chunk number (for 14 chunks?)

chr=20
#total number of chromosomes, since all "chunks" will be simulated together 


# Make script
# chunk gets set when you run slim (based on SGE task id) # or something?? how to do this part? I don't really want to make a separte slim script each time? have it be a -d thing maybe?
# have to figure out the chunks/replicates situation.

cat > $scriptdir/slim_elut_${model}_${pop}_hs.job << EOM

// changes to make: apparently 1e-03 is reasonable between-gene recomb rate
// and then want to make separate chromosomes with 0.5 between them (or just simulate them separately)
// okay I think I am going to model 1.5Mb stretches of seqeunce, each containing 1000 genes of size 1500bp
// recomb w/in each gene will be 1e-08, between genes will be 1e-3
// will either do 14x1.5mb so they are independent, or will simulate them all together if want to think about overall genetic load (then would have to set up chromosomes within the simulation)
initialize() {
	defineConstant("g",$g); //number of genes; starting with 1000 (AB)
	defineConstant("geneLength", $gLen); // length of each gene
	defineConstant("seqLength", g*geneLength); // total chromosome length starting with 1.5Mb (AB)
	//defineConstant("outdir",\"$outdir\"); -- set in command line
	//defineConstant("v_CHUNK",$chunk); // portion of genome  -- set in command line
	//defineConstant("v_REP",$rep); // overall replicate number -- set in command line
	defineConstant("v_h",$h); // dominance coefficient
	defineConstant("v_SS_CA",$ss_CA); // sample size
	defineConstant("v_MU",$mu);
	defineConstant("v_NANC_CA",$nanc_CA); // ancestral size for CA
	defineConstant("v_NU_CA",$nu_CA); // contraction size California
	defineConstant("v_NREC_CA",$nrec_CA); // recovery size CA
	defineConstant("v_SS_AK",$ss_AK); // sample size
	defineConstant("v_NANC_AK",$nanc_AK); // ancestral size
	defineConstant("v_NU_AK",$nu_AK); // contraction size
	defineConstant("v_NREC_AK",$nrec_AK); // contraction size
	defineConstant("v_NCOMBO",$nanc_combo); // pre split time
	defineConstant("migAK_CA",$migAtoC); // migration to CA from AK
	defineConstant("migCA_AK",$migCtoA); // migration to AK from CA
	

	//cat("Exome portion length:"+seqLength+"\n");
	initializeMutationRate(v_MU);
	// m1 mutation type: neutral *synonymous*
	initializeMutationType("m1", 0.5, "f", 0.0);
	// m2 mutation type: missense(?) -- This is from Chris; ask where he got params (Kim et al?)
	initializeMutationType("m2", v_h, "g",-0.01314833, 0.186); // set H in array (0 for recessive, 0.5 for additive) 
	m2.convertToSubstitution = F; // okay cool if you have this, then fixed sites are in the vcf file
	m1.convertToSubstitution = F; // keeps fixed sites in vcf file (do I want?)
	// initialize exon: g1, has both neutral (m1) and misense (m2) mutations
	// syn happens at rate 1, mis at rate 2.31:1 since there are more missense sites (Christian/Chris); note chris names neutral as m2 and missense at m1 so I reversed things here 
	initializeGenomicElementType("g1", c(m1,m2), c(1.0,2.31)); // 2.31 is the NS:S ratio
	
	for (i in 0:(g-1)){
		initializeGenomicElement(g1, ((i)*geneLength+1), ((i+1)*geneLength));
	}
	// figured out a way for recomb to be 0.5 between blocks, but 1e-08 within blocks
	// so that they aren't actually linked; simulating  independent genes/exons
	//initializeRecombinationRate(1e-8);
	// this sets up rates that alternate between 1e-08 and 1e-03
	// and ends that are in pattern 999 1000 1999 2000 2999 ...
	// so that the pattern is that r is 1e-8 for 0-999, then 1e-03 (a reasonable between-gene recomb rate) between 999 and 1000 (each gene), then is 1e-8 through next gene, and so on.
	// making 5000 independent blocks.
	rates=c(rep(c(1e-08,1e-3),g));
	ends=NULL;
	for (index in 0:(g-1))
	{
		ends=c(ends,index*geneLength+(geneLength-1),index*geneLength+geneLength);
	}
	initializeRecombinationRate(rates,ends);


}
1: fitness(m2) {
// this is from Deng and Lynch: 
h = 0.5 * exp(-13*abs(mut.selectionCoeff))
if (homozygous) {
    // 20210107: this was Bernard's code to deal with slight excess heterosis in his Plos Genet paper; we don't need this; initial set of revisions was run with this, but then we re-ran without it for final submission --> return ((1.0 + 0.5*mut.selectionCoeff)*(1.0 + 0.5*mut.selectionCoeff));
    // 20210107: now am calculating homozgyous derived fitness the same as in the other simulations: 
    return (1.0 + mut.selectionCoeff)
} else {
	return (1.0 + mut.selectionCoeff * h);
}
}
// create a population of variable v_NANC individuals
1 {
	sim.addSubpop("p1", v_NCOMBO);
}

// output generation number so I can track progress

1:${t} late() {
	if (sim.generation % 1000 == 0){
		cat(sim.generation+"\n");
	}
}

//Sample After burn in 

${t} late() {
			// set up outfile: 
	writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),"replicate,chunk,generation,mutid,type,s,age,originpop,subpop,numhet,numhom,popsizeDIP\n",append=F); // open fresh file
	//file header
	//mutation id
	//mutation type
	//selection coefficient
	//age of mutation in generations
	//subpopulation it arose in
	//number of heterozygote derived in p1
	//number of homozygote derived in p1
	//number of heterozygote derived in p2
	//number of homozygote derived in p2
	//these are genotype counts not allele counts

	
	//for every mutation in the simulation
	//pops=sim.subpopulations;
	for (pop in sim.subpopulations){
		for (mut in sim.mutations){
			id = mut.id;
			s = mut.selectionCoeff;
			generation= sim.generation - 50000;
			originpop = mut.subpopID;
			age = sim.generation - mut.originGeneration;
			type = mut.mutationType;
			popsize = size(pop.individuals);
			popID= pop.id;
			//initialize genotype counts
			pnumhet = 0;
			pnumhom = 0;
			
			//count hom and het derived in p1
			for (p1i in pop.individuals){
				gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
				if (gt == 1){
					pnumhet = pnumhet + 1;
				} else if (gt == 2){
					pnumhom = pnumhom + 1;
				}
			}
					// string for mutation type. add m3, m4, etc. if you have multiple types
			if (type == m1){
				type = "m1";
			} else if (type == m2){
				type = "m2";
			}
			//print results
		writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),paste(c(v_REP,v_CHUNK,generation,id,type,s,age,originpop,popID,pnumhet,pnumhom,popsize),sep=","),append=T);
		}
	}
	}

//Burn in 50,000 generations, then split pops. Use CA params for p2 Resize p1 to the AK parameters
$((${t} + 1))  early() {
	sim.addSubpopSplit("p2",v_NANC_CA,p1);
	p1.setSubpopulationSize(v_NANC_AK);

}

//after burn in, calculate load every 1000 generations until contraction 
$((${t} + 1)):$((${t} + 1+ ${tdiv})) late() {
	if (sim.generation % 1000 == 0){
	//file header
	//mutation id
	//mutation type
	//selection coefficient
	//age of mutation in generations
	//subpopulation it arose in
	//number of heterozygote derived in p1
	//number of homozygote derived in p1
	//number of heterozygote derived in p2
	//number of homozygote derived in p2
	//these are genotype counts not allele counts

	
	//for every mutation in the simulation
	//pops=sim.subpopulations;
	for (pop in sim.subpopulations){
		for (mut in sim.mutations){
			id = mut.id;
			s = mut.selectionCoeff;
			generation= sim.generation - 50000;
			originpop = mut.subpopID;
			age = sim.generation - mut.originGeneration;
			type = mut.mutationType;
			popsize = size(pop.individuals);
			popID= pop.id;
			//initialize genotype counts
			pnumhet = 0;
			pnumhom = 0;
			
			//count hom and het derived in p1
			for (p1i in pop.individuals){
				gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
				if (gt == 1){
					pnumhet = pnumhet + 1;
				} else if (gt == 2){
					pnumhom = pnumhom + 1;
				}
			}
					// string for mutation type. add m3, m4, etc. if you have multiple types
			if (type == m1){
				type = "m1";
			} else if (type == m2){
				type = "m2";
			}
			//print results
		writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),paste(c(v_REP,v_CHUNK,generation,id,type,s,age,originpop,popID,pnumhet,pnumhom,popsize),sep=","),append=T);
		}
	}
}
}


//Sample before contraction
$((${t} + 1+ ${tdiv})) late() {
	p1.outputVCFSample(v_SS_AK, F,filePath=paste(c(outdir,"/slim.output.PreContraction.p1.",v_CHUNK,".vcf"),sep=""));
	p2.outputVCFSample(v_SS_CA, F,filePath=paste(c(outdir,"/slim.output.PreContraction.p2.",v_CHUNK,".vcf"),sep=""));
}


//after burn in and precontraction, calculate load every 2 generations until endsim  
$((${t} + 2+ ${tdiv})): late() {
	if (sim.generation % 2 == 0){
	//file header
	//mutation id
	//mutation type
	//selection coefficient
	//age of mutation in generations
	//subpopulation it arose in
	//number of heterozygote derived in p1
	//number of homozygote derived in p1
	//number of heterozygote derived in p2
	//number of homozygote derived in p2
	//these are genotype counts not allele counts

	
	//for every mutation in the simulation
	//pops=sim.subpopulations;
	for (pop in sim.subpopulations){
		for (mut in sim.mutations){
			id = mut.id;
			s = mut.selectionCoeff;
			generation= sim.generation - 50000;
			originpop = mut.subpopID;
			age = sim.generation - mut.originGeneration;
			type = mut.mutationType;
			popsize = size(pop.individuals);
			popID= pop.id;
			//initialize genotype counts
			pnumhet = 0;
			pnumhom = 0;
			
			//count hom and het derived in p1
			for (p1i in pop.individuals){
				gt = sum(c(p1i.genomes[1].containsMutations(mut), p1i.genomes[0].containsMutations(mut)));
				if (gt == 1){
					pnumhet = pnumhet + 1;
				} else if (gt == 2){
					pnumhom = pnumhom + 1;
				}
			}
					// string for mutation type. add m3, m4, etc. if you have multiple types
			if (type == m1){
				type = "m1";
			} else if (type == m2){
				type = "m2";
			}
			//print results
		writeFile(paste(c(outdir,"/slim.output.",v_CHUNK,".summary.txt"),sep=""),paste(c(v_REP,v_CHUNK,generation,id,type,s,age,originpop,popID,pnumhet,pnumhom,popsize),sep=","),append=T);
		}
	}
}
}


//Contract after 4000 generations 

$((${t} + 2+ ${tdiv})) late() {
	p1.setSubpopulationSize(v_NU_AK);
	p2.setSubpopulationSize(v_NU_CA);
}


//Sample before recovery
$((${t} + 2+ ${tdiv} + ${tcontract})) late() {
	p1.outputVCFSample(v_SS_AK, F,filePath=paste(c(outdir,"/slim.output.PreRecovery.p1.",v_CHUNK,".vcf"),sep=""));
	p2.outputVCFSample(v_SS_CA, F,filePath=paste(c(outdir,"/slim.output.PreRecovery.p2.",v_CHUNK,".vcf"),sep=""));
}

//Recover populatoin sizes after (35/25) generations. Tcontract is 35, but CA started 10 generations later

$((${t} + 3+ ${tdiv} + ${tcontract})) late() {
	p1.setSubpopulationSize(v_NREC_AK);
	p2.setSubpopulationSize(v_NREC_CA);
	}

//Present Day is here! Should be around 14 gens later

//Sample before setting up migration rates
$((${t} + 3+ ${tdiv} + ${tcontract}+ ${trecovery})) late() {
	p1.outputVCFSample(v_SS_AK, F,filePath=paste(c(outdir,"/slim.output.PreMigration.p1.",v_CHUNK,".vcf"),sep=""));
	p2.outputVCFSample(v_SS_CA, F,filePath=paste(c(outdir,"/slim.output.PreMigration.p2.",v_CHUNK,".vcf"),sep=""));
}
$((${t} + 4 + ${tdiv} + ${tcontract}+ ${trecovery})) late() {
		p1.setMigrationRates(c(p2), c(migCA_AK));
		p2.setMigrationRates(c(p1), c(migAK_CA));
		}
//Let the simulation run into the future, then sample at the ends
$((${t} + 4+ ${tdiv} + ${tcontract}+ ${trecovery}+${tfuture})) late() {
	p1.outputVCFSample(v_SS_AK, F,filePath=paste(c(outdir,"/slim.output.EndSimulation.p1.",v_CHUNK,".vcf"),sep=""));
	p2.outputVCFSample(v_SS_CA, F,filePath=paste(c(outdir,"/slim.output.EndSimulation.p2.",v_CHUNK,".vcf"),sep=""));
}

$((${t} + 5 + ${tdiv} + ${tcontract}+ ${trecovery}+${tfuture})) late() 
{
sim.simulationFinished();
} 
EOM



