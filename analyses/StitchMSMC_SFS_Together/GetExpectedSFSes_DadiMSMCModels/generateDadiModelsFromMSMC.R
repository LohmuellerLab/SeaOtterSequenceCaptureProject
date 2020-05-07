require(ggplot2)
require(scales)
require(RColorBrewer)
require(scales)
elutCol="#377EB8"
elut2Col="#C77CFF"
mu=8.64e-09
elutTimeIndices=c(39,33,27,22)
elut2TimeIndex=c(39,31,27,20) # i like 20 for this. 
# these were chosen to match the pbra curve , so maybe want different ones now. Can make different models in a for loop for any different trim points which is nice. (set up script that way so it's relatively easy.)
# start with SSO and let's see what happens.
################ generate southern sea otter (elut) model #############
for(elutTimeIndex in elutTimeIndices){

sso = read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingGidgetReadsToFerret/MSMC_DemographicInference_WholeGenome/output_20180209_50_250DPFilter/elut.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",header=T)
# northern sea otter (AK/AL)

# trim: 

# rescale for dadi (sizes scaled by Nanc, times by 2*Nanc)
sso_trim = sso[sso$time_index<elutTimeIndex,]
# select the columsn that correspond to desired 'reasonable' mu (8.64e-09)
sso_trim_reasonable = sso_trim[,c("time_index","left_time_boundary","right_time_boundary","lambda_00","Ne.reasonable","Left_generations.reasonable","Right_generations.reasonable")]
# need to make dadi time intervals first:
sso_trim_reasonable$time_interval_gen <- sso_trim_reasonable$Right_generations.reasonable - sso_trim_reasonable$Left_generations.reasonable

# Choose Nanc (oldest)
sso_Nanc=sso_trim_reasonable[sso_trim_reasonable$time_index==max(sso_trim_reasonable$time_index),]$Ne.reasonable
# will need to reverse their order as well.
# Divide by 2*Nanc
sso_trim_reasonable$time_interval_dadiUnits <- sso_trim_reasonable$time_interval_gen / (2*sso_Nanc)
# divide by Nanc to get NE in dadi units:
# checkthat last interval has size = 1 (correct)
sso_trim_reasonable$Ne.DadiUnits <- sso_trim_reasonable$Ne.reasonable/sso_Nanc
# think more of how I want to do this . want nuA (ancestral) to be 0. 

#### count up intervals: will be same as time index (but since my numbering is going to zero I want index-1 as highest value)
# put together list of interval names
# want them to be nubmered backward where most recent interval is nuHighestNumber and oldest interval is nuA (ancestral) 
### note that nu0 is ancestral size (not calling it nuA anymore.)
variableNames_nu = c(paste("nu",seq(elutTimeIndex-1,0,-1),sep="")) # should be length time index in length (matched nu and T pairs) (counting from 0)
variableNames_T = c(paste("T",seq(elutTimeIndex-1,0,-1),sep=""))
# want to write out function for dadi.
outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/StitchPSMC_SFS_Together/"
###### write out dadi model (can eventually do for multiple trim points) ############
# header:
sso_filename=paste(outdir,"southern.sea.otter.msmc2dadi.model.trim.",elutTimeIndex,".Nanc.",round(sso_Nanc),".py",sep="")
write(file = sso_filename,paste("# dadi function for sso msmc model; trim point = ",elutTimeIndex,"; Nanc=",round(sso_Nanc),sep="")) # initialize file
write(file = sso_filename,"# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.\n# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.\n",append=T)
write(file=sso_filename,"from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference\nimport numpy as np\nimport matplotlib.pyplot as pyplot\n\n",append=T)
######### need to write out params in same order that they shoudl go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0)
# add Nanc so I can use it to scale my contraction?
write(file=sso_filename,paste("Nanc = ",round(sso_Nanc),sep=""),append=T)
write(file=sso_filename,"# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.",append=T)
write(file=sso_filename,paste("params = (",toString(sso_trim_reasonable$Ne.DadiUnits),", ",toString(sso_trim_reasonable$time_interval_dadiUnits),")",sep=""),append=T)
write(file=sso_filename,"\n\n",append=T)
write(file = sso_filename,paste("def sso_model_trim_",elutTimeIndex,"((",toString(variableNames_nu),", ",toString(variableNames_T),"),ns,pts):",sep=""),append=T)
### write rest of model:
write(file=sso_filename,"\txx = Numerics.default_grid(pts)",append=T)
write(file=sso_filename,"\t# intialize phi with ancestral pop: (nu0 = 1)",append=T)
write(file=sso_filename,"\tphi = PhiManip.phi_1D(xx,nu=nu0)",append=T)
write(file=sso_filename,"\t# stays at nu0 for T0 duration of time:",append=T)
write(file=sso_filename,"\tphi = Integration.one_pop(phi,xx,T0,nu0)",append=T)
write(file=sso_filename,"\t# followed by a number of time steps, with associated pop changes:",append=T)
# write out each pair of variables, going in reverse order from how they're listed so that nu0 (ANCESTRAL) comes first, ends with nu32 or whatever the highest is.
# note that this one line writes out each pair as its own line going through the vector (works nicely)
write(file=sso_filename,paste("\tphi = Integration.one_pop(phi, xx, ",rev(variableNames_T),", ",rev(variableNames_nu),")",sep=""),append=T)
# finally, get expected sfs:
write(file=sso_filename,"\t# get expected SFS:",append=T)
write(file=sso_filename,"\tfs = Spectrum.from_phi(phi,ns,(xx,))",append=T)
write(file=sso_filename,"\treturn fs\n\n",append = T)
# make extrapoloation version of function
write(file=sso_filename,"# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function",append=T)
write(file=sso_filename,paste("extrap_sso_model_trim_",elutTimeIndex,"_function = Numerics.make_extrap_log_func(sso_model_trim_",elutTimeIndex,")",sep=""),append=T)


######### write version of model that has values for nu and T plugged in and an extra contraction epoch for dadi to infer: #######

write(file=sso_filename,"\n\n# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize ",append=T)
write(file = sso_filename,paste("def sso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):",sep=""),append=T)
### write rest of model:
write(file=sso_filename,"\txx = Numerics.default_grid(pts)",append=T)
write(file=sso_filename,"\t# intialize phi with ancestral pop: (nu0 = 1)",append=T)
# set as 1 because nu0 will always be one (everything in terms of it)
write(file=sso_filename,"\tphi = PhiManip.phi_1D(xx,nu=1)",append=T)
write(file=sso_filename,"\t# stays at nu0=1 for T0 duration of time:",append=T)
write(file=sso_filename,"\t# followed by a number of time steps, with associated pop changes:",append=T)
# write out each pair of variables, going in reverse order from how they're listed so that nu0 (ANCESTRAL) comes first, ends with nu32 or whatever the highest is.
# note that this one line writes out each pair as its own line going through the vector (works nicely)
write(file=sso_filename,paste("\tphi = Integration.one_pop(phi, xx, ",rev(sso_trim_reasonable$time_interval_dadiUnits),", ",rev(sso_trim_reasonable$Ne.DadiUnits),")",sep=""),append=T)
### then add a contraction period:
write(file=sso_filename,"\t### add contraction:",append=T)
write(file=sso_filename,"\tphi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)",append=T)
# finally, get expected sfs:
write(file=sso_filename,"\t# get expected SFS:",append=T)
write(file=sso_filename,"\tfs = Spectrum.from_phi(phi,ns,(xx,))",append=T)
write(file=sso_filename,"\treturn fs\n\n",append = T)
# make extrapoloation version of function
write(file=sso_filename,"# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function",append=T)
write(file=sso_filename,paste("extrap_sso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization)",sep=""),append=T)
}

################################## northern sea otter ###################
nso=read.table("/Users/annabelbeichman/Documents/UCLA/Otters/RawReadsToGenotypes_Dec2016/mappingNorthernSeaOtterReadsToFerret/MSMC_50_250DPFilter_simsInSppDIrs/output_20190215/elut2.MSMC.results.allScalings.10_14_20_MyaDivergence.txt",header=T)

for(elutTimeIndex in elut2TimeIndices){

  # northern sea otter (AK/AL)
  
  # trim: 
  
  # rescale for dadi (sizes scaled by Nanc, times by 2*Nanc)
  nso_trim = nso[nso$time_index<elutTimeIndex,]
  # select the columsn that correspond to desired 'reasonable' mu (8.64e-09)
  nso_trim_reasonable = nso_trim[,c("time_index","left_time_boundary","right_time_boundary","lambda_00","Ne.reasonable","Left_generations.reasonable","Right_generations.reasonable")]
  # need to make dadi time intervals first:
  nso_trim_reasonable$time_interval_gen <- nso_trim_reasonable$Right_generations.reasonable - nso_trim_reasonable$Left_generations.reasonable
  
  # Choose Nanc (oldest)
  nso_Nanc=nso_trim_reasonable[nso_trim_reasonable$time_index==max(nso_trim_reasonable$time_index),]$Ne.reasonable
  # will need to reverse their order as well.
  # Divide by 2*Nanc
  nso_trim_reasonable$time_interval_dadiUnits <- nso_trim_reasonable$time_interval_gen / (2*nso_Nanc)
  # divide by Nanc to get NE in dadi units:
  # checkthat last interval has size = 1 (correct)
  nso_trim_reasonable$Ne.DadiUnits <- nso_trim_reasonable$Ne.reasonable/nso_Nanc
  # think more of how I want to do this . want nuA (ancestral) to be 0. 
  
  #### count up intervals: will be same as time index (but since my numbering is going to zero I want index-1 as highest value)
  # put together list of interval names
  # want them to be nubmered backward where most recent interval is nuHighestNumber and oldest interval is nuA (ancestral) 
  ### note that nu0 is ancestral size (not calling it nuA anymore.)
  variableNames_nu = c(paste("nu",seq(elutTimeIndex-1,0,-1),sep="")) # should be length time index in length (matched nu and T pairs) (counting from 0)
  variableNames_T = c(paste("T",seq(elutTimeIndex-1,0,-1),sep=""))
  # want to write out function for dadi.
  outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/scripts/scripts_20180521/analyses/StitchPSMC_SFS_Together/"
  ###### write out dadi model (can eventually do for multiple trim points) ############
  # header:
  nso_filename=paste(outdir,"northern.sea.otter.msmc2dadi.model.trim.",elutTimeIndex,".Nanc.",round(nso_Nanc),".py",sep="")
  write(file = nso_filename,paste("# dadi function for nso msmc model; trim point = ",elutTimeIndex,"; Nanc=",round(nso_Nanc),sep="")) # initialize file
  write(file = nso_filename,"# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.\n# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.\n",append=T)
  write(file=nso_filename,"from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference\nimport numpy as np\nimport matplotlib.pyplot as pyplot\n\n",append=T)
  ######### need to write out params in same order that they shoudl go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0)
  # add Nanc so I can use it to scale my contraction?
  write(file=nso_filename,paste("Nanc = ",round(nso_Nanc),sep=""),append=T)
  write(file=nso_filename,"# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.",append=T)
  write(file=nso_filename,paste("params = (",toString(nso_trim_reasonable$Ne.DadiUnits),", ",toString(nso_trim_reasonable$time_interval_dadiUnits),")",sep=""),append=T)
  write(file=nso_filename,"\n\n",append=T)
  write(file = nso_filename,paste("def nso_model_trim_",elutTimeIndex,"((",toString(variableNames_nu),", ",toString(variableNames_T),"),ns,pts):",sep=""),append=T)
  ### write rest of model:
  write(file=nso_filename,"\txx = Numerics.default_grid(pts)",append=T)
  write(file=nso_filename,"\t# intialize phi with ancestral pop: (nu0 = 1)",append=T)
  write(file=nso_filename,"\tphi = PhiManip.phi_1D(xx,nu=nu0)",append=T)
  write(file=nso_filename,"\t# stays at nu0 for T0 duration of time:",append=T)
  write(file=nso_filename,"\tphi = Integration.one_pop(phi,xx,T0,nu0)",append=T)
  write(file=nso_filename,"\t# followed by a number of time steps, with associated pop changes:",append=T)
  # write out each pair of variables, going in reverse order from how they're listed so that nu0 (ANCESTRAL) comes first, ends with nu32 or whatever the highest is.
  # note that this one line writes out each pair as its own line going through the vector (works nicely)
  write(file=nso_filename,paste("\tphi = Integration.one_pop(phi, xx, ",rev(variableNames_T),", ",rev(variableNames_nu),")",sep=""),append=T)
  # finally, get expected sfs:
  write(file=nso_filename,"\t# get expected SFS:",append=T)
  write(file=nso_filename,"\tfs = Spectrum.from_phi(phi,ns,(xx,))",append=T)
  write(file=nso_filename,"\treturn fs\n\n",append = T)
  # make extrapoloation version of function
  write(file=nso_filename,"# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function",append=T)
  write(file=nso_filename,paste("extrap_nso_model_trim_",elutTimeIndex,"_function = Numerics.make_extrap_log_func(nso_model_trim_",elutTimeIndex,")",sep=""),append=T)
  
  
  ######### write version of model that has values for nu and T plugged in and an extra contraction epoch for dadi to infer: #######
  
  write(file=nso_filename,"\n\n# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize ",append=T)
  write(file = nso_filename,paste("def nso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):",sep=""),append=T)
  ### write rest of model:
  write(file=nso_filename,"\txx = Numerics.default_grid(pts)",append=T)
  write(file=nso_filename,"\t# intialize phi with ancestral pop: (nu0 = 1)",append=T)
  # set as 1 because nu0 will always be one (everything in terms of it)
  write(file=nso_filename,"\tphi = PhiManip.phi_1D(xx,nu=1)",append=T)
  write(file=nso_filename,"\t# stays at nu0=1 for T0 duration of time:",append=T)
  write(file=nso_filename,"\t# followed by a number of time steps, with associated pop changes:",append=T)
  # write out each pair of variables, going in reverse order from how they're listed so that nu0 (ANCESTRAL) comes first, ends with nu32 or whatever the highest is.
  # note that this one line writes out each pair as its own line going through the vector (works nicely)
  write(file=nso_filename,paste("\tphi = Integration.one_pop(phi, xx, ",rev(nso_trim_reasonable$time_interval_dadiUnits),", ",rev(nso_trim_reasonable$Ne.DadiUnits),")",sep=""),append=T)
  ### then add a contraction period:
  write(file=nso_filename,"\t### add contraction:",append=T)
  write(file=nso_filename,"\tphi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)",append=T)
  # finally, get expected sfs:
  write(file=nso_filename,"\t# get expected SFS:",append=T)
  write(file=nso_filename,"\tfs = Spectrum.from_phi(phi,ns,(xx,))",append=T)
  write(file=nso_filename,"\treturn fs\n\n",append = T)
  # make extrapoloation version of function
  write(file=nso_filename,"# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function",append=T)
  write(file=nso_filename,paste("extrap_nso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(nso_model_trim_",elutTimeIndex,"_plusContraction_forOptimization)",sep=""),append=T)
}
