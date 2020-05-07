### want to better automate these scripts for different trim values (can do in R easily enough I think)

# dadi function for sso msmc model; trim point = 33; Nanc=22759
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot
import math
############### LL funcs ##########
########## Multinomial Log  Likelihood ##########
# writing numBins here instead of numHaps because for folded SFS you want
# number of bins not number of haps.
def LhoodCalc(model_SFS_freq,obs_SNP_counts,numBins):
    llSum=0
    for i in range(1,numBins):
        llSum += obs_SNP_counts[i]*math.log(model_SFS_freq[i])
    return llSum

########## POISSON Log LIKELIHOODS (see Lohmueller 2008) ##########

def LhoodCalcPoisson(model_SFS_count,obs_SNP_counts,numBins):
    llSum=0
    for i in range(1,numBins):
        llSum += obs_SNP_counts[i]*math.log(model_SFS_count[i])
    ll = -np.sum(model_SFS_count) + llSum
    return ll
    
################# set up dirs #####################
pop="AL" # want to try with two pops
sfsdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/"
outdir="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/analysisResults/StitchPSMC_SFS_Together/"+pop # outdir
mu=8.64e-09
CA_L=5989967 # from /Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/CA-13.totalSiteCount.L.withMonomorphic.txt


############### Input data ####################################
sfs=sfsdir+"CA-13.plusMonomorphic.sfs" # get CA SFS 
fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS
# fold the fs:
#fs=fs.fold() # folded
# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
    

############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10

########## get LL of data:data based on frequency of data (fs/sum(fs) and fs)
multinom_data2data=LhoodCalc(fs/sum(fs),fs,ns/2) # ns/2 is diploid SS (number of filled bins) 
################### Copy in parameters from the model file you output in R ############

Nanc = 5467
##
scalingTheta=4*Nanc*mu*CA_L # this is the theta to scale the SFS by
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.36675344664418, 0.141422067280315, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968, 0.0881712197999992)



def sso_model_trim_27((nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, T0, nu0)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	phi = Integration.one_pop(phi, xx, T2, nu2)
	phi = Integration.one_pop(phi, xx, T3, nu3)
	phi = Integration.one_pop(phi, xx, T4, nu4)
	phi = Integration.one_pop(phi, xx, T5, nu5)
	phi = Integration.one_pop(phi, xx, T6, nu6)
	phi = Integration.one_pop(phi, xx, T7, nu7)
	phi = Integration.one_pop(phi, xx, T8, nu8)
	phi = Integration.one_pop(phi, xx, T9, nu9)
	phi = Integration.one_pop(phi, xx, T10, nu10)
	phi = Integration.one_pop(phi, xx, T11, nu11)
	phi = Integration.one_pop(phi, xx, T12, nu12)
	phi = Integration.one_pop(phi, xx, T13, nu13)
	phi = Integration.one_pop(phi, xx, T14, nu14)
	phi = Integration.one_pop(phi, xx, T15, nu15)
	phi = Integration.one_pop(phi, xx, T16, nu16)
	phi = Integration.one_pop(phi, xx, T17, nu17)
	phi = Integration.one_pop(phi, xx, T18, nu18)
	phi = Integration.one_pop(phi, xx, T19, nu19)
	phi = Integration.one_pop(phi, xx, T20, nu20)
	phi = Integration.one_pop(phi, xx, T21, nu21)
	phi = Integration.one_pop(phi, xx, T22, nu22)
	phi = Integration.one_pop(phi, xx, T23, nu23)
	phi = Integration.one_pop(phi, xx, T24, nu24)
	phi = Integration.one_pop(phi, xx, T25, nu25)
	phi = Integration.one_pop(phi, xx, T26, nu26)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs

func=sso_model_trim_27
modelName="sso_model_trim_27"
params=params # change for each model 
# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
func_ex = Numerics.make_extrap_log_func(func)


model=func_ex(params,ns,pts_l) # this is relative to theta =1
### get proportional sfs 
model_freq = model/(sum(model))
model_freq_fold =model_freq.fold()


############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(modelName)+".PROPORTIONAL.FOLDED.expSFS.txt"

model_freq_fold.to_file(outputSFS)


######### get and output LLs ##################
##### get multinomial ll from proportional model SFS and count obs sfs:
### do my own LL calc, not dadis (which optimizes model to fit data)
outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()
########## plot an image: ############
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(modelName)+".expSFS.DadiScaling.figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
pyplot.title((modelName))
plt.savefig(outputFigure)



###################### Try without that semi-recent dip ###############
# converted that nu31 interval that was 0.141422067280315 to  0.386870175412936 (same as previous interval)
params_noDip = (0.36675344664418, 0.386870175412936, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968, 0.0881712197999992)


modelName="sso_model_trim_27_noDip"
params=params_noDip
###### put in params no dip
model=func_ex(params,ns,pts_l) # this is relative to theta =1
### get proportional sfs 
model_freq = model/(sum(model))
model_freq_fold =model_freq.fold()


############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(modelName)+".PROPORTIONAL.FOLDED.expSFS.txt"

model_freq_fold.to_file(outputSFS)


######### get and output LLs ##################
##### get multinomial ll from proportional model SFS and count obs sfs:
### do my own LL calc, not dadis (which optimizes model to fit data)
outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()
########## plot an image: ############
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(modelName)+".expSFS.DadiScaling.figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
pyplot.title((modelName))
plt.savefig(outputFigure)


########################## Try to add the contraction ###########
 

contractionGen=35 # 35 gen
contractionSize=195 # diploids 
# rescale by nanc:
contractionGen_RescBy2Nanc = float(contractionGen)/(2*Nanc)
contractionSize_RescByNanc = float(contractionSize)/Nanc

params_plusContraction=(0.36675344664418, 0.141422067280315, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968,0.0881712197999992,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

## adding them to end of parameter vector:
def sso_model_trim_27_plusContraction((nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0,contractionSize_RescByNanc,contractionGen_RescBy2Nanc),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, T0, nu0)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	phi = Integration.one_pop(phi, xx, T2, nu2)
	phi = Integration.one_pop(phi, xx, T3, nu3)
	phi = Integration.one_pop(phi, xx, T4, nu4)
	phi = Integration.one_pop(phi, xx, T5, nu5)
	phi = Integration.one_pop(phi, xx, T6, nu6)
	phi = Integration.one_pop(phi, xx, T7, nu7)
	phi = Integration.one_pop(phi, xx, T8, nu8)
	phi = Integration.one_pop(phi, xx, T9, nu9)
	phi = Integration.one_pop(phi, xx, T10, nu10)
	phi = Integration.one_pop(phi, xx, T11, nu11)
	phi = Integration.one_pop(phi, xx, T12, nu12)
	phi = Integration.one_pop(phi, xx, T13, nu13)
	phi = Integration.one_pop(phi, xx, T14, nu14)
	phi = Integration.one_pop(phi, xx, T15, nu15)
	phi = Integration.one_pop(phi, xx, T16, nu16)
	phi = Integration.one_pop(phi, xx, T17, nu17)
	phi = Integration.one_pop(phi, xx, T18, nu18)
	phi = Integration.one_pop(phi, xx, T19, nu19)
	phi = Integration.one_pop(phi, xx, T20, nu20)
	phi = Integration.one_pop(phi, xx, T21, nu21)
	phi = Integration.one_pop(phi, xx, T22, nu22)
	phi = Integration.one_pop(phi, xx, T23, nu23)
	phi = Integration.one_pop(phi, xx, T24, nu24)
	phi = Integration.one_pop(phi, xx, T25, nu25)
	phi = Integration.one_pop(phi, xx, T26, nu26)
 	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


func=sso_model_trim_27_plusContraction
modelName="sso_model_trim_27_plusContraction"
params=params_plusContraction # change for each model 
# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
func_ex = Numerics.make_extrap_log_func(func)


model=func_ex(params,ns,pts_l) # this is relative to theta =1
### get proportional sfs 
model_freq = model/(sum(model))
model_freq_fold =model_freq.fold()


############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(modelName)+".PROPORTIONAL.FOLDED.expSFS.txt"

model_freq_fold.to_file(outputSFS)


######### get and output LLs ##################
##### get multinomial ll from proportional model SFS and count obs sfs:
### do my own LL calc, not dadis (which optimizes model to fit data)
outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()
########## plot an image: ############
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(modelName)+".expSFS.DadiScaling.figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
pyplot.title((modelName))
plt.savefig(outputFigure)

########################## try with No dip and dadi contraction ######
# converted that nu31 interval that was 0.0339730033926217 to  0.0929355795355233 (same as previous interval)
# converted that nu31 interval that was 0.141422067280315 to  0.386870175412936 (same as previous interval)
params_plusContraction_noDIP = (0.36675344664418, 0.386870175412936, 0.386870175412936, 0.609784494435127, 0.774150824662306, 0.91011303614748, 1.01928801965413, 1.09213414017289, 1.12514450221579, 1.12274330685933, 1.07023813979649, 1.07023813979649, 0.968579772770968, 0.968579772770968, 0.869604885707006, 0.869604885707006, 0.795226999398677, 0.795226999398677, 0.750983120146509, 0.750983120146509, 0.738213027247672, 0.738213027247672, 0.761765214135537, 0.761765214135537, 0.836551538729166, 0.836551538729166, 1, 0.030119771118, 0.030902139933, 0.0317262985630002, 0.0325958441059994, 0.0335133156900004, 0.0344855901200002, 0.0355149949300004, 0.0366078779399999, 0.0377705869699996, 0.03900946984, 0.0403308743700001, 0.0417474961999999, 0.0432646251800001, 0.0448991888299998, 0.04665859294, 0.0485650548800005, 0.0506312702899994, 0.0528826304500005, 0.05534346867, 0.0580444660799995, 0.0610215936600004, 0.0643224600599991, 0.0679989058100032, 0.072122872869994, 0.0767768829000064, 0.0820773125999968,0.0881712197999992,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)


func=sso_model_trim_27_plusContraction
modelName="sso_model_trim_27_noDip_plusContraction"
params=params_plusContraction_noDIP # change for each model 
# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
func_ex = Numerics.make_extrap_log_func(func)


model=func_ex(params,ns,pts_l) # this is relative to theta =1
### get proportional sfs 
model_freq = model/(sum(model))
model_freq_fold =model_freq.fold()


############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(modelName)+".PROPORTIONAL.FOLDED.expSFS.txt"

model_freq_fold.to_file(outputSFS)


######### get and output LLs ##################
##### get multinomial ll from proportional model SFS and count obs sfs:
### do my own LL calc, not dadis (which optimizes model to fit data)
outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()
########## plot an image: ############
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(modelName)+".expSFS.DadiScaling.figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
pyplot.title((modelName))
plt.savefig(outputFigure)

############## pulling from grid search, this is the best-fit sfs from dadi with T = 35 for CA ##########
modelName="bestFitDadiModel.T35.fromGridSearch"
model = dadi.Spectrum([0,0.7390820 ,0.5346896, 0.4275530 ,0.3692772, 0.3400618 ,0.1655816,0,0,0,0,0,0]).fold() # this is from R in my grid search for CA ; is w/in 1 pt of MLE with T = 35 gen
model_freq_fold = model/sum(model)

outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

####################### highly simplified MSMC model: just two epochs ###############
### going to go from nanc =1 (4500~) to nu = 2000/4500 = 0.4, for 1000 gen (1000/(2*4500)) = 0.1 [4500 is weighted Ne average of MSMC curve ]
# what if I fix T1 and do inference?

params_2Epoch=(0.4,0.1)
def sso_model_trim_27_simplify((nu1, T1),ns,pts):
	nu1,T1=params
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs

func=sso_model_trim_27_simplify
modelName="sso_model_trim_27_simplify"
params=params_2Epoch # change for each model 
# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
func_ex = Numerics.make_extrap_log_func(func)


model=func_ex(params,ns,pts_l) # this is relative to theta =1
### get proportional sfs 
model_freq = model/(sum(model))
model_freq_fold =model_freq.fold()
############### Output SFS ########################
print('Writing out SFS **************************************************')                                   

outputSFS=str(outdir)+"/"+str(modelName)+".PROPORTIONAL.FOLDED.expSFS.txt"

model_freq_fold.to_file(outputSFS)


######### get and output LLs ##################
##### get multinomial ll from proportional model SFS and count obs sfs:
### do my own LL calc, not dadis (which optimizes model to fit data)
outputFile=open(str(outdir)+"/"+str(modelName)+".LLs.andOptimalTheta.txt","w")
multinom_LL_AB= LhoodCalc(model_freq_fold,fs,ns/2)
dadi_ll_msmc_model = dadi.Inference.ll_multinom(model, fs )
optimalthetaFromDadi = dadi.Inference.optimal_sfs_scaling(model, fs) # 
header='\t'.join(str(x) for x in ("dadiLL","AnnabelLL","NancTheta","dadiOptimalTheta"))
output='\t'.join(str(x) for x in (dadi_ll_msmc_model,multinom_LL_AB,scalingTheta,optimalthetaFromDadi))
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()
########## plot an image: ############
#import pylab
import matplotlib.pyplot as plt 
#fig=plt.figure(1)
#pylab.ion()
outputFigure=str(str(outdir)+"/"+str(modelName)+".expSFS.DadiScaling.figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)
pyplot.title((modelName))
plt.savefig(outputFigure)