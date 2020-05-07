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
pop="CA"
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


Nanc = 4036
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
##
scalingTheta=4*Nanc*mu*CA_L # this is the theta to scale the SFS by

params = (0.496812482450453, 0.191573518835868, 0.52406305650732, 0.82602781572228, 1.0486821501222, 1.23285962527743, 1.38075051784768, 1.47942951405878, 1.52414609426596, 1.52089338093277, 1.44976869859196, 1.44976869859196, 1.31205998407016, 1.31205998407016, 1.17798637196801, 1.17798637196801, 1.07723241130487, 1.07723241130487, 1.01729865557432, 1.01729865557432, 1, 1, 0.04080092061, 0.041860735035, 0.0429771588850002, 0.0441550648699991, 0.0453978925500004, 0.0467149574000002, 0.0481094123500004, 0.0495898562999998, 0.0511648881499993, 0.0528431067999999, 0.0546331111500001, 0.0565520989999997, 0.0586072361000001, 0.0608214528499997, 0.0632047812999999, 0.0657873176000006, 0.0685862595499991, 0.0716360027500006, 0.0749695096499999, 0.0786283415999992, 0.0826612257000005, 0.0871326536999987)



def sso_model_trim_22((nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	phi = Integration.one_pop(phi,xx,T0,nu0)
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
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs

func=sso_model_trim_22
modelName="sso_model_trim_22"
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
# converted that nu31 interval that was 0.191573518835868 to  0.52406305650732 (same as previous interval)
params_noDip =  (0.496812482450453, 0.52406305650732, 0.52406305650732, 0.82602781572228, 1.0486821501222, 1.23285962527743, 1.38075051784768, 1.47942951405878, 1.52414609426596, 1.52089338093277, 1.44976869859196, 1.44976869859196, 1.31205998407016, 1.31205998407016, 1.17798637196801, 1.17798637196801, 1.07723241130487, 1.07723241130487, 1.01729865557432, 1.01729865557432, 1, 1, 0.04080092061, 0.041860735035, 0.0429771588850002, 0.0441550648699991, 0.0453978925500004, 0.0467149574000002, 0.0481094123500004, 0.0495898562999998, 0.0511648881499993, 0.0528431067999999, 0.0546331111500001, 0.0565520989999997, 0.0586072361000001, 0.0608214528499997, 0.0632047812999999, 0.0657873176000006, 0.0685862595499991, 0.0716360027500006, 0.0749695096499999, 0.0786283415999992, 0.0826612257000005, 0.0871326536999987)


modelName="sso_model_trim_22_noDip"
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
 
###### you are here ########
contractionGen=35 # 35 gen
contractionSize=195 # diploids 
# rescale by nanc:
contractionGen_RescBy2Nanc = float(contractionGen)/(2*Nanc)
contractionSize_RescByNanc = float(contractionSize)/Nanc

params_plusContraction=(0.496812482450453, 0.191573518835868, 0.52406305650732, 0.82602781572228, 1.0486821501222, 1.23285962527743, 1.38075051784768, 1.47942951405878, 1.52414609426596, 1.52089338093277, 1.44976869859196, 1.44976869859196, 1.31205998407016, 1.31205998407016, 1.17798637196801, 1.17798637196801, 1.07723241130487, 1.07723241130487, 1.01729865557432, 1.01729865557432, 1, 1, 0.04080092061, 0.041860735035, 0.0429771588850002, 0.0441550648699991, 0.0453978925500004, 0.0467149574000002, 0.0481094123500004, 0.0495898562999998, 0.0511648881499993, 0.0528431067999999, 0.0546331111500001, 0.0565520989999997, 0.0586072361000001, 0.0608214528499997, 0.0632047812999999, 0.0657873176000006, 0.0685862595499991, 0.0716360027500006, 0.0749695096499999, 0.0786283415999992, 0.0826612257000005, 0.0871326536999987,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)

## adding them to end of parameter vector:
def sso_model_trim_22_plusContraction((nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0,contractionSize_RescByNanc,contractionGen_RescBy2Nanc,  ),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=nu0)
	# stays at nu0 for T0 duration of time:
	phi = Integration.one_pop(phi,xx,T0,nu0)
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
 	### add contraction:
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


func=sso_model_trim_22_plusContraction
modelName="sso_model_trim_22_plusContraction"
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
# converted that nu31 interval that was0.191573518835868 to  0.52406305650732 (same as previous interval)

params_plusContraction_noDIP = (0.496812482450453, 0.52406305650732, 0.52406305650732, 0.82602781572228, 1.0486821501222, 1.23285962527743, 1.38075051784768, 1.47942951405878, 1.52414609426596, 1.52089338093277, 1.44976869859196, 1.44976869859196, 1.31205998407016, 1.31205998407016, 1.17798637196801, 1.17798637196801, 1.07723241130487, 1.07723241130487, 1.01729865557432, 1.01729865557432, 1, 1, 0.04080092061, 0.041860735035, 0.0429771588850002, 0.0441550648699991, 0.0453978925500004, 0.0467149574000002, 0.0481094123500004, 0.0495898562999998, 0.0511648881499993, 0.0528431067999999, 0.0546331111500001, 0.0565520989999997, 0.0586072361000001, 0.0608214528499997, 0.0632047812999999, 0.0657873176000006, 0.0685862595499991, 0.0716360027500006, 0.0749695096499999, 0.0786283415999992, 0.0826612257000005, 0.0871326536999987,contractionSize_RescByNanc,contractionGen_RescBy2Nanc)


func=sso_model_trim_22_plusContraction
modelName="sso_model_trim_22_noDip_plusContraction"
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


####################### highly simplified MSMC model: just two epochs ###############
### going to go from nanc =1 (4000~) to nu = 2000/4000 = 0.5, for 1000 gen (1000/(2*400)) = 0.125 [4500 is weighted Ne average of MSMC curve ]
# what if I fix T1 and do inference?

params_2Epoch=(0.5,0.125)
def sso_model_trim_22_simplify((nu1, T1),ns,pts):
	nu1,T1=params
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	phi = Integration.one_pop(phi, xx, T1, nu1)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs

func=sso_model_trim_22_simplify
modelName="sso_model_trim_22_simplify"
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