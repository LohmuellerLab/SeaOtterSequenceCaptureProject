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


Nanc = 411055
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
scalingTheta=4*Nanc*mu*CA_L # this is the theta to scale the SFS by

params = (0.00487809781986972, 0.00188102029958802, 0.00514566550504812, 0.00811059429737348, 0.0102967906221188, 0.0121051907162397, 0.0135573004479985, 0.014526208865738, 0.014965271610991, 0.0149333338993196, 0.0142349755244366, 0.0142349755244366, 0.0128828424685752, 0.0128828424685752, 0.01156640172282, 0.01156640172282, 0.0105771196632592, 0.0105771196632592, 0.00998864265534717, 0.00998864265534717, 0.00981879077556433, 0.00981879077556433, 0.0101320526482532, 0.0101320526482532, 0.0111267672454692, 0.0111267672454692, 0.0133007552199023, 0.0133007552199023, 0.0179070175840453, 0.0179070175840453, 0.0283778878420236, 0.0283778878420236, 0.0553680897108008, 0.0553680897108008, 0.137643054169846, 0.137643054169846, 0.447459632029814, 0.447459632029814, 1, 0.000400615702920001, 0.00041102179902, 0.000421983731220002, 0.000433549343639992, 0.000445752408600005, 0.000458684392800002, 0.000472376254200005, 0.000486912423599999, 0.000502377331799994, 0.0005188554096, 0.000536431087800002, 0.000555273227999998, 0.000575452189200002, 0.000597193120199998, 0.0006205945236, 0.000645951907200007, 0.000673434132599993, 0.000703378923000007, 0.0007361099298, 0.000772035235199994, 0.000811633280400006, 0.000855537296399988, 0.000904436801400043, 0.00095928867779992, 0.00102119052600009, 0.00109169024399996, 0.00117274381199999, 0.00126646200000007, 0.00137692562999991, 0.00150807480600003, 0.00166722686399998, 0.00186366919200004, 0.00211288077000004, 0.00243920581199994, 0.00288500043600006, 0.00353089605600002, 0.00455208658199998, 0.00641589649199996, 0.0109679830740001)



def sso_model_trim_39((nu38, nu37, nu36, nu35, nu34, nu33, nu32, nu31, nu30, nu29, nu28, nu27, nu26, nu25, nu24, nu23, nu22, nu21, nu20, nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T38, T37, T36, T35, T34, T33, T32, T31, T30, T29, T28, T27, T26, T25, T24, T23, T22, T21, T20, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
	phi = Integration.one_pop(phi, xx, T22, nu22)
	phi = Integration.one_pop(phi, xx, T23, nu23)
	phi = Integration.one_pop(phi, xx, T24, nu24)
	phi = Integration.one_pop(phi, xx, T25, nu25)
	phi = Integration.one_pop(phi, xx, T26, nu26)
	phi = Integration.one_pop(phi, xx, T27, nu27)
	phi = Integration.one_pop(phi, xx, T28, nu28)
	phi = Integration.one_pop(phi, xx, T29, nu29)
	phi = Integration.one_pop(phi, xx, T30, nu30)
	phi = Integration.one_pop(phi, xx, T31, nu31)
	phi = Integration.one_pop(phi, xx, T32, nu32)
	phi = Integration.one_pop(phi, xx, T33, nu33)
	phi = Integration.one_pop(phi, xx, T34, nu34)
	phi = Integration.one_pop(phi, xx, T35, nu35)
	phi = Integration.one_pop(phi, xx, T36, nu36)
	phi = Integration.one_pop(phi, xx, T37, nu37)
	phi = Integration.one_pop(phi, xx, T38, nu38)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs
func=sso_model_trim_39
modelName="sso_model_trim_39"
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

