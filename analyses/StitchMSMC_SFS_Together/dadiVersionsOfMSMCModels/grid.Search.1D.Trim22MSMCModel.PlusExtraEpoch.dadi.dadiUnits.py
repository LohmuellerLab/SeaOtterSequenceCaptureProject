# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 15:32:23 2018

@author: annabelbeichman
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018

@author: annabelbeichman
"""
import math # to get log10
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg') # so it doesn't pop up graphics on hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum, Demographics1D
# note: module  Demographics1D has the following models file:///Users/annabelbeichman/Documents/UCLA/dadi/dadi/doc/api/dadi.Demographics1D-module.html
# has 2 epoch, and 3 epoch (bottleneck)

from numpy import array # don't comment this out
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')

modelName="1D.Trim22.MSMCModel"
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Carry out a grid search for a '+ modelName +' model ')
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easysfs (mask optional)")
parser.add_argument("--nu_Low",required=True,help="input lower bound on nu in dadi units -- distance between low and high will be searched evenly along a log scale")
parser.add_argument("--nu_High",required=True,help="input upper bound on nu in dadi units -- distance points between nu_Low and nu_High will be searched evenly along a log scale")
parser.add_argument("--T_Low",required=True,help="input lower bound on T in dadi units ; distance between T_Low and T_High will be searched evenly along a log scale")
parser.add_argument("--T_High",required=True,help="input upper bound on T in dadi units ; distance between low and high will be searched evenly along a log scale")
parser.add_argument("--numGridPoints",required=True,help="number of grid points per parameter you want. Keep in mind this can drastically affect run time (e.g. 10 --> 100 calculations; 25 -->  625 calculations")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:

# instead: put input as a range and do the grid search within python'
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --sfs [path to sfs] --Nanc 3500 --nu_Low 1 --nu_High 300 --T_Low 1 --T_High 70 --numGridPoints 25 --outdir [path to outdir]
args = parser.parse_args()
pop=str(args.pop)

outdir=str(args.outdir)
sfs=str(args.sfs)
numGridPoints=float(args.numGridPoints)

### 2019011 update: input is in dadi units: ### don't need to do any rescaling because already input as dadi units
nu_High_rescaled=float(args.nu_High)
nu_Low_rescaled=float(args.nu_Low)
T_High_rescaled=float(args.T_High)
T_Low_rescaled=float(args.T_Low)
maxiter=100

############### Input data ####################################
fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS

# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############ set up grid for specific model  ##########
######## set up model #######
    # this is from CA MSMC model , trimmed at index 22
def sso_model_trim_22_plusContraction_forOptimization((T,nu),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0871326536999987, 1)
	phi = Integration.one_pop(phi, xx, 0.0826612257000005, 1)
	phi = Integration.one_pop(phi, xx, 0.0786283415999992, 1.01729865557432)
	phi = Integration.one_pop(phi, xx, 0.0749695096499999, 1.01729865557432)
	phi = Integration.one_pop(phi, xx, 0.0716360027500006, 1.07723241130487)
	phi = Integration.one_pop(phi, xx, 0.0685862595499991, 1.07723241130487)
	phi = Integration.one_pop(phi, xx, 0.0657873176000006, 1.17798637196801)
	phi = Integration.one_pop(phi, xx, 0.0632047812999999, 1.17798637196801)
	phi = Integration.one_pop(phi, xx, 0.0608214528499997, 1.31205998407016)
	phi = Integration.one_pop(phi, xx, 0.0586072361000001, 1.31205998407016)
	phi = Integration.one_pop(phi, xx, 0.0565520989999997, 1.44976869859196)
	phi = Integration.one_pop(phi, xx, 0.0546331111500001, 1.44976869859196)
	phi = Integration.one_pop(phi, xx, 0.0528431067999999, 1.52089338093277)
	phi = Integration.one_pop(phi, xx, 0.0511648881499993, 1.52414609426596)
	phi = Integration.one_pop(phi, xx, 0.0495898562999998, 1.47942951405878)
	phi = Integration.one_pop(phi, xx, 0.0481094123500004, 1.38075051784768)
	phi = Integration.one_pop(phi, xx, 0.0467149574000002, 1.23285962527743)
	phi = Integration.one_pop(phi, xx, 0.0453978925500004, 1.0486821501222)
	phi = Integration.one_pop(phi, xx, 0.0441550648699991, 0.82602781572228)
	phi = Integration.one_pop(phi, xx, 0.0429771588850002, 0.52406305650732)
	phi = Integration.one_pop(phi, xx, 0.041860735035, 0.191573518835868)
	phi = Integration.one_pop(phi, xx, 0.04080092061, 0.496812482450453)
	### add contraction:
	phi = Integration.one_pop(phi, xx,T,  nu)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs

func = sso_model_trim_22_plusContraction_forOptimization #

param_names= ("nu","T")
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func) # this will give you the exp SFS



from sklearn.model_selection import ParameterGrid
# from numpy use linspace where you tell it how many points you want between two numbers
# note that log10(nu_Low_rescaled) gives you the exponent to put in as the start
# this is correct! so for example if you wanted to search 10 points between
# 0.0001 and 0.001
# you would do np.logspace(-4,-3,10) where -4 and -3 are the exponents
# so to get those exponents you take the base10 log of 0.0001 to get -4 and put that in: *make sure to use log10!!! (not natural log) ***
nus= np.logspace(math.log10(nu_Low_rescaled),math.log10(nu_High_rescaled),numGridPoints) # return values evenly spaced along log scale,from base**start to base**stop base =10, get 10 data points
Ts= np.logspace(math.log10(T_Low_rescaled),math.log10(T_High_rescaled),numGridPoints) # goes from approx 1 gen to ~70 gen
# set up your set of parameters
param_grid = {'nu': nus, 'T' : Ts}
# use the sklearn module to make a list with every pair of parameters in it:
# similar to R expand
grid = ParameterGrid(param_grid)

# nu: Ratio of contemporary to ancient population size
# T: Time in the past at which size change happened (in units of 2*Na
# don't need upper/lower bounds because not optimizing


#############3 set up outfile #########
outputFile=open(str(outdir)+"/"+"dadi.grid.search."+str(pop)+"."+str(modelName)+".LL.output.txt","w")
param_names_str='\t'.join(str(x) for x in param_names)
header=param_names_str+"\ttheta\tLL_model\tLL_data\tmodelFunction\tsampleSize\texpectedSFS_fold_Theta1\tobservedSFS_folded" # add additional parameters theta, log-likelihood, model name, run number and rundate
# set Nanc externally from past runs for each population
outputFile.write(header)
outputFile.write("\n")
# function to get the exp SFS and LL based on parameters YOU provide

def provideParams_2EpochModel(fs,nu,T):
    # get sample size from sfs
    ns= fs.sample_sizes
    # get pts from sample size
    pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]
    # get the expected sfs for this set of parameters:
    model=func_ex([nu,T],ns,pts_l)
    # get the LL of model
    ll_model = dadi.Inference.ll_multinom(model, fs)
    # also get the LL of the data to itself (best possible ll)
    ll_data=dadi.Inference.ll_multinom(fs, fs)
    # Note Nanc is not calculated from the SFS, it is set externally based on past runs (could change that if I wanted )
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)

    model_folded=model.fold()
    # note use of sfs.tolist() has to have () at the end
    # otherwise you get weird newlines in the mix
    output=[nu,T,theta,ll_model,ll_data,func.func_name,ns[0],model_folded.tolist(),fs.tolist()] # put all the output terms together
    output='\t'.join(str(x) for x in output)
    return(output)


# run the function on the grid of all parameter pairs:
for params in grid:
    output=provideParams_2EpochModel(fs=fs, nu=params['nu'],T=params['T'])
    outputFile.write(output)
    outputFile.write("\n")

outputFile.close()

###### exit #######
sys.exit()
