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

modelName="1D.3Epoch"
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Carry out a grid search for a '+ modelName +' model ')
#parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
#parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
#parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easysfs (mask optional)")
#parser.add_argument("--Nanc",required=True,help="input ancestral size parameter in diploids") # fix this
parser.add_argument("--nuB_Low",required=True,help="input lower bound on nuB (bottleneck size) in dadi units -- distance between low and high will be searched evenly along a log scale")
parser.add_argument("--nuB_High",required=True,help="input upper bound on nuB (bottleneck size) in dadi units -- distance points between nu_Low and nu_High will be searched evenly along a log scale")
parser.add_argument("--nuF_Low",required=True,help="input lower bound on nuF (recovery size) in dadi units -- distance between low and high will be searched evenly along a log scale")
parser.add_argument("--nuF_High",required=True,help="input upper bound on nuF (recovery size) in dadi units -- distance points between nu_Low and nu_High will be searched evenly along a log scale")
parser.add_argument("--TF_Low",required=True,help="input lower bound on TF (duration of recovery period; bneck duration is fixed at 30 gen) in dadi units ; distance between T_Low and T_High will be searched evenly along a log scale")
parser.add_argument("--TF_High",required=True,help="input upper bound on TF (duration of recovery period; bneck duration is fixed at 30 gen) in dadi units; distance between low and high will be searched evenly along a log scale")
parser.add_argument("--numGridPoints",required=True,help="number of grid points per parameter you want. Keep in mind this can drastically affect run time (e.g. 10 --> 100 calculations; 25 -->  625 calculations")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:

# instead: put input as a range and do the grid search within python'
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --sfs [path to sfs] --Nanc 3500 --nu_Low 1 --nu_High 300 --T_Low 1 --T_High 70 --numGridPoints 25 --outdir [path to outdir]
args = parser.parse_args()
#runNum=str(args.runNum)
pop=str(args.pop)
#mu=float(args.mu)
#L=float(args.L)
outdir=str(args.outdir)
sfs=str(args.sfs)
#Nanc=float(args.Nanc) # no longer needed -- calculated from theta
numGridPoints=float(args.numGridPoints)
# not needed any more because you're already in dadi units:
#nu_High_rescaled=float(args.nu_High)/Nanc # rescale by Nanc to get into dadi units
#nu_Low_rescaled=float(args.nu_Low)/Nanc # rescale by Nanc to get into dadi units
#T_High_rescaled=float(args.T_High)/(2*Nanc) # rescale by 2*Nanc to get into dadi units
#T_Low_rescaled=float(args.T_Low)/(2*Nanc) # rescale by 2*Nanc to get into dadi units

### 2019011 update: input is in dadi units: ### don't need to do any rescaling because already input as dadi units
nuB_High_rescaled=float(args.nuB_High)
nuB_Low_rescaled=float(args.nuB_Low)
nuF_High_rescaled=float(args.nuF_High)
nuF_Low_rescaled=float(args.nuF_Low)
TF_High_rescaled=float(args.TF_High)
TF_Low_rescaled=float(args.TF_Low)
maxiter=100

############### Input data ####################################
#for testing
#sfs="/Users/annabelbeichman/Documents/UCLA/Otters/OtterExomeProject/results/datafiles/SFS/20181119/easySFS_projection/neutral/projection-20181221-hetFilter-0.75/dadi-plusMonomorphic/CA-13.plusMonomorphic.sfs"
fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS

# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############ set up grid for specific model  ##########

func = Demographics1D.three_epoch
#     params = (nuB,nuF,TB,TF)
#     ns = (n1,)
#
#     nuB: Ratio of bottleneck population size to ancient pop size
#     nuF: Ratio of contemporary to ancient pop size
#     TB: Length of bottleneck (in units of 2*Na generations)
#     TF: Time since bottleneck recovery (in units of 2*Na generations)
#
#     n1: Number of samples in resulting Spectrum
#     pts: Number of grid points to use in integration.
param_names= ("nuB","nuF","TB","TF")
# want to fix TB as ? ~30 generations? 30/(2*5000) = 0.003
#scaled_param_names=("nu_scaledByGivenNanc_dip"," T_scaledByGivenNanc_gen")
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func) # this will give you the exp SFS


# Nanc is an input parameter
# nu and T are searched through a grid
# but should be relative to Nanc
# https://stackoverflow.com/questions/13370570/elegant-grid-search-in-python-numpy
from sklearn.model_selection import ParameterGrid
# from numpy use linspace where you tell it how many points you want between two numbers
# note that log10(nu_Low_rescaled) gives you the exponent to put in as the start
# this is correct! so for example if you wanted to search 10 points between
# 0.0001 and 0.001
# you would do np.logspace(-4,-3,10) where -4 and -3 are the exponents
# so to get those exponents you take the base10 log of 0.0001 to get -4 and put that in: *make sure to use log10!!! (not natural log) ***
nuBs= np.logspace(math.log10(nuB_Low_rescaled),math.log10(nuB_High_rescaled),numGridPoints) # return values evenly spaced along log scale,from base**start to base**stop base =10, get 10 data points
TB_fixed=0.003 # fixed
nuFs= np.logspace(math.log10(nuF_Low_rescaled),math.log10(nuF_High_rescaled),numGridPoints) # return values evenly spaced along log scale,from base**start to base**stop base =10, get 10 data points
TFs= np.logspace(math.log10(TF_Low_rescaled),math.log10(TF_High_rescaled),numGridPoints) # goes from approx 1 gen to ~70 gen
# set up your set of parameters
param_grid = {'nuB': nuBs, 'nuF' : nuFs, 'TF' : TFs,}
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

def provideParams_3EpochModel(fs,nuB,nuF,TB,TF):
    # get sample size from sfs
    ns= fs.sample_sizes
    # get pts from sample size
    pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]
    # get the expected sfs for this set of parameters:
    # parameters in order: (nuB,nuF,TB,TF)
    model=func_ex([nuB,nuF,TB,TF],ns,pts_l)
    # get the LL of model
    ll_model = dadi.Inference.ll_multinom(model, fs)
    # also get the LL of the data to itself (best possible ll)
    ll_data=dadi.Inference.ll_multinom(fs, fs)
    # Note Nanc is not calculated from the SFS, it is set externally based on past runs (could change that if I wanted )
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)
    #Nanc=theta / (4*mu*L)
    # Set Nanc without theta (maybe? see how this works)
    #nu_scaled_dip=nu*Nanc
    #T_scaled_gen=T*2*Nanc
    # fold exp sfs?
    model_folded=model.fold()
    # note use of sfs.tolist() has to have () at the end
    # otherwise you get weird newlines in the mix
    output=[nuB,nuF,TB,TF,theta,ll_model,ll_data,func.func_name,ns[0],model_folded.tolist(),fs.tolist()] # put all the output terms together
    output='\t'.join(str(x) for x in output)
    return(output)


# run the function on the grid of all parameter pairs:
for params in grid:
    #print(params)
    output=provideParams_3EpochModel(fs=fs,nuB=params['nuB'],nuF=params['nuF'],TB=TB_fixed,TF=params['TF'])
    outputFile.write(output)
    outputFile.write("\n")

outputFile.close()

###### exit #######
sys.exit()
