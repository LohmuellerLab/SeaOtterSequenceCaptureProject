# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018

@author: annabelbeichman
"""
import matplotlib
matplotlib.use('Agg') # so graphics show up on hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array # don't comment this out
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')


modelName="1D.AL.PSMC.Trim20"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easySFS (mask optional)")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --mu 8.64411385098638e-09 --L 4193488 --sfs [path to sfs] --outdir [path to outdir]
args = parser.parse_args()
runNum=str(args.runNum)
pop=str(args.pop)
mu=float(args.mu)
L=float(args.L)
outdir=str(args.outdir)
sfs=str(args.sfs)
maxiter=100
############### Input data ####################################
fs=dadi.Spectrum.from_file(sfs) # this is folded from easy SFS
# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################
### going to go from nanc =1 (4500~) to nu = 4000/4500 = 0.89, for 1000 gen (1000/(2*4500)) = 0.1 [4500 is weighted Ne average of MSMC curve ]
# what if I fix T1 and do inference?

## adding them to end of parameter vector:
def nso_model_trim_20_plusContraction_forOptimization((nu1,T1),ns,pts):
	xx = Numerics.default_grid(pts)
	# intialize phi with ancestral pop: (nu0 = 1)
	phi = PhiManip.phi_1D(xx,nu=1)
	# stays at nu0=1 for T0 duration of time:
	# followed by a number of time steps, with associated pop changes:
	phi = Integration.one_pop(phi, xx, 0.0855562252000041, 1)
	phi = Integration.one_pop(phi, xx, 0.0815775310999996, 1)
	phi = Integration.one_pop(phi, xx, 0.0779488967000003, 1.05963716704008)
	phi = Integration.one_pop(phi, xx, 0.0746321831999993, 1.05963716704008)
	phi = Integration.one_pop(phi, xx, 0.0715851654999999, 1.16537332842806)
	phi = Integration.one_pop(phi, xx, 0.0687765152999997, 1.16537332842806)
	phi = Integration.one_pop(phi, xx, 0.0661803527000002, 1.31334850354829)
	phi = Integration.one_pop(phi, xx, 0.0637735219999999, 1.31334850354829)
	phi = Integration.one_pop(phi, xx, 0.0615355917000003, 1.47126650327662)
	phi = Integration.one_pop(phi, xx, 0.0594488544999996, 1.47126650327662)
	phi = Integration.one_pop(phi, xx, 0.0574996894000009, 1.54650880264499)
	phi = Integration.one_pop(phi, xx, 0.0556744753999996, 1.53009528131691)
	phi = Integration.one_pop(phi, xx, 0.0539609536000003, 1.43652643457535)
	phi = Integration.one_pop(phi, xx, 0.052349589299999, 1.259955414543)
	phi = Integration.one_pop(phi, xx, 0.0508308478, 1.02558503749661)
	phi = Integration.one_pop(phi, xx, 0.0494006428000006, 0.79005829297294)
	phi = Integration.one_pop(phi, xx, 0.0480467154000004, 0.619685630445168)
	phi = Integration.one_pop(phi, xx, 0.0467642982499995, 0.571544862138562)
	phi = Integration.one_pop(phi, xx, 0.0455501223099999, 0.663969387506398)
	phi = Integration.one_pop(phi, xx, 0.04439683224, 0.904484906437176)
	### add contraction:
	phi = Integration.one_pop(phi, xx,T1,  nu1)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs
param_names=("nu1","T1")

lower_bound = [1e-4,  1e-5]
upper_bound = [10,  0.1]
lower_bound = [1e-4, 1e-5]
p0 = [0.01,0.001] # initial parameters


func=nso_model_trim_20_plusContraction_forOptimization # set the function

############### Carry out optimization (same for any model) ########################
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
# perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,                                  lower_bound=lower_bound)
# optimize:
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######

Nanc=theta / (4*mu*L)
nu1_scaled_dip=popt[0]*Nanc
T1_scaled_gen=popt[1]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nu1_scaled_dip","T1_scaled_gen")
scaled_popt=(Nanc,nu1_scaled_dip,T1_scaled_gen)
############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=param_names_str+"\t"+scaled_param_names_str+"\ttheta\tLL\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[popt_str,scaled_popt_str,theta,ll_model,func.func_name,mu,L,maxiter,runNum,todaysdate,p0,upper_bound,lower_bound] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".expSFS"

model.to_file(outputSFS)

############### Output plot ########################
print('Making plots **************************************************')


import matplotlib.pyplot as plt
fig=plt.figure(1)

outputFigure=str(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)

plt.savefig(outputFigure)



###### exit #######
sys.exit()
