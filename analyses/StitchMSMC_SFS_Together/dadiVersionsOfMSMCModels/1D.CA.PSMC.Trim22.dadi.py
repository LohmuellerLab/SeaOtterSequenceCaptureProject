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


modelName="1D.CA.PSMC.Trim22"

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
param_names=("nu","T")



upper_bound = [10, 2]
lower_bound = [1e-4, 1e-4]
p0 = [0.01,0.1] # initial parameters


func=sso_model_trim_22_plusContraction_forOptimization # set the function

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
nu_scaled_dip=popt[0]*Nanc
T_scaled_gen=popt[1]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nu_scaled_dip","T_scaled_gen")
scaled_popt=(Nanc,nu_scaled_dip,T_scaled_gen)

############ what if I force Nanc to be what it is in psmc? how does that affect fit? Could get theta from N*mu*L*4 and then multiply up to get counts. Or could divide to get prportional.
############### Write out output (same for any model)
 ########################
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
