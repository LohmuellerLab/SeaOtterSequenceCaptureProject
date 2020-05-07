# dadi function for nso msmc model; trim point = 20; Nanc=4247
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 4247
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
params = (0.904484906437176, 0.663969387506398, 0.571544862138562, 0.619685630445168, 0.79005829297294, 1.02558503749661, 1.259955414543, 1.43652643457535, 1.53009528131691, 1.54650880264499, 1.47126650327662, 1.47126650327662, 1.31334850354829, 1.31334850354829, 1.16537332842806, 1.16537332842806, 1.05963716704008, 1.05963716704008, 1, 1, 0.04439683224, 0.0455501223099999, 0.0467642982499995, 0.0480467154000004, 0.0494006428000006, 0.0508308478, 0.052349589299999, 0.0539609536000003, 0.0556744753999996, 0.0574996894000009, 0.0594488544999996, 0.0615355917000003, 0.0637735219999999, 0.0661803527000002, 0.0687765152999997, 0.0715851654999999, 0.0746321831999993, 0.0779488967000003, 0.0815775310999996, 0.0855562252000041)



def nso_model_trim_20((nu19, nu18, nu17, nu16, nu15, nu14, nu13, nu12, nu11, nu10, nu9, nu8, nu7, nu6, nu5, nu4, nu3, nu2, nu1, nu0, T19, T18, T17, T16, T15, T14, T13, T12, T11, T10, T9, T8, T7, T6, T5, T4, T3, T2, T1, T0),ns,pts):
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
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_20_function = Numerics.make_extrap_log_func(nso_model_trim_20)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def nso_model_trim_20_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
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
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_nso_model_trim_20_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(nso_model_trim_20_plusContraction_forOptimization)
