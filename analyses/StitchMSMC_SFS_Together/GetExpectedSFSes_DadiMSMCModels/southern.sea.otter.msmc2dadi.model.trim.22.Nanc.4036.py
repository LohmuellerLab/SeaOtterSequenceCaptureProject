# dadi function for sso msmc model; trim point = 22; Nanc=4036
# note: order of variables goes from oldest --> youngest; so nu0 is Nanc, T0 is time for Nanc. nu32 would be most recent if there were 33 time indices.
# This is reverse order for time index numbering in MSMC so don't get confused. They are listed in descending order because that's how they come from MSMC.

from dadi import Numerics, PhiManip, Integration, Spectrum, Plotting, Inference
import numpy as np
import matplotlib.pyplot as pyplot


Nanc = 4036
# params in dadi units (nu = Nx/Nanc; T = gen/(2*Nanc)) in same order that they should go into function (all nus from nu32 --> nu0 followed by Ts from T32 --> T0) -- note that order is from recent --> old because that is how it comes out of msmc. But then dadi function starts with oldest (nu0) first.
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


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_22_function = Numerics.make_extrap_log_func(sso_model_trim_22)


# Version of model with nu and T values plugged in plus one more epoch for dadi to optimize 
def sso_model_trim_22_plusContraction_forOptimization((contractionGen_RescBy2Nanc,contractionSize_RescByNanc),ns,pts):
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
	phi = Integration.one_pop(phi, xx,contractionGen_RescBy2Nanc,  contractionSize_RescByNanc)
	# get expected SFS:
	fs = Spectrum.from_phi(phi,ns,(xx,))
	return fs


# wrap your function in Numerics.make_extrap_log_func to make extrapolation version of your function
extrap_sso_model_trim_22_plusContraction_forOptimization_function = Numerics.make_extrap_log_func(sso_model_trim_22_plusContraction_forOptimization)
