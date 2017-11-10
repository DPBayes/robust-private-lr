#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Generate synthetic data sets and noise samples used in deciding the optimal budget split

import numpy as np
import diffpri as dp

# Program parameters
n = 500		# number of samples in synthetic data sets
d = 10		# dimensionality of synthetic data sets
eps = 2.0	# privacy budget used in tests
l0 = 1		# precision parameters ..
l = 1		# .. (correspond to the means of the gamma hyperpriors)
ln = 5		# average results over ln synthetic data sets and ln noise samples for each data set

datapath = ''	# TODO: set output file path

np.random.seed(42) # for reproducibility

# Generate synthetic data (ln)
xdata = []
ydata = []
sxdata = []
sydata = []
for i in range(ln):
	x = np.random.normal(0.0,1.0,(n,d))
	x = dp.xnormalise(x)
	sx = np.std(x,ddof=1)
	b = np.random.normal(0.0,1.0/np.sqrt(l0),d)
	y = np.random.normal(np.dot(x,b),1.0/np.sqrt(l)).reshape(n,1)
	y = dp.ynormalise(y)
	sy = np.std(y,ddof=1)
	xdata.append(x)
	ydata.append(y)
	sxdata.append(sx)
	sydata.append(sy)

# Generate noise for each data (ln^2)
wnoise = []
lnoise = []
vnoise = []
for jj in range(ln*ln):
	W = dp.symmetric(np.random.laplace(scale=(d**2+d)/eps,size=int((d**2+d)/2)),d)
	L = np.random.laplace(scale=(2.0*d)/eps,size=d).reshape(d,1)
	V = np.random.laplace(scale=1.0/eps,size=1)
	wnoise.append(W)
	lnoise.append(L)
	vnoise.append(V)

# Save into .npz-file
np.savez(datapath+'budgetsplit-data.npz',xdata=xdata,ydata=ydata,sxdata=sxdata,sydata=sydata)
np.savez(datapath+'budgetsplit-noise.npz',wnoise=wnoise,lnoise=lnoise,vnoise=vnoise)
