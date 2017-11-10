#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Test each sensible split with given synthetic data and noise sample
# that were generated using budgetsplit_pre.py (ln data sets and ln^2 noise samples)

# Run (on a computer cluster): python3 budgetsplit_test.py i j k
# where
# - i = 0..ln-1 (specifies data)
# - j = 0..ln-1 (together with i specifies noise)
# - k = 0..17 (specifies budget split for nxx)
# Should be run for each tuple (i,j,k) in {0..ln-1}x{0..ln-1}x{0..17}
# Here ln = 5

import sys
import os

# Create a new Theano compile directory on local node (used in cluster computing)
if len(sys.argv) > 3:
	v1 = int(sys.argv[1])
	v2 = int(sys.argv[2])
	v3 = int(sys.argv[3])
	mypath1 = "theano"
	mypath2 = mypath1+"/theano-tmp-"+str(v1)+"-"+str(v2)+"-"+str(v3)
	if not os.path.exists(mypath1):
		os.mkdir(mypath1)
	if not os.path.exists(mypath2):
		os.mkdir(mypath2)
	os.environ["THEANO_FLAGS"] = "base_compiledir="+mypath1+",compiledir="+mypath2

import diffpri as dp
import numpy as np
import csv

# Same parameters as in budgetsplit_pre.py:
n = 500		# number of samples in synthetic data sets
d = 10		# dimensionality of synthetic data sets
eps = 2.0	# privacy budget used in tests
l0 = 1		# precision parameters ..
l = 1		# .. (correspond to the means of the gamma hyperpriors)
ln = 5		# average results over ln synthetic data sets and ln noise samples for each data set

datapath = ''	# TODO: set input/output file path

np.random.seed(42) # for reproducibility

# Arguments
if len(sys.argv) > 3:
	i = int(sys.argv[1])
	j = int(sys.argv[2])
	k = int(sys.argv[3])
else: # default
	i = 0
	j = 0
	k = 0

# Import ith data
f = np.load(datapath+'budgetsplit-data.npz')
x = f['xdata'][i]
y = f['ydata'][i]
sx = f['sxdata'][i]
sy = f['sydata'][i]
f.close()

# Import jth noise for ith data
jj = ln*i+j
f = np.load(datapath+'budgetsplit-noise.npz')
W = f['wnoise'][jj]
L = f['lnoise'][jj]
V = f['vnoise'][jj]
f.close()

# Define and test splits
p = np.arange(0.05,0.95,0.05)
lenp = len(p)
minp = 0.03
maxp = 0.97

err = np.zeros((lenp,lenp),dtype=np.float64)

for m in range(lenp):

	# Budget split
	p1 = p[k]
	p2 = p[m]
	p3 = 1.0-p1-p2
	
	# Check if split is sensible
	t1 = minp <= p1 and p1 <= maxp
	t2 = minp <= p2 and p2 <= maxp
	t3 = minp <= p3 and p3 <= maxp

	if not all([t1,t2,t3]):
		continue

	# Clipping omega and thresholds
	w_x,w_y = dp.omega(n,d,eps,True,ln=ln,p1=p1,p2=p2,p3=p3)
	c1 = sx * w_x
	c2 = sy * w_y

	# Clip data
	xc,yc = dp.clip(x,y,c1,c2)

	# Perturbed suff.stats.
	xx = dp.nxx(xc) + W*(c1**2.0)/p1
	xy = dp.nxy(xc,yc) + L*c1*c2/p2
	yy = dp.nyy(yc) + V*(c2**2.0)/p3

	# Prediction
	pred = dp.doADVI(n,xx,xy,yy,x)

	# Precision
	err[k,m] = dp.precision(pred,y)

# Save result
csvpath = datapath+'budgetsplit-result-'+str(i)+'-'+str(j)+'-'+str(k)+'.csv'
np.savetxt(csvpath,err,delimiter=',')
