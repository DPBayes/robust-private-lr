#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Compute drug statistics for wpc-index computation
# Run: python3 drugstats.py drugid
# where drugid in 0..265 specifies drug

import numpy as np
from scipy.special import erf
import csv
import sys
import os

def pcindex(y,o,sd):
	n = len(y)
	pc = 0.0
	for i in range(n):
		if not np.isnan(y[i]):
			for j in range(i+1,n):
				ci = 0.5
				if not np.isnan(y[j]):
					if o[i] > o[j]:
						ci = 0.5*(1+erf((y[j]-y[i])/(2.0*sd)))
					if o[i] < o[j]:
						ci = 0.5*(1+erf((y[i]-y[j])/(2.0*sd)))
				pc = pc+ci
	return (2.0/float(n*(n-1)))*pc

if len(sys.argv) > 1:
	drugid = int(sys.argv[1])
else: # default
	drugid = 0
	
niter = 1000
datapath = '' # TODO: set path for input and output data
f = open(datapath+'DrugResponse.csv','rt')
reader = csv.reader(f,delimiter=',')
y = np.array(list(reader)).astype(float)
f.close()

y = y[:,drugid]
y = -1*y # take negative
n = y.shape[0]	# 985 samples
sd = np.nanstd(y,ddof=1)	# standard deviation
pc = np.zeros(niter,dtype=np.float)
for iter in range(niter):
	#print('iteration',iter+1,'/',niter)
	g = np.random.rand(n)
	ig = np.argsort(g)	# ascending
	pc[iter] = pcindex(y,ig,sd)
m = np.mean(pc)
s = np.std(pc,ddof=1)
iy = np.argsort(y)[::-1]	# descending
pcd = pcindex(y,iy,sd)
wd = (pcd-m)/s

# Save into file
stats = np.array([sd,wd])
np.savetxt(datapath+'drugstats-'+str(drugid)+'.csv',stats,delimiter=',')
