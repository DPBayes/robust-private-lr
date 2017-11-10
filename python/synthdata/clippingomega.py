#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# Choose omega parameters for clipping in each test case using auxiliary data

# Run: python3 clippingomega.py test
# where 
# test = 0 for test cases with fixed n and varying eps
# test = 1 for test cases with fixed eps and varying n
# Here it is assumed the optimal privacy budget split is already found and defined in diffpri.py.

import sys
import os

# Create a new Theano compile directory on local node
# (used when running tests on computer cluster)
if len(sys.argv) > 1:
	v1 = int(sys.argv[1])
	mypath1 = "theano"
	mypath2 = mypath1+"/theano-tmp-"+str(v1)
	if not os.path.exists(mypath1):
		os.mkdir(mypath1)
	if not os.path.exists(mypath2):
		os.mkdir(mypath2)
	os.environ["THEANO_FLAGS"] = "base_compiledir="+mypath1+",compiledir="+mypath2

import diffpri as dp
import numpy as np
import csv

# Test cases
eps = [0.1,0.2,0.5,1.0,2.0,5.0,10.0]
pv_size = [0,100,200,300,400,500,600,700,800]
d = 10
csvpath = '/scratch/work/niemina7/cliptest/' 	# TODO: set path for output csv files
mcmc = True	# True -> use uneven privacy budget split

if len(sys.argv) > 1:
	test = int(sys.argv[1])
else:
	print('No test specified as a command line argument. Choose 0 or 1.')
	sys.exit()

if test == 0:
	nx = len(eps)
	WX = np.zeros(nx,dtype=np.float)
	WY = np.zeros(nx,dtype=np.float)
	for i in range(len(eps)):
		e = eps[i]
		n_pv = 500
		n = n_pv	
		w_x,w_y = dp.omega(n,d,e,mcmc)
		WX[i] = w_x
		WY[i] = w_y
	np.savetxt(csvpath+'A-WX.csv',WX,delimiter=',')
	np.savetxt(csvpath+'A-WY.csv',WY,delimiter=',')

if test == 1:
	nx = len(pv_size)
	WX = np.zeros(nx,dtype=np.float)
	WY = np.zeros(nx,dtype=np.float)
	for i in range(len(pv_size)):
		n_pv = pv_size[i]
		n = n_pv	
		e = 2.0
		if i == 0: # no pv data -> no clipping
			WX[i] = 0.0
			WY[i] = 0.0
		else:
			w_x,w_y = dp.omega(n,d,e,mcmc)
			WX[i] = w_x
			WY[i] = w_y
	np.savetxt(csvpath+'B-WX.csv',WX,delimiter=',')
	np.savetxt(csvpath+'B-WY.csv',WY,delimiter=',')
