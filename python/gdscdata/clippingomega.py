#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# Choose omega parameters for clipping in each test case using auxiliary data

# Run: python3 clippingomega.py test
# where 
# test = 0 for test cases in tensor A (private data size vs. dimensionality)
# test = 2 for test cases in tensor C (private data size vs. privacy parameter)
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
pv_size = [0,100,200,300,400,500,600,700,800]
dim = [5,10,15,20,25,30,35,40]	# A)
eps = [1.0,1.5,2.0,2.5,3.0]	# C)
ny = len(pv_size)
csvpath = '' 	# TODO: set path for output csv files
mcmc = True	# True -> use uneven privacy budget split

if len(sys.argv) > 1:
	test = int(sys.argv[1])
else:
	print('No test specified as a command line argument. Choose 0 or 2.')
	sys.exit()

if test == 0:
	# Test A
	nx = len(dim)
	WX = np.zeros((ny,nx),dtype=np.float)
	WY = np.zeros((ny,nx),dtype=np.float)
	for i in range(len(pv_size)):
		for j in range(len(dim)):
			n_pv = pv_size[i]
			n = n_pv	
			d = dim[j]
			e = 2.0
			if i == 0: # no pv data -> no clipping
				WX[i,j] = 0.0
				WY[i,j] = 0.0
			else:
				w_x,w_y = dp.omega(n,d,e,mcmc)
				WX[i,j] = w_x
				WY[i,j] = w_y
	np.savetxt(csvpath+'A-WX.csv',WX,delimiter=',')
	np.savetxt(csvpath+'A-WY.csv',WY,delimiter=',')

if test == 2:
	# Test C
	nx = len(eps)
	WX = np.zeros((ny,nx),dtype=np.float)
	WY = np.zeros((ny,nx),dtype=np.float)
	for i in range(len(pv_size)):
		for j in range(len(eps)):
			n_pv = pv_size[i]
			n = n_pv	
			d = 10
			e = eps[j]
			if i == 0: # no pv data -> no clipping
				WX[i,j] = 0.0
				WY[i,j] = 0.0
			else:
				w_x,w_y = dp.omega(n,d,e,mcmc)
				WX[i,j] = w_x
				WY[i,j] = w_y
	np.savetxt(csvpath+'C-WX.csv',WX,delimiter=',')
	np.savetxt(csvpath+'C-WY.csv',WY,delimiter=',')
