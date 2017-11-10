#!/bin/env python3
# Differentially private Bayesian linear regression
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# Synthetic data
# Precision measure: MSE

# Run: python3 test.py param
# where
# - param = 0 runs test cases with fixed n and varying eps
# - param = 1 runs test cases with fixed eps and varying n
# This program does 50-fold cross-validation.

import sys
import os

# Create a new Theano compile directory on local node (used in cluster computing)
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

# Import data
datapath = '/scratch/work/niemina7/cliptest/' # TODO: set path for input and output data
f = open(datapath+'x_data.csv','rt')
reader = csv.reader(f,delimiter=',')
x = np.array(list(reader)).astype(float)
f.close()
f = open(datapath+'y_data.csv','rt')
reader = csv.reader(f,delimiter=',')
y = np.array(list(reader)).astype(float)
f.close()

# Arguments
if len(sys.argv) > 1:
        param = int(sys.argv[1])
else: # default
        param = 0

# Test cases
eps = [0.1,0.2,0.5,1.0,2.0,5.0,10.0]
pv_size = [0,100,200,300,400,500,600,700,800]
pv_max = max(pv_size)
d = 10
n_npv = 10
n_test = 100
mcmc = True # use priors instead of fixed values for precision parameter lambda,lambda_0
n_cv = 50
drugid = 0

# Fetch clipping threshold omegas
if param == 0:
	t = 'A'
else:
	t = 'B'
f = open(datapath+t+'-WX.csv','rt')
reader = csv.reader(f,delimiter=',')
WX = np.array(list(reader)).astype(float)
f.close()
f = open(datapath+t+'-WY.csv','rt')
reader = csv.reader(f,delimiter=',')
WY = np.array(list(reader)).astype(float)
f.close()

if param == 0:
	# Cross-validation
	for seed in range(n_cv):

		S = np.zeros(len(eps),dtype=np.float64)

		for i in range(len(eps)):

			e = eps[i]
			n_pv = 500
			w_x = WX[i]
			w_y = WY[i]

			# Process data
			nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

			private = False # modification: rlr

			# Fit model
			if mcmc:
				pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
			else:
				pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

			# Evaluate
			S[i] = dp.precision(pred,y_test)

		# Save results into file
		csvpath = datapath+'synth-rlr-'+t+'-'+str(seed)+'.csv'
		np.savetxt(csvpath,S,delimiter=',')

if param == 1:
	# Cross-validation
	for seed in range(n_cv):

		S = np.zeros(len(pv_size),dtype=np.float64)

		for i in range(len(pv_size)):

			n_pv = pv_size[i]
			e = 2.0
			w_x = WX[i]
			w_y = WY[i]

			# Process data
			nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

			private = False # modification: rlr

			# Fit model
			if mcmc:
				pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
			else:
				pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

			# Evaluate
			S[i] = dp.precision(pred,y_test)

		# Save results into file
		csvpath = datapath+'synth-rlr-'+t+'-'+str(seed)+'.csv'
		np.savetxt(csvpath,S,delimiter=',')
