#!/bin/env python3
# Differentially private Bayesian linear regression
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Precision measures: 
# - Spearman's rank correlation coefficient
# - probabilistic concordance index (afterwards weighted average over drugs (wpc-index) should be computed)

# Run: python3 test.py drugid
# where
# - drugid is an integer in [0,1,...,264] (specifies drug)
# This program does 50-fold cross-validation for given drug for all test cases.

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
datapath = '' # TODO: set path for input and output data
f = open(datapath+'GeneExpressionReducted.csv','rt')
reader = csv.reader(f,delimiter=',')
x = np.array(list(reader)).astype(float)
f.close()
f = open(datapath+'DrugResponse.csv','rt')
reader = csv.reader(f,delimiter=',')
y = np.array(list(reader)).astype(float)
f.close()

# Arguments
if len(sys.argv) > 1:
        drugid = int(sys.argv[1])
else: # default
        drugid = 123

# Test cases
pv_size = [0,100,200,300,400,500,600,700,800]
pv_max = max(pv_size)
d = 10
n_npv = 10
# TODO: set privacy parameter, corresponding string identifier for output files, and index in clipping omegas
e = 1.0
stre = 'e1'
w_ind = 0 # ind = 0 corresp.to e = 1.0
#e = 2.0
#stre = 'e2'
#w_ind = 2 # ind = 2 corresp.to e = 2.0
n_test = 100
mcmc = True # use priors instead of fixed values for precision parameter lambda,lambda_0

# Fetch clipping threshold omegas
f = open(datapath+'C-WX.csv','rt')
reader = csv.reader(f,delimiter=',')
WX = np.array(list(reader)).astype(float)
f.close()
f = open(datapath+'C-WY.csv','rt')
reader = csv.reader(f,delimiter=',')
WY = np.array(list(reader)).astype(float)
f.close()

# Fetch sd from drug stats for pc-index computation
f = open(datapath+'drugstats.csv','rt')
reader = csv.reader(f,delimiter=',')
stats = np.array(list(reader)).astype(float)
f.close()
sd = stats[0,drugid]

# Cross-validation
n_cv = 50
for seed in range(n_cv):

	S = np.zeros((len(pv_size),1),dtype=np.float64)
	R = np.zeros((len(pv_size),1),dtype=np.float64)

	for i in range(len(pv_size)):

		n_pv = pv_size[i]
		w_x = WX[i,w_ind]
		w_y = WY[i,w_ind]

		# Process data
		nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

		# Fit model
		if mcmc:
			pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
		else:
			pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

		# Evaluate
		S[i,0] = dp.precision(pred,y_test)
		R[i,0] = dp.pc(pred,y_test,sd)	

	# Save results into file
	csvpath = datapath+'cliptest-drugsens-corr-'+stre+'-'+str(drugid)+'-'+str(seed)+'.csv'
	np.savetxt(csvpath,S,delimiter=',')
	csvpath = datapath+'cliptest-drugsens-wpc-'+stre+'-'+str(drugid)+'-'+str(seed)+'.csv'
	np.savetxt(csvpath,R,delimiter=',')
