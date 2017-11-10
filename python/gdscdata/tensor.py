#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data

# Run: python3 tensor.py drugid seed
# where 
# - drugid is an integer in [0,1,...,264] (specifies drug)
# - seed is an integer (specifies cv fold)
# This program does 1-fold cv for given drug for all test tensors.
# The cv split is defined by the given random seed.

import sys
import os

# Create a new Theano compile directory on local node (used in cluster computing)
if len(sys.argv) > 2:
	v1 = int(sys.argv[1])
	v2 = int(sys.argv[2])
	mypath1 = "theano"
	mypath2 = mypath1+"/theano-tmp-"+str(v1)+"-"+str(v2)
	if not os.path.exists(mypath1):
		os.mkdir(mypath1)
	if not os.path.exists(mypath2):
		os.mkdir(mypath2)
	os.environ["THEANO_FLAGS"] = "base_compiledir="+mypath1+",compiledir="+mypath2

import diffpri as dp
import numpy as np
import csv

# Import data
datapath = '' # TODO: set path for input and output data files
f = open(datapath+'GeneExpressionReducted.csv','rt')
reader = csv.reader(f,delimiter=',')
x = np.array(list(reader)).astype(float)
f.close()
f = open(datapath+'DrugResponse.csv','rt')
reader = csv.reader(f,delimiter=',')
y = np.array(list(reader)).astype(float)
f.close()

# Arguments
if len(sys.argv) > 2:
	drugid = int(sys.argv[1])
	seed = int(sys.argv[2])
else: # default
	drugid = 123
	seed = 42

# Test cases
pv_size = [0,100,200,300,400,500,600,700,800]
pv_max = max(pv_size)
dim = [5,10,15,20,25,30,35,40] 	# A)
npv_size = [0,5,10,15,20,25,30] # B)
eps = [1.0,1.5,2.0,2.5,3.0]	# C)
n_test = 100
mcmc = True # use priors instead of fixed values for precision parameter lambda,lambda_0

# All tensors
for test in range(3):

	# A) Dimensionality
	if test == 0:

		# Fetch clipping threshold omegas for test A
		f = open(datapath+'A-WX.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WX = np.array(list(reader)).astype(float)
		f.close()
		f = open(datapath+'A-WY.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WY = np.array(list(reader)).astype(float)
		f.close()

		R = np.zeros((len(pv_size),len(dim)),dtype=np.float64)
		for i in range(len(pv_size)):

			n_pv = pv_size[i]
			n_npv = 10
			e = 2.0
			for j in range(len(dim)):
			
				d = dim[j]
				w_x = WX[i,j]
				w_y = WY[i,j]
			
				# Process data
				nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

				# Fit model
				if mcmc:
					pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
				else:
					pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

				# Evaluate
				R[i,j] = dp.precision(pred,y_test)

		# Save results into file
		csvpath = datapath+'cliptest-drugsens-A-'+str(drugid)+'-'+str(seed)+'.csv'
		np.savetxt(csvpath,R,delimiter=',')

	# B) Non-private data size
	if test == 1:

		# Fetch clipping threshold omegas for test B (in A-W*[:,1])
		f = open(datapath+'A-WX.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WX = np.array(list(reader)).astype(float)
		f.close()
		f = open(datapath+'A-WY.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WY = np.array(list(reader)).astype(float)
		f.close()
	
		R = np.zeros((len(pv_size),len(npv_size)),dtype=np.float64)
		for i in range(len(pv_size)):

			n_pv = pv_size[i]
			d = 10
			e = 2.0
			w_x = WX[i,1]
			w_y = WY[i,1]
			for j in range(len(npv_size)):
			
				n_npv = npv_size[j]

				# If no training data
				if n_pv == 0 and n_npv == 0:
					R[i,j] = 0.0
					continue
			
				# Process data
				nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

				# Fit model
				if mcmc:
					pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
				else:
					pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

				# Evaluate
				R[i,j] = dp.precision(pred,y_test)

		# Save results into file
		csvpath = datapath+'cliptest-drugsens-B-'+str(drugid)+'-'+str(seed)+'.csv'
		np.savetxt(csvpath,R,delimiter=',')

	# C) Privacy parameter
	if test == 2:

		# Fetch clipping threshold omegas for test C
		f = open(datapath+'C-WX.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WX = np.array(list(reader)).astype(float)
		f.close()
		f = open(datapath+'C-WY.csv','rt')
		reader = csv.reader(f,delimiter=',')
		WY = np.array(list(reader)).astype(float)
		f.close()

		R = np.zeros((len(pv_size),len(eps)),dtype=np.float64)
		for i in range(len(pv_size)):

			n_pv = pv_size[i]
			n_npv = 10
			d = 10
			for j in range(len(eps)):
		
				e = eps[j]
				w_x = WX[i,j]
				w_y = WY[i,j]
			
				# Process data
				nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)

				# Fit model
				if mcmc:
					pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
				else:
					pred = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

				# Evaluate
				R[i,j] = dp.precision(pred,y_test)

		# Save results into file
		csvpath = datapath+'cliptest-drugsens-C-'+str(drugid)+'-'+str(seed)+'.csv'
		np.savetxt(csvpath,R,delimiter=',')
