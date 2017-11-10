#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Compute trade-off tensor results produced by tensor.py

import numpy as np
import csv
import os.path
import sys

yd = 265
cv = 50
pv_size = [0,100,200,300,400,500,600,700,800] 
dim = [5,10,15,20,25,30,35,40]		# A)
npv_size = [0,5,10,15,20,25,30] 	# B)
eps = [1.0,1.5,2.0,2.5,3.0]		# C)
ny = len(pv_size)

# TODO: set:
inpath = 'tensor-results/' 		# set path for individual files from different drugs and folds
outpath = 'resultsdata/' 		# set path for computed final results
inprefix = 'cliptest-drugsens-' 	# set input file prefix
outprefix = 'tensor-'		 	# set output file prefix

indatapath = inpath+inprefix
outdatapath = outpath+outprefix

# Tests
for t in range(3):

	if t == 0:
		test = 'A'
		nx = len(dim)
	if t == 1:
		test = 'B'
		nx = len(npv_size)
	if t == 2:
		test = 'C'
		nx = len(eps)

	datapath = indatapath+test+'-'
	aod_tensors = []

	# Folds
	for j in range(cv):

		sum_dp = np.zeros((ny,nx),dtype=np.float)

		# Drugs
		for i in range(yd):

			# Open file
			filename = datapath+str(i)+'-'+str(j)+'.csv'
			if os.path.isfile(filename):
				f = open(filename)
				reader = csv.reader(f,delimiter=',')
				r = np.array(list(reader)).astype(float)
				f.close()
				if r.shape[0] != ny or r.shape[1] != nx:
					print('Incomplete file: '+filename)
			else:
				sys.exit('Missing file: '+filename)

			# Update sum over drugs
			sum_dp = sum_dp + r

		# Average over drugs
		aod_tensors.append(sum_dp/float(yd))

	# Compute mean, median and std of 'averages over drugs' over cv
	mean_dp = np.mean(aod_tensors,axis=0)
	std_dp = np.std(aod_tensors,axis=0,ddof=1)		

	# Save data
	np.savetxt(outdatapath+test+'-mean.csv',mean_dp,delimiter=',')
	np.savetxt(outdatapath+test+'-std.csv',std_dp,delimiter=',')
