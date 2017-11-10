#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Compute the final results from intermediate results from different drugs and folds
# as produced by test.py (or a modified version test_lr.py, test_rlr.py)

import numpy as np
import csv
import os.path
import sys

cv = 50 # number of cv folds
eps = [0.1,0.2,0.5,1.0,2.0,5.0,10.0]
pv_size = [0,100,200,300,400,500,600,700,800]

# TODO: set:
test = 0				# SET: 0 for varying eps, 1 for varying n_pv
method = 'rlr'				# SET: method: rplr, lr, rlr
inpath = 'rlr-results/' 		# SET: path for individual files from different drugs and folds
outpath = 'resultsdata/' 		# path for computed final results
inprefix = 'synth-'+method+'-' 		# input file prefix
outprefix = 'synth-'+method+'-'		# output file prefix

if test == 0:
	t = 'A'
	nx = len(eps)
else:
	t = 'B'
	nx = len(pv_size)

indatapath = inpath+inprefix+t+'-'
outdatapath = outpath+outprefix+t

aod = []

# Folds
for j in range(cv):

	sum_dp = np.zeros(nx,dtype=np.float)

	# Open file
	filename = indatapath+str(j)+'.csv'
	if os.path.isfile(filename):
		f = open(filename)
		reader = csv.reader(f,delimiter=',')
		r = np.array(list(reader)).astype(float)
		f.close()
				
		if r.shape[0] != nx:
			print('Incomplete file: '+filename)
	else:
		sys.exit('Missing file: '+filename)
		
	# Average over drugs
	aod.append(r)

# Compute mean and std of results over cv
mean_dp = np.mean(aod,axis=0)
std_dp = np.std(aod,axis=0,ddof=1)		

# Save data
np.savetxt(outdatapath+'-mean.csv',mean_dp.T,delimiter=',')
np.savetxt(outdatapath+'-std.csv',std_dp.T,delimiter=',')
