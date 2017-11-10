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

yd = 265 # number of drugs
cv = 50 # number of cv folds
pv_size = [0,100,200,300,400,500,600,700,800]
nx = len(pv_size)
ny = 1

# TODO: set:
wpc = False	# set True to compute wpc-index results
corr = True	# set True to compute Spearman's rank correlation coefficient results
inpath = 'ts-e2-mcmc-results/' 		# set path for individual files from different drugs and folds
outpath = 'resultsdata/' 		# set path for computed final results
inprefix = 'cliptest-drugsens-' 			# set input file prefix
outprefix = 'corr-e2-' 			# set output file prefix
method = 'mcmc'		# set name of used method for output file

indatapath = inpath+inprefix+'corr-e2-'
outdatapath = outpath+outprefix+method

if wpc:
	# Fetch drug stats for wpc computation
	f = open('drugstats.csv','rt')
	reader = csv.reader(f,delimiter=',')
	stats = np.array(list(reader)).astype(float)
	f.close()

	aod = []

	# Folds
	for j in range(cv):

		prec_num = np.zeros((nx,ny),dtype=np.float)
		prec_den = np.zeros((nx,ny),dtype=np.float)

		# Drugs
		for i in range(yd):
		
			# Drug weight
			wd = stats[1,i]
			if wd < 0:
				wd = -1.0*wd

			# Open file
			filename = indatapath+str(i)+'-'+str(j)+'.csv'
			if os.path.isfile(filename):
				f = open(filename)
				reader = csv.reader(f,delimiter=',')
				r = np.array(list(reader)).astype(float)
				f.close()
					
				if r.shape[0] != nx or r.shape[1] != ny:
					print('Incomplete file: '+filename)
			else:
				sys.exit('Missing file: '+filename)

			prec_num = prec_num + wd*r
			prec_den = prec_den + wd
				
		# Average over drugs
		aod.append(prec_num/prec_den)

	# Compute mean and std of 'averages over drugs' over cv
	mean_dp = np.mean(aod,axis=0)
	std_dp = np.std(aod,axis=0,ddof=1)		

	# Save data
	np.savetxt(outdatapath+'-mean.csv',mean_dp.T,delimiter=',')
	np.savetxt(outdatapath+'-std.csv',std_dp.T,delimiter=',')

if corr:
	aod = []

	# Folds
	for j in range(cv):

		sum_dp = np.zeros((nx,ny),dtype=np.float)

		# Drugs
		for i in range(yd):

			# Open file
			filename = indatapath+str(i)+'-'+str(j)+'.csv'
			if os.path.isfile(filename):
				f = open(filename)
				reader = csv.reader(f,delimiter=',')
				r = np.array(list(reader)).astype(float)
				f.close()
					
				if r.shape[0] != nx or r.shape[1] != ny:
					print('Incomplete file: '+filename)
			else:
				sys.exit('Missing file: '+filename)

			sum_dp = sum_dp + r
				
		# Average over drugs
		aod.append(sum_dp/float(yd))

	# Compute mean and std of 'averages over drugs' over cv
	mean_dp = np.mean(aod,axis=0)
	std_dp = np.std(aod,axis=0,ddof=1)		

	# Save data
	np.savetxt(outdatapath+'-mean.csv',mean_dp.T,delimiter=',')
	np.savetxt(outdatapath+'-std.csv',std_dp.T,delimiter=',')
