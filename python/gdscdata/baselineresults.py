#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data

# Compute final results from intermediate results from different drugs and folds

import numpy as np
import csv
import os.path
import sys

yd = 265 # number of drugs
cv = 50 # number of cv folds
nx = 1
ny = 1

# TODO: set:
wpc = False	# set True to compute wpc-index results
corr = True	# set True to compute Spearman's rank correlation coefficient results
inpath = 'baseline-results/' 	# set path for individual files from different drugs and folds
outpath = 'resultsdata/' 		# set path for computed final results
inprefix = 'base-' 			# set input file prefix
outprefix = 'corr-e0-' 			# set output file prefix
method = 'base-fixd-d64'		# set name of used method for output file

indatapath = inpath+inprefix+'corr-fixd-d64-'
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

	res = np.mean(aod,axis=0)		

	# Save data
	np.savetxt(outdatapath+'.csv',res,delimiter=',')

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

	res = np.mean(aod,axis=0)		

	# Save data
	np.savetxt(outdatapath+'.csv',res,delimiter=',')
