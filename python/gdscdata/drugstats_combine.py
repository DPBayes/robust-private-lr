#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Combine drug statistics computed with drugstats.py for wpc-index computation

import numpy as np
import csv

yd = 265 # number of drugs
inpath = 'drugstats-results/' # TODO: set path for input and output data
outpath = ''
prefix = 'drugstats-'
postfix = '.csv'

stats = np.zeros((2,yd),dtype=np.float64)
for i in range(yd):
	filename = inpath+prefix+str(i)+postfix
	f = open(filename,'rt')
	reader = csv.reader(f,delimiter=',')
	s = np.array(list(reader)).astype(float)
	f.close()
	stats[0,i] = s[0]
	stats[1,i] = s[1]
	
np.savetxt(outpath+'drugstats.csv',stats,delimiter=',')
