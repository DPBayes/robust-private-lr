#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# Generate synthetic data

import numpy as np
import diffpri as dp

# Program parameters
n = 1000	# number of samples in synthetic data sets
d = 10		# dimensionality of synthetic data sets
l0 = 1		# precision parameters
l = 1
datapath = ''	# TODO: set output file path
np.random.seed(42) # for reproducibility

# Generate synthetic data
x = np.random.normal(0.0,1.0,(n,d))
b = np.random.normal(0.0,1.0/np.sqrt(l0),d)
y = np.random.normal(np.dot(x,b),1.0/np.sqrt(l)).reshape(n,1)
np.savetxt('y_data.csv',y,delimiter=',')
np.savetxt('x_data.csv',x,delimiter=',')
