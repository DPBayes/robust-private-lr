#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Combine the results from budgetsplit_test.py and output the optimal privacy budget split

import numpy as np
import csv
import os.path
import sys
import matplotlib.pyplot as plt 
from matplotlib import cm
from matplotlib import ticker

# Same parameters as in budgetsplit_pre.py:
n = 500		# number of samples in synthetic data sets
d = 10		# dimensionality of synthetic data sets
eps = 2.0	# privacy budget used in tests
l0 = 1		# precision parameters ..
l = 1		# .. (correspond to the means of the gamma hyperpriors)
ln = 5		# average results over ln synthetic data sets and ln noise samples for each data set
p = np.arange(0.05,0.95,0.05)
lenp = len(p)
plotting_on = True

# TODO: set:
inpath = 'bs-results/'
outpath = ''
prefix = 'budgetsplit-result-'
postfix = '.csv'

err = np.zeros((lenp,lenp),dtype=np.float64)

for i in range(ln):
	for j in range(ln):
		for k in range(lenp):
		
			# Open file
			filename = inpath+prefix+str(i)+'-'+str(j)+'-'+str(k)+postfix
			f = open(filename,'rt')
			reader = csv.reader(f,delimiter=',')
			r = np.array(list(reader)).astype(float)
			f.close()
			
			# Sum
			err = err + r

# Average
err = err/float(ln*ln)

# Save results
np.savetxt(outpath+'budgetsplit.csv',err,delimiter=',')

# Choose best
ind = np.unravel_index(err.argmax(),err.shape)
p1 = p[ind[0]]
p2 = p[ind[1]]
p3 = 1.0-p1-p2

print('Optimal privacy budget split (n='+str(n)+', d='+str(d)+', e='+str(eps)+'):')
print('p1 =',p1,' for X\'X')
print('p2 =',p2,' for X\'y')
print('p3 =',p3,' for y\'y')

if plotting_on:
	cmap = cm.plasma
	interp = 'spline16'
	cmin = np.min(err)
	cmax = np.max(err)
	ax = plt.subplot(111)
	plt.imshow(err,cmap=cmap,vmin=cmin,vmax=cmax,interpolation=interp,origin='lower',extent=[0,err.shape[1]-1,0,err.shape[0]-1])
	plt.autoscale(False)
	plt.plot(ind[1],ind[0],'rx',ms=10,mew=2)

	cb = plt.colorbar()
	cb.ax.set_ylabel('Spearman rank correlation coefficient')
	tick_locator = ticker.MaxNLocator(nbins=5)
	cb.locator = tick_locator
	cb.update_ticks()

	plt.xlabel('p2 (for X\'y)')
	plt.ylabel('p1 (for X\'X)')
	tix = list(p)
	plt.xticks(range(lenp))
	a = ax.get_xticks().tolist() 
	a = tix
	ax.set_xticklabels(a,rotation='vertical')
	i = 1
	for label in ax.get_xticklabels():
		if i%2 != 0:
			label.set_visible(False)
		i = i+1
	plt.yticks(range(lenp))
	a = ax.get_yticks().tolist() 
	a = tix
	ax.set_yticklabels(a)
	i = 1
	for label in ax.get_yticklabels():
		if i%2 != 0:
			label.set_visible(False)
		i = i+1
	plt.title('Optimal privacy budget split')

	plt.show()
