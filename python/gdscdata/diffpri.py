#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Tools for executing the robust private linear regression algorithm and studying its performance

import sys, os
import numpy as np
from numpy import linalg as la
from math import erf
from scipy import optimize
from scipy.stats import spearmanr
import theano.tensor as th
from theano import shared
from pymc3 import Model, Normal, Gamma, MvNormal, find_MAP, NUTS, sample, DensityDist, Deterministic
from pymc3.variational.advi import advi, sample_vp
import warnings
import matplotlib.pyplot as plt 
from matplotlib import cm

# Centers and L2-normalises x-data (removes columnwise mean, normalises rows to norm 1)
def xnormalise(x):
	n = x.shape[0]
	d = x.shape[1]
	if n == 0:
		return x
	else:
		z = x-np.dot(np.ones((n,1),dtype=np.float),np.nanmean(x,0).reshape(1,d))
		return np.divide(z,np.dot(np.sqrt(np.nansum(np.power(z,2.0),1)).reshape(n,1),np.ones((1,d),dtype=np.float)))

# Centers y-data (removes columnwise mean, except for columns where all samples have / all but one sample has missing drug response(s))
def ynormalise(y):
	n = y.shape[0]
	d = y.shape[1]
	if n == 0:
		return y
	else:
		with warnings.catch_warnings():
			warnings.simplefilter("ignore", category=RuntimeWarning)
			m = np.nanmean(y,0)
		ind = np.where(np.sum(~np.isnan(y),0)<=1)[0]
		m[ind] = 0.0 # don't center samples of size <= 1
		return y-np.dot(np.ones((n,1),dtype=np.float),m.reshape(1,d))

# Clip data
def clip(x,y,B_x,B_y):
	C = np.multiply(np.sign(x),np.minimum(np.abs(x),B_x))
	with np.errstate(invalid='ignore'):
		D = np.multiply(np.sign(y),np.minimum(np.abs(y),B_y))
	return C,D

# Selects drug based on drugid, removes cell lines with missing drug response
def ignoreNaN(xx,yy,drugid):
	ind = np.where(np.isnan(yy[:,drugid]))
	y = np.delete(yy[:,drugid],ind,axis=0)
	x = np.delete(xx,ind,axis=0)
	return x,y

# Non-private sufficient statistics
def nxx(x):
	return np.dot(x.T,x)
def nxy(x,y):
	return np.dot(x.T,y)
def nyy(y):
	return np.dot(y.T,y)
	
# Return a symmetric matrix formed from a given vector (noise sample)
# Should be: len(v) = (d^2+d)/2
def symmetric(v,d):
	n = int((d**2+d)/2)
	if len(v) != n:
		sys.exit('Error in symmetric: given vector is not of length (d^2+d)/2')
	m = np.zeros((d,d),dtype=np.float64)
	c = 0
	for i in range(d):
		for j in range(d):
			if i == j:	# diagonal
				m[i,j] = v[c]
				c = c+1
			elif i > j:	# lower triangle
				m[i,j] = m[j,i]
			else:		# upper triangle
				m[i,j] = v[c]
				c = c+1
	return m

# Force matrix to be positive definite
def fixposdef(x):
	if x.shape[0] != x.shape[1]:
		sys.exit('Non-square matrix in fixposdef!')
	w,v = la.eig(x)
	u = np.diag(np.abs(w))
	r = np.abs(np.dot(np.dot(v,u),la.inv(v)))
	return r	

# Check if matrix is positive definite
def checkposdef(x):
	if not np.allclose(x,x.T):
		#print('Non-symmetric matrix found in checkposdef!')
		return False
	else:
		try:
			la.cholesky(x)
		except np.linalg.linalg.LinAlgError:
			#print('Non-positive definite matrix found in checkposdef!')
			return False
	return True

# Precision measure: Spearman's rank correlation coefficient
def precision(y_pred,y_real):
	r = spearmanr(y_pred,y_real)[0]
	if np.isnan(r):
		return 0.0
	else:
		return r

# Precision measure: probabilistic concordance index
def pc(pred,real,sd):
	n = real.shape[0]
	pred = -1.0*pred
	real = -1.0*real
	rho = 0.0
	if sd < 0:
		sd = -1.0*sd
	for i in range(n):
		for j in range(n):
			if i<j:
				c = 0.5
				if pred[i] > pred[j]:
					c = 0.5 * (1+erf((real[i]-real[j])/(2*sd)))
				if pred[i] < pred[j]:
					c = 0.5 * (1+erf((real[j]-real[i])/(2*sd)))
				rho = rho + c
	rho = (2.0/(n*(n-1)))*rho
	return rho

# Choose optimal w_x,w_y for clipping thresholds
# The optimal privacy budget split p1,p2,p3 was first studied with budgetsplit-*.py scripts 
# and then hardcoded here. Alternative splits can also be used.
def omega(n,d,eps,mcmc,ln=20,p1=0.35,p2=0.6,p3=0.05):
	
	plotting_on = False
	np.random.seed(42)	# for reproducibility

	# Precision parameters (correspond to the means of the gamma hyperpriors)
	l = 1.0
	l0 = 1.0

	l1 = ln
	l2 = ln
	st = np.arange(0.1,2.1,0.1)
	lenC1 = len(st)
	lenC2 = lenC1
	err = np.zeros((lenC1,lenC2),dtype=np.float64)
	
	for i in range(l1):

		# Create synthetic data
		x = np.random.normal(0.0,1.0,(n,d))
		x = xnormalise(x)
		sx = np.std(x,ddof=1)
		b = np.random.normal(0.0,1.0/np.sqrt(l0),d)
		y = np.random.normal(np.dot(x,b),1.0/np.sqrt(l)).reshape(n,1)
		y = ynormalise(y)
		sy = np.std(y,ddof=1)
		
		# Thresholds to be tested
		cs1 = st*sx
		cs2 = st*sy

		for j in range(l2):

			# Generate noise
			if mcmc:
				U = symmetric(np.random.laplace(scale=(d**2+d)/(p1*eps),size=int((d**2+d)/2)),d)
				V = np.random.laplace(scale=(2.0*d)/(p2*eps),size=d).reshape(d,1)
			else:
				sys.exit('Error in omega: non-mcmc version not currently implemented.')

			for ci1 in range(lenC1):
				c1 = cs1[ci1]
				for ci2 in range(lenC2):
					c2 = cs2[ci2]

					# Clip data
					xc,yc = clip(x,y,c1,c2)

					# Perturbed suff.stats.
					xx = nxx(xc) + U*(c1**2.0)
					xy = nxy(xc,yc) + V*c1*c2

					# Prediction
					prec = l0*np.identity(d) + l*xx
					mean = np.linalg.solve(prec,l*xy)
					pred = np.dot(x,mean)

					# Precision
					rho = precision(pred,y)
					err[ci1,ci2] = err[ci1,ci2] + rho

	# Average
	err = err/float(l1*l2)

	# Choose best
	ind = np.unravel_index(err.argmax(),err.shape)
	w_x = st[ind[0]]
	w_y = st[ind[1]]

	if plotting_on:
		cmap = cm.viridis
		interp = 'spline16'
		cmin = np.min(err)
		cmax = np.max(err)
		ax = plt.subplot(111)
		plt.imshow(err,cmap=cmap,vmin=cmin,vmax=cmax,interpolation=interp,origin='lower',extent=[0,err.shape[1]-1,0,err.shape[0]-1])
		plt.autoscale(False)
		plt.plot(ind[1],ind[0],'rx')
		plt.colorbar()
		plt.xlabel('w_y')
		plt.ylabel('w_x')
		tix = list(st)
		plt.xticks(range(len(st)))
		a = ax.get_xticks().tolist() 
		a = tix
		ax.set_xticklabels(a)
		i = 1
		for label in ax.get_xticklabels():
			if i%2 != 0:
				label.set_visible(False)
			i = i+1
		plt.yticks(range(len(st)))
		a = ax.get_yticks().tolist() 
		a = tix
		ax.set_yticklabels(a)
		i = 1
		for label in ax.get_yticklabels():
			if i%2 != 0:
				label.set_visible(False)
			i = i+1
		plt.title('Average Spearman rank correlation coefficient on a synthetic data set')
		plt.show()

	return w_x,w_y

# Compute prediction using ADVI (a faster alternative to 'doMCMC')
def doADVI(n,xx,xy,yy,x):
	
	# Optional setting for reproducibility
	use_seed = False

	d = xx.shape[0]
	ns = 5000
	if use_seed:	# optional setting for reproducibility
		seed = 42

	# Disable printing
	sys.stdout = open(os.devnull, 'w')

	# Sufficient statistics
	NXX = shared(xx)
	NXY = shared(xy)
	NYY = shared(yy)

	# Define model and perform MCMC sampling
	with Model() as model:

		# Fixed hyperparameters for priors
		b0 = Deterministic('b0',th.zeros((d),dtype='float64'))
		ide = Deterministic('ide',th.eye(d,m=d,k=0,dtype='float64'))

		# Priors for parameters
		l0 = Gamma('l0',alpha=2.0,beta=2.0)
		l = Gamma('l',alpha=2.0,beta=2.0)
		b = MvNormal('b',mu=b0,tau=l0*ide,shape=d)

		# Custom log likelihood
		def logp(xtx,xty,yty):	
			return (n/2.0)*th.log(l/(2*np.pi))+(-l/2.0)*(th.dot(th.dot(b,xtx),b)-2*th.dot(b,xty)+yty)

		# Likelihood
		delta = DensityDist('delta',logp,observed={'xtx':NXX,'xty':NXY,'yty':NYY})

		# Inference
		if use_seed:
			v_params = advi(n=ns,random_seed=seed)
			trace = sample_vp(v_params, draws=ns,random_seed=seed)
		else:
			v_params = advi(n=ns)
			trace = sample_vp(v_params, draws=ns)		
	
	# Enable printing
	sys.stdout = sys.__stdout__

	# Compute prediction over posterior
	return np.mean([np.dot(x,trace['b'][i]) for i in range(ns)],0)

# Compute prediction by generating MCMC samples from the posterior distribution (n=n_train)
def doMCMC(n,nxx,nxy,nyy,x):

	# Optional setting for reproducibility
	use_seed = False

	d = nxx.shape[0]
	ns = 2000
	if use_seed:	# optional setting for reproducibility
		seed = 42

	# Disable printing
	sys.stdout = open(os.devnull, 'w')

	# Sufficient statistics
	NXX = shared(nxx)
	NXY = shared(nxy)
	NYY = shared(nyy)

	# Define model and perform MCMC sampling
	with Model() as model:

		# Fixed hyperparameters for priors
		b0 = Deterministic('b0',th.zeros((d),dtype='float64'))
		ide = Deterministic('ide',th.eye(d,m=d,k=0,dtype='float64'))

		# Priors for parameters
		l0 = Gamma('l0',alpha=2.0,beta=2.0)
		l = Gamma('l',alpha=2.0,beta=2.0)
		b = MvNormal('b',mu=b0,tau=l0*ide,shape=d)

		# Custom log likelihood
		def logp(xtx,xty,yty):	
			return (n/2.0)*th.log(l/(2*np.pi))+(-l/2.0)*(th.dot(th.dot(b,xtx),b)-2*th.dot(b,xty)+yty)

		# Likelihood
		delta = DensityDist('delta',logp,observed={'xtx':NXX,'xty':NXY,'yty':NYY})

		# Inference
		print('doMCMC: start NUTS')
		step = NUTS()
		if use_seed:
			trace = sample(ns,step,progressbar=True,random_seed=seed)
		else:
			trace = sample(ns,step,progressbar=True)
	
	# Enable printing
	sys.stdout = sys.__stdout__

	# Compute prediction over posterior
	return np.mean([np.dot(x,trace['b'][i]) for i in range(ns)],0)

# Fit (private or non-private) model and sample posterior using MCMC/ADVI, return prediction on test data.
# The optimal privacy budget split p1,p2,p3 was first studied with budgetsplit-*.py scripts 
# and then hardcoded here. Alternative splits can also be used.
def predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private,p1=0.35,p2=0.6,p3=0.05):
	
	d = nxx_pv.shape[0]

	# Generate noise
	if private:
		W = symmetric(np.random.laplace(scale=((d**2+d)*(B_x**2))/(p1*e),size=int((d**2+d)/2)),d)
		L = np.random.laplace(scale=(2*d*B_x*B_y)/(p2*e),size=d)
		V = np.random.laplace(scale=(B_y**2)/(p3*e),size=1)
	else:
		W = 0.0
		L = 0.0
		V = 0.0

	# Perturbate sufficient statistics
	nxx = nxx_npv+nxx_pv+W
	nxy = nxy_npv+nxy_pv+L
	nyy = nyy_npv+nyy_pv+V

	# Select sampling method and return prediction
	use_advi = True	
	if use_advi:
		return doADVI(n_train,nxx,nxy,nyy,x_test)
	else:
		if private:
			if not checkposdef(nxx): # NB: NUTS MCMC sampler suffers from non-positive definiteness
				nxx = fixposdef(nxx)
		return doMCMC(n_train,nxx,nxy,nyy,x_test)		

# Fit (private or non-private) model using fixed values for lambdas and an even budget split, return prediction on test data.
def predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private):
	
	l = 1.0
	l0 = 1.0
	d = nxx_pv.shape[0]

	# Generate noise
	if private:
		W = symmetric(np.random.laplace(scale=(2*B_x*B_x*(d*d+d))/e,size=int((d**2+d)/2)),d)
		L = np.random.laplace(scale=(4*d*B_x*B_y)/e,size=d)
	else:
		W = 0.0
		L = 0.0

	# Posterior distribution
	prec = l*(nxx_npv + nxx_pv + W) + l0*np.identity(d)
	mean = np.linalg.solve(prec,l*(nxy_npv + nxy_pv + L))

	# Compute prediction
	return np.dot(x_test,mean)

# Process data for model fitting: splits, dimensionality reduction, 
# normalisation, clipping, dropping missing values.
# Return sufficient statistics, test data, clipping thresholds, 
# number of trainings samples, and private = True if private data size > 0
def processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed,clipdata=True):
	
	n_train = n_pv + n_npv

	# Set rng seed
	np.random.seed(seed)
	np.random.seed(int(np.floor(np.random.rand()*5000)))

	# Test/training split + dimensionality reduction
	ind = np.random.permutation(x.shape[0])
	x_test = x[ind[0:n_test],0:d]
	y_test = y[ind[0:n_test],:]
	x_train = x[ind[n_test:],0:d]
	y_train = y[ind[n_test:],:]

	# Training data: private/non-private split
	x_pv = x_train[0:n_pv,:]
	y_pv = y_train[0:n_pv,:]
	x_npv = x_train[pv_max:pv_max+n_npv,:]
	y_npv = y_train[pv_max:pv_max+n_npv,:]
	
	# Normalise x-data
	x_test = xnormalise(x_test)
	x_npv = xnormalise(x_npv)
	x_pv = xnormalise(x_pv)

	# Normalise y-data
	y_test = ynormalise(y_test)
	y_npv = ynormalise(y_npv)
	y_pv = ynormalise(y_pv)

	# Clip data
	if clipdata:
		n = np.sum(~np.isnan(y_pv[:,drugid])) # number of private data
		if n == 1: # std not possible to compute => no clipping
			B_x = np.max(np.abs(x_pv))
			B_y = np.nanmax(np.abs(y_pv))
			x_pv,y_pv = clip(x_pv,y_pv,B_x,B_y)
			x_npv,y_npv = clip(x_npv,y_npv,B_x,B_y)
		elif n > 1:
			B_x = w_x * np.nanstd(x_pv,ddof=1) 
			B_y = w_y * np.nanstd(y_pv,ddof=1)
			x_pv,y_pv = clip(x_pv,y_pv,B_x,B_y)
			x_npv,y_npv = clip(x_npv,y_npv,B_x,B_y)
		else: # no pv data => no clipping
			B_x = 0.0
			B_y = 0.0
	else: # no clipping
		B_x = 0.0
		B_y = 0.0

	# Select drug and drop cell lines with missing response
	x_pv,y_pv = ignoreNaN(x_pv,y_pv,drugid)
	x_npv,y_npv = ignoreNaN(x_npv,y_npv,drugid)
	x_test,y_test = ignoreNaN(x_test,y_test,drugid)
	n_train = x_pv.shape[0] + x_npv.shape[0] # update n_train
	
	# Compute suff.stats
	nxx_pv = nxx(x_pv)
	nxy_pv = nxy(x_pv,y_pv)
	nyy_pv = nyy(y_pv)
	nxx_npv = nxx(x_npv)
	nxy_npv = nxy(x_npv,y_npv)
	nyy_npv = nyy(y_npv)
	
	if x_pv.shape[0] == 0: # no private data after dropping missing values
		private = False
	else:
		private = True

	return nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private
