#!/bin/env python3
# Differentially private Bayesian linear regression
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Precision measures: 
# - Spearman's rank correlation coefficient
# - probabilistic concordance index (afterwards weighted average over drugs (wpc-index) should be computed)

# Run: python3 baseline.py drugid
# where
# - drugid is an integer in [0,1,...,264] (specifies drug)
# This program does 50-fold cv for given drug for all test cases.
# selectgenes.py should be run first in order to create the needed data files.

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
n_pv = 0
pv_max = max(pv_size)
n_npv = 10
n_test = 100
e = 0.0
w_x = 0.0
w_y = 0.0

# Fetch sd from drug stats for pc-index computation
f = open(datapath+'drugstats.csv','rt')
reader = csv.reader(f,delimiter=',')
stats = np.array(list(reader)).astype(float)
f.close()
sd = stats[0,drugid]

# Cross-validation
n_cv = 50
for seed in range(n_cv):

	### BASELINE: d=10
	d = 10

	# Process data
	nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)
	private = False

	# Fit model
	pred_mcmc_d10 = pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
	pred_fixd_d10 = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

	# Evaluate & save
	base_mcmc_d10_corr = dp.precision(pred_mcmc_d10,y_test)
	tmp = np.array([base_mcmc_d10_corr])
	np.savetxt(datapath+'base-corr-mcmc-d10'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_mcmc_d10_wpc = dp.pc(pred_mcmc_d10,y_test,sd)
	tmp = np.array([base_mcmc_d10_wpc])
	np.savetxt(datapath+'base-wpc-mcmc-d10'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_fixd_d10_corr = dp.precision(pred_fixd_d10,y_test)
	tmp = np.array([base_fixd_d10_corr])
	np.savetxt(datapath+'base-corr-fixd-d10'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_fixd_d10_wpc = dp.pc(pred_fixd_d10,y_test,sd)
	tmp = np.array([base_fixd_d10_wpc])
	np.savetxt(datapath+'base-wpc-fixd-d10'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	### BASELINE: d=64
	d = 64

	# Process data
	nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,x_test,y_test,B_x,B_y,n_train,private = dp.processData(x,y,d,n_test,n_pv,n_npv,pv_max,w_x,w_y,drugid,seed)
	private = False

	# Fit model
	pred_mcmc_d64 = pred = dp.predictMCMC(n_train,nxx_pv,nxx_npv,nxy_pv,nxy_npv,nyy_pv,nyy_npv,B_x,B_y,e,x_test,private)
	pred_fixd_d64 = dp.predictL(nxx_pv,nxx_npv,nxy_pv,nxy_npv,B_x,B_y,e,x_test,private)

	# Evaluate & save
	base_mcmc_d64_corr = dp.precision(pred_mcmc_d64,y_test)
	tmp = np.array([base_mcmc_d64_corr])
	np.savetxt(datapath+'base-corr-mcmc-d64'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_mcmc_d64_wpc = dp.pc(pred_mcmc_d64,y_test,sd)
	tmp = np.array([base_mcmc_d64_wpc])
	np.savetxt(datapath+'base-wpc-mcmc-d64'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_fixd_d64_corr = dp.precision(pred_fixd_d64,y_test)
	tmp = np.array([base_fixd_d64_corr])
	np.savetxt(datapath+'base-corr-fixd-d64'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')

	base_fixd_d64_wpc = dp.pc(pred_fixd_d64,y_test,sd)
	tmp = np.array([base_fixd_d64_wpc])
	np.savetxt(datapath+'base-wpc-fixd-d64'+'-'+str(drugid)+'-'+str(seed)+'.csv',tmp,delimiter=',')
