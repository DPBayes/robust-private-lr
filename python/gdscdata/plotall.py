#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data

# Plot and save results: tensors, correlation (eps=1,eps=2), wpc-index (eps=2), optimal budget split

import numpy as np
import csv
import os.path
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

### GENERAL
outpath = 'resultsdata/'
postmean = '-mean.csv'
poststd = '-std.csv'
yd = 265
cv = 50
pv_size = [0,100,200,300,400,500,600,700,800]
dim = [5,10,15,20,25,30,35,40]		# TENSOR A
npv_size = [0,5,10,15,20,25,30] 	# TENSOR B
eps = [1.0,1.5,2.0,2.5,3.0]		# TENSOR C
n_npv = 10
y_lower_corr = -0.02
y_lower_corr_e1 = -0.05
y_upper_corr = 0.295
y_lower_wpc = 0.495
y_upper_wpc = 0.568
showtitle = False
cols = 2	# figure width: 1 or 2 columns

if cols == 1: # figure width is 1 column
	figpath = 'figures/'
	postfix = '_col1'
	matplotlib.rcParams['axes.linewidth'] = 0.2 
	figsize = (3.1,2.6)	# figure size: comparison plot
	figsizet = (3.1,2.2)	# figure size: tensor
	figsizebs = (3.1,3.0)	# figure size: budget split
	dpi = 1200
	ms = 4		# marker size
	mew = 0		# marker edge width
	mec = 'w'	# marker edge color
	msbs = 7	# marker size: best budgetsplit
	mewbs = 2	# marker edge width: best budgetsplit
	lw1 = 1		# line weight: baselines
	lw2 = 1		# line weight: methods
	lw3 = 1		# line weight: budget splits
	lw = 1		# line weight: error bars
	cs = 2		# size of error bar caps
	fst = 7		# font size: tensor labels and titles
	fs = 7 		# font size: axis labels
	lgdfs = 7	# font size: legend
	lgdncol = 2	# number of columns in legend
	lgdw = 0.2	# legend border width
	ncol = 1	# number of columns in legend
	hl = 3.0	# handlelength in legend labels
	lgdloc = 'upper center' # legend location
	bbta = (0.5,-0.2)	# bbox_to_anchor value for legend
	dashesnpv = (5,3)	# npv methods linestyle
	xlim = 30		# extra space in x-axis
	wspace = 20.0		# tensor: wspace between subfigs
	hspace = 0.6		# tensor: hspace between subfigs
	width_ratios = [1,1,1,1] # tensor: width ratios
	height_ratios = [1,1] 	# tensor: height ratios
	top = 0.8		# tensor: height to width (=1.0) ratio
	lb = 6			# tensor: label font size for tensors b,c
	tytitle = 'Size of\nprivate dataset' # tensor: y-axis label in b,c
	x = np.array(pv_size)
	xlabels = [str(n_npv)+'\n+'+s for s in list(map(str,pv_size))]
elif cols == 2: # arXiv version
	figpath = 'figures/'
	postfix = '_col2'
	matplotlib.rcParams['axes.linewidth'] = 0.2 
	figsize = (6.7,4.7)	# figure size: comparison plot
	figsizet = (6.0,3.5)	# figure size: tensor
	figsizebs = (6.7,5.0)	# figure size: budget split
	dpi = 1200
	ms = 5		# marker size
	mew = 0		# marker edge width
	mec = 'w'	# marker edge color
	msbs = 8	# marker size: best budgetsplit
	mewbs = 2	# marker edge width: best budgetsplit
	lw1 = 1		# line weight: baselines
	lw2 = 1		# line weight: methods
	lw3 = 1		# line weight: budget splits
	lw = 1		# line weight: error bars
	cs = 3		# size of error bar caps
	fst = 7		# font size: tensor labels and titles
	fs = 7 		# font size: axis labels
	lgdfs = 7	# font size: legend
	lgdncol = 2	# number of columns in legend
	lgdw = 0.2	# legend border width
	ncol = 1	# number of columns in legend
	hl = 3.7	# handlelength in legend labels
	lgdloc = 'upper center' # legend location
	bbta = (0.5,-0.1)	# bbox_to_anchor value for legend
	dashesnpv = (6,4)	# npv methods linestyle
	xlim = 30		# extra space in x-axis
	wspace = 1.0		# tensor: wspace between subfigs
	hspace = 0.3		# tensor: hspace between subfigs
	width_ratios = [1.5,1.5,1,1] # tensor: width ratios
	height_ratios = [1,1]	# tensor: height ratios
	top = 1.0		# tensor: height to width (=1.0) ratio
	tytitle = 'Size of private dataset' # tensor: y-axis label in b,c
	x = np.array(pv_size)
	xlabels = [str(n_npv)+'+'+s for s in list(map(str,pv_size))]
else: # another option
	print('todo')

### TENSORS
matplotlib.rcParams.update({'font.size': fst})
relative = True
if relative:
	print('\nPlotting tensors (relative improvement over baseline)...')
else:
	print('\nPlotting tensors (performance)...')
pre = 'tensor-'
f = open(outpath+pre+'A'+postmean)
reader = csv.reader(f,delimiter=',')
A = np.array(list(reader)).astype(float)
f.close()
f = open(outpath+pre+'B'+postmean)
reader = csv.reader(f,delimiter=',')
B = np.array(list(reader)).astype(float)
f.close()
f = open(outpath+pre+'C'+postmean)
reader = csv.reader(f,delimiter=',')
C = np.array(list(reader)).astype(float)
f.close()

if relative:
	base = np.max(A[0,:])
	print('Tensors baseline: d='+str(dim[np.argmax(A[0,:])])+', n_pv=0, n_npv=10, eps=2.0')
	A = A/base
	B = B/base
	C = C/base
cmin = 0
m = max(np.max(A),np.max(B),np.max(C))
if relative:
	from math import trunc
	cmax = trunc(np.floor(m*10))/10
else:
	cmax = m

cmap = cm.jet
interp = 'spline16'
fig = plt.figure(figsize=figsizet,dpi=dpi)
gs = gridspec.GridSpec(2,4,width_ratios=width_ratios,height_ratios=height_ratios,wspace=wspace,hspace=hspace,left=0.0,right=1.0,bottom=0.0,top=top)
if showtitle:
	if relative:
		fig.suptitle('Relative improvement over baseline')
	else:
		fig.suptitle('Spearman rank correlation coefficient')

# TENSOR A
ax1 = plt.subplot(gs[:,:-2])
plt.imshow(A,cmap=cmap,vmin=cmin,vmax=cmax,interpolation=interp,origin='lower',extent=[0,A.shape[1]-1,0,A.shape[0]-1],aspect='auto')
plt.xlabel('a) Reduced dimensionality')
plt.ylabel('Size of private dataset')
##
ax1.xaxis.set_ticks(np.arange(0,len(dim),1))
ax1.set_xticklabels(dim)
ax1.set_yticklabels(pv_size)
if cols == 1:
	ax1.tick_params(axis=u'both', which=u'both',length=0,labelsize=lb)
else:
	ax1.tick_params(axis=u'both', which=u'both',length=0)

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="8%", pad="3%")
if relative:
	cbar = plt.colorbar(ticks=[cmin,1,cmax],ax=ax1,cax=cax1)
	cbar.set_label('Relative improvement on rank correlation',rotation=90)
else:
	cbar = plt.colorbar(ax=ax1,cax=cax1)
	cbar.set_label('Spearman\'s rank correlation coefficient',rotation=90)

# TENSOR B
ax2 = plt.subplot(gs[0,2:])
plt.imshow(B,cmap=cmap,vmin=cmin,vmax=cmax,interpolation=interp,origin='lower',extent=[0,B.shape[1]-1,0,B.shape[0]-1],aspect='auto')
plt.xlabel('b) Size of non-private data')
plt.ylabel(tytitle)
ax2.xaxis.set_ticks(np.arange(0,len(npv_size),1))
ax2.set_xticklabels(npv_size)
ax2.yaxis.set_ticks(np.arange(0,len(pv_size),1))
ax2.set_yticklabels(pv_size)
if cols == 1:
	ax2.tick_params(axis=u'both', which=u'both',length=0,labelsize=lb)
else:
	ax2.tick_params(axis=u'both', which=u'both',length=0)

# TENSOR C
ax3 = plt.subplot(gs[1,2:])
plt.imshow(C,cmap=cmap,vmin=cmin,vmax=cmax,interpolation=interp,origin='lower',extent=[0,C.shape[1]-1,0,C.shape[0]-1],aspect='auto')
plt.xlabel('c) Privacy parameter')
plt.ylabel(tytitle)
ax3.xaxis.set_ticks(np.arange(0,len(eps),1))
ax3.set_xticklabels(eps)
ax3.set_xticklabels(['1','1.5','2','2.5','3'])
ax3.yaxis.set_ticks(np.arange(0,len(pv_size),1))
ax3.set_yticklabels(pv_size)
if cols == 1:
	ax3.tick_params(axis=u'both', which=u'both',length=0,labelsize=lb)
else:
	ax3.tick_params(axis=u'both', which=u'both',length=0)

plt.savefig(figpath+'laplace_mcmc_tensor_corr'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')


### COMPARISON PLOTS GENERAL
matplotlib.rcParams.update({'font.size': fs})
methods = ['base-fixd-d10','base-fixd-d64','mcmc-lr','plr','mcmc','oplr','fmlr','mcmc-rlr'] # order: baselines, non-private methods, private methods, rlr last
methodsnames = ['LR non-private data (d=10)','LR non-private data (d=64)','Linear Regression (LR)','Private LR','Robust private LR','Output perturbed LR','Functional mechanism LR','Robust LR']
numbase = 2 # number of baselines
numnpv = 3 # number of non-private methods (including baselines, excluding rlr)
color = ['k','gray','magenta','b','r','lime','cyan','darkcyan']

### CORRELATION EPS = 1
print('\nPlotting performance comparison (Spearman\'s rank correlation, eps=1)...')
pre_npv = 'corr-e0-'
pre = 'corr-e1-'
meandata = []
stddata = []
for k in range(len(methods)-1):
	method = methods[k]
	if k < numbase:
		f = open(outpath+pre_npv+method+'.csv')
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][0]*np.ones_like(x))
		f.close()
		stddata.append(np.zeros_like(x))
	else:
		if k < numnpv:
			f = open(outpath+pre_npv+method+postmean)
		else:
			f = open(outpath+pre+method+postmean)
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()
		if k < numnpv:
			f = open(outpath+pre_npv+method+poststd)
		else:
			f = open(outpath+pre+method+poststd)
		reader = csv.reader(f,delimiter=',')
		stddata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()

fig = plt.figure(figsize=figsize,dpi=dpi)
ax = fig.add_subplot(111)
ax.set_xticks(pv_size)
ax.set_xticklabels(xlabels)

for k,c in zip(range(len(methods)-1),color):
	if k < numbase:
		ax.plot(x,meandata[k],'--',ms=ms,lw=lw1,label=methodsnames[k],color=c,dashes=dashesnpv)
	else:
		if k < numnpv:
			ax.plot(x,meandata[k],'s--',ms=ms,lw=lw1,label=methodsnames[k],color=c,mew=mew,mec=mec,dashes=dashesnpv)
		else:		
			ax.plot(x,meandata[k],'s-',ms=ms,lw=lw2,label=methodsnames[k],color=c,mew=mew,mec=mec)
		ax.errorbar(x,meandata[k],lw=lw,yerr=[stddata[k],stddata[k]],fmt='none',ecolor=c,capthick=lw,capsize=cs)
lgd = ax.legend(loc=lgdloc,numpoints=1,prop={'size':lgdfs},ncol=ncol,handlelength=hl,bbox_to_anchor=bbta)
lgd.get_frame().set_linewidth(lgdw)
lgd.get_frame().set_edgecolor("black")
plt.xlabel('Size of dataset (internal+external)')
plt.ylabel('Spearman\'s rank correlation coefficient')
if showtitle:
	plt.title('Comparison of different methods ('+r'$\epsilon$'+'=1.0)')
plt.xlim([min(pv_size)-xlim,max(pv_size)+xlim])
plt.ylim([y_lower_corr_e1,y_upper_corr])

plt.savefig(figpath+'laplace_mcmc_corr1'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')

### CORRELATION EPS = 2
print('\nPlotting performance comparison (Spearman\'s rank correlation, eps=2)...')
pre_npv = 'corr-e0-'
pre = 'corr-e2-'
meandata = []
stddata = []
for k in range(len(methods)-1):
	method = methods[k]
	if k < numbase:
		f = open(outpath+pre_npv+method+'.csv')
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][0]*np.ones_like(x))
		f.close()
		stddata.append(np.zeros_like(x))
	else:
		if k < numnpv:
			f = open(outpath+pre_npv+method+postmean)
		else:
			f = open(outpath+pre+method+postmean)
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()
		if k < numnpv:
			f = open(outpath+pre_npv+method+poststd)
		else:
			f = open(outpath+pre+method+poststd)
		reader = csv.reader(f,delimiter=',')
		stddata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()

fig = plt.figure(figsize=figsize,dpi=dpi)
ax = fig.add_subplot(111)
ax.set_xticks(pv_size)
ax.set_xticklabels(xlabels)

for k,c in zip(range(len(methods)-1),color):
	if k < numbase:
		ax.plot(x,meandata[k],'--',ms=ms,lw=lw1,label=methodsnames[k],color=c,dashes=dashesnpv)
	else:
		if k < numnpv:
			ax.plot(x,meandata[k],'s--',ms=ms,lw=lw1,label=methodsnames[k],color=c,mew=mew,mec=mec,dashes=dashesnpv)
		else:		
			ax.plot(x,meandata[k],'s-',ms=ms,lw=lw2,label=methodsnames[k],color=c,mew=mew,mec=mec)
		ax.errorbar(x,meandata[k],lw=lw,yerr=[stddata[k],stddata[k]],fmt='none',ecolor=c,capthick=lw,capsize=cs)
lgd = ax.legend(loc=lgdloc,numpoints=1,prop={'size':lgdfs},ncol=ncol,handlelength=hl,bbox_to_anchor=bbta)
lgd.get_frame().set_linewidth(lgdw)
lgd.get_frame().set_edgecolor("black")
plt.xlabel('Size of dataset (internal+external)')
plt.ylabel('Spearman\'s rank correlation coefficient')
if showtitle:
	plt.title('Comparison of different methods ('+r'$\epsilon$'+'=2.0)')
plt.xlim([min(pv_size)-xlim,max(pv_size)+xlim])
plt.ylim([y_lower_corr,y_upper_corr])

plt.savefig(figpath+'laplace_mcmc_corr'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')


### CORRELATION EPS = 2, INCLUDE RLR
print('\nPlotting performance comparison (Spearman\'s rank correlation, eps=2)...')
pre_npv = 'corr-e0-'
pre = 'corr-e2-'
meandata = []
stddata = []
for k in range(len(methods)):
	method = methods[k]
	if k < numbase:
		f = open(outpath+pre_npv+method+'.csv')
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][0]*np.ones_like(x))
		f.close()
		stddata.append(np.zeros_like(x))
	else:
		if k < numnpv or k == len(methods)-1:
			f = open(outpath+pre_npv+method+postmean)
		else:
			f = open(outpath+pre+method+postmean)
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()
		if k < numnpv or k == len(methods)-1:
			f = open(outpath+pre_npv+method+poststd)
		else:
			f = open(outpath+pre+method+poststd)
		reader = csv.reader(f,delimiter=',')
		stddata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()

fig = plt.figure(figsize=figsize,dpi=dpi)
ax = fig.add_subplot(111)
ax.set_xticks(pv_size)
ax.set_xticklabels(xlabels)

for k,c in zip(range(len(methods)),color):
	if k < numbase:
		ax.plot(x,meandata[k],'--',ms=ms,lw=lw1,label=methodsnames[k],color=c,dashes=dashesnpv)
	else:
		if k < numnpv:
			ax.plot(x,meandata[k],'s--',ms=ms,lw=lw1,label=methodsnames[k],color=c,mew=mew,mec=mec,dashes=dashesnpv)
		elif k < len(methods)-1:		
			ax.plot(x,meandata[k],'s-',ms=ms,lw=lw2,label=methodsnames[k],color=c,mew=mew,mec=mec)
		else:
			ax.plot(x,meandata[k],'s--',ms=ms,lw=lw2,label=methodsnames[k],color=c,mew=mew,mec=mec,dashes=dashesnpv)
		ax.errorbar(x,meandata[k],lw=lw,yerr=[stddata[k],stddata[k]],fmt='none',ecolor=c,capthick=lw,capsize=cs)
lgd = ax.legend(loc=lgdloc,numpoints=1,prop={'size':lgdfs},ncol=ncol,handlelength=hl,bbox_to_anchor=bbta)
lgd.get_frame().set_linewidth(lgdw)
lgd.get_frame().set_edgecolor("black")
plt.xlabel('Size of dataset (internal+external)')
plt.ylabel('Spearman\'s rank correlation coefficient')
if showtitle:
	plt.title('Comparison of different methods ('+r'$\epsilon$'+'=2.0)')
plt.xlim([min(pv_size)-xlim,max(pv_size)+xlim])
plt.ylim([y_lower_corr,y_upper_corr])

plt.savefig(figpath+'laplace_mcmc_corr_rlr'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')

### WPC-INDEX EPS = 2
print('\nPlotting performance comparison (wpc-index, eps=2)...')
pre_npv = 'wpc-e0-'
pre = 'wpc-e2-'
meandata = []
stddata = []
for k in range(len(methods)-1):
	method = methods[k]
	if k < numbase:
		f = open(outpath+pre_npv+method+'.csv')
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][0]*np.ones_like(x))
		f.close()
		stddata.append(np.zeros_like(x))
	else:
		if k < numnpv:
			f = open(outpath+pre_npv+method+postmean)
		else:
			f = open(outpath+pre+method+postmean)
		reader = csv.reader(f,delimiter=',')
		meandata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()
		if k < numnpv:
			f = open(outpath+pre_npv+method+poststd)
		else:
			f = open(outpath+pre+method+poststd)
		reader = csv.reader(f,delimiter=',')
		stddata.append(np.array(list(reader)).astype(float)[0][:len(x)])
		f.close()

fig = plt.figure(figsize=figsize,dpi=dpi)
ax = fig.add_subplot(111)
ax.set_xticks(pv_size)
ax.set_xticklabels(xlabels)

for k,c in zip(range(len(methods)-1),color):
	if k < numbase:
		ax.plot(x,meandata[k],'--',ms=ms,lw=lw1,label=methodsnames[k],color=c,dashes=dashesnpv)
	else:
		if k < numnpv:
			ax.plot(x,meandata[k],'s--',ms=ms,lw=lw1,label=methodsnames[k],color=c,mew=mew,mec=mec,dashes=dashesnpv)
		else:		
			ax.plot(x,meandata[k],'s-',ms=ms,lw=lw2,label=methodsnames[k],color=c,mew=mew,mec=mec)
		ax.errorbar(x,meandata[k],lw=lw,yerr=[stddata[k],stddata[k]],fmt='none',ecolor=c,capthick=lw,capsize=cs)
lgd = ax.legend(loc=lgdloc,numpoints=1,prop={'size':lgdfs},ncol=ncol,handlelength=hl,bbox_to_anchor=bbta)
lgd.get_frame().set_linewidth(lgdw)
lgd.get_frame().set_edgecolor("black")
plt.xlabel('Size of dataset (internal+external)')
plt.ylabel('Wpc-index')
if showtitle:
	plt.title('Comparison of different methods ('+r'$\epsilon$'+'=2.0)')
plt.xlim([min(pv_size)-xlim,max(pv_size)+xlim])
plt.ylim([y_lower_wpc,y_upper_wpc])

plt.savefig(figpath+'laplace_mcmc_wpc'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')

### OPTIMAL BUDGET SPLIT
print('\nPlotting budget split comparison...')
n = 500
d = 10
e = 2.0
f = open(outpath+'budgetsplit.csv')
reader = csv.reader(f,delimiter=',')
err = np.array(list(reader)).astype(float)
f.close()

p = np.arange(0.05,0.95,0.05)
p = [round(float(i),2) for i in p]
lenp = len(p)
bestind = np.unravel_index(err.argmax(),err.shape)
p1_best = p[bestind[0]]
p2_best = p[bestind[1]]
p3_best = round(1.0 - p1_best - p2_best,2)
s_best = err[bestind[0],bestind[1]]
print('Optimal privacy budget split (n='+str(n)+', d='+str(d)+', e='+str(e)+'):')
print('p1 =',p1_best,' for nxx')
print('p2 =',p2_best,' for nxy')
print('p3 =',p3_best,' for nyy')

table = []
for i in range(lenp):
	for j in range(lenp):
		p1 = p[i]
		p2 = p[j]
		p3 = round(1.0-p1-p2,2)
		if p3 <= 0.01:
			continue
		s = err[i,j]
		table.append((p1,p2,p3,s))
p3list = p
p2list = []
slist = []
for i in range(lenp):
	p3 = p3list[i]
	p2 = [t[1] for t in table if t[2] == p3]
	s = [t[3] for t in table if t[2] == p3]
	p2list.append(p2)
	slist.append(s)

fig = plt.figure(figsize=figsizebs,dpi=dpi)
ax = fig.add_subplot(111)
color = cm.jet(np.linspace(0,1,lenp))
for i,c in zip(range(lenp),color):
	if i%2 == 0:
		plt.plot(p2list[i],slist[i],'s-',lw=lw2,ms=ms,mew=mew,mec=mec,label=str(p3list[i]),color=c)
plt.plot(p2_best,s_best,'kx',ms=msbs,mew=mewbs,label='best')
plt.xlim([min(p)-0.05,max(p)+0.05])
plt.ylim([0.2,0.72])
lgd = ax.legend(loc='lower right',numpoints=1,prop={'size':lgdfs},title='privacy budget\n' +r'$\mathregular{share\; for\; n\overline{yy}}$',ncol=lgdncol)
plt.setp(lgd.get_title(), multialignment='center')
ax.get_legend().get_title().set_fontsize(lgdfs)
lgd.get_frame().set_linewidth(lgdw)
lgd.get_frame().set_edgecolor("black")
plt.xlabel(r'$\mathregular{privacy\; budget\; share\; for\; n\overline{xy}}$')
plt.ylabel('Spearman\'s rank correlation coefficient')
if showtitle:
	plt.title('Optimal privacy budget split')

plt.savefig(figpath+'laplace_mcmc_budgetsplit'+postfix+'.eps',format='eps',dpi=dpi,bbox_inches='tight')


