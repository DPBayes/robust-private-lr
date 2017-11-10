#!/bin/env python3
# Differentially private Bayesian linear regression 
# Arttu Nieminen 2016-2017
# University of Helsinki Department of Computer Science
# Helsinki Institute of Information Technology HIIT

# GDSC/drug sensitivity data
# Selects the most relevant genes from gene expression based on expert knowledge.
# Saves the reduced gene expression and drug response data into csv files.

import scipy.io as sio
import numpy as np
from operator import itemgetter
import pandas

datapath = 'drugsens_data/'
SelGenes70 = sio.loadmat(datapath+'SelGenes70.mat')['SelGenes70']
GenesMutations = sio.loadmat(datapath+'GenesMutations.mat')['GenesMutations']
MutationCnts = sio.loadmat(datapath+'MutationCnts.mat')['MutationCnts']
#GeneNames = sio.loadmat(datapath+'GeneNames.mat')['GeneNames']
geneexpr = pandas.read_hdf("data/GDSC_geneexpr.h5",
                           'rma_gene_expressions')


# Rank genes based on mutation counts from highest to lowest
#GeneNamesList = [n[0][0] for n in GeneNames]
GenesMutationsList = [n[0][0] for n in GenesMutations]
GeneIndex = []
for row in SelGenes70:
	genes = row[0][0].split('_')
	for gene in genes:
		#if gene in GeneNamesList:
		if gene in geneexpr.columns:
			#i = GeneNamesList.index(gene)
			i = geneexpr.columns.get_loc(gene)
			if gene in GenesMutationsList:
				j = GenesMutationsList.index(gene)
				GeneIndex.append((i,MutationCnts[j][0]))
			else:
				GeneIndex.append((i,0))
GeneIndex_sorted = sorted(GeneIndex,key=itemgetter(1),reverse=True)
RankedGenesInd = [t[0] for t in GeneIndex_sorted]
print('Found',len(RankedGenesInd),'/ 70 preselected genes in gene expression.')
print('Out of them,',np.count_nonzero([t[1] for t in GeneIndex]),'had mutations.')

# Pick out the most important genes from gene expression, in order of importance from high to low
x = geneexpr.as_matrix()[:,RankedGenesInd]

# Save into csv
np.savetxt('GeneExpressionReducted.csv',x,delimiter=',')

# Drugsens data
drugres = pandas.read_hdf("data/GDSC_drugres.h5", 'drug_responses')
y = drugres.as_matrix()

# Save into csv
np.savetxt('DrugResponse.csv',y,delimiter=',')
