# Convert data from hdf5 files to csv files
# GeneNames, GeneExpression, DrugResponse
# Afterwards run matlab csvtomat.m to convert csv data into mat

import pandas
import numpy as np
import csv

datapath = '../Data/ver3/'

# Import data
g = pandas.read_hdf(datapath+"GDSC_geneexpr.h5",'rma_gene_expressions')
d = pandas.read_hdf(datapath+"GDSC_drugres.h5", 'drug_responses')

# Gene names as a single column list of length 17490
l = g.columns.tolist()
with open(datapath+'GeneNames.csv','w') as f:
	writer = csv.writer(f)
	for gene in l:
		writer.writerow([gene])

# Gene expression as a 985 x 17490 matrix
x = g.as_matrix()
np.savetxt(datapath+'GeneExpression.csv',x,delimiter=',')

# Drug response
y = d.as_matrix()
np.savetxt(datapath+'DrugResponse.csv',y,delimiter=',')
