# Script for GDSC data conversion
# Teppo Niinimäki 2017

import pandas
import numpy as np

geneexpr_file = "data/GDSC/sanger1018_brainarray_ensemblgene_rma.txt"
ensembl_name_file = "data/GDSC/ensembl_gene_ids.csv"
drugres_file = "data/GDSC/v17_fitted_dose_response.csv"
data_output_dir = "data"

# load gene expression data
# note: Some (4) samples appear twice, for instance 1503362.
#       The second appearance is renamed to 1503362.1, which is then later ignored
#       when finding common samples with the gene expression data.
geneexpr = pandas.read_csv(geneexpr_file, sep='\t', header=0, index_col=0)

# load ensembl gene ids to gene names conversion table
gene_names = pandas.read_csv(ensembl_name_file, header=0, index_col=0)

# find index genes by their names instead of ensembl ids,
# drop those that do not have a name
common_gene_ids = geneexpr.index.intersection(gene_names.index)
geneexpr = geneexpr.loc[common_gene_ids]
geneexpr.index = pandas.Index(gene_names.loc[common_gene_ids, "Gene name"].values, name='gene_name')

# transpose so that columns=genes and rows=samples
geneexpr = geneexpr.transpose()
geneexpr.index.name = 'sample_id'

# load drug response data
drugres_data = pandas.read_csv(drugres_file)
drugres = drugres_data.pivot(index='COSMIC_ID', columns='DRUG_ID',
                             values='LN_IC50')
drugres.index = pandas.Index(drugres.index.map(str), name='sample_id')
drugres.columns.name = 'drug_id'

# find and select commmon sample ids
common_samples = drugres.index.intersection(geneexpr.index)
drugres = drugres.loc[common_samples]
geneexpr = geneexpr.loc[common_samples]


#geneexpr = geneexpr_data.as_matrix()
#genes = geneexpr_data.columns

# export to hdf5 format
geneexpr.to_hdf(data_output_dir + "/GDSC_geneexpr.h5",
                'rma_gene_expressions', mode='w')
drugres.to_hdf(data_output_dir + "/GDSC_drugres.h5",
               'drug_responses', mode='w')


#drugres_data['DRUG_ID'].unique()
#drugres_data['COSMIC_ID'].unique()
#geneexpr_data.index.values

#geneexpr = pandas.read_hdf("GDSC_geneexpr.h5", 'rma_gene_expressions')
