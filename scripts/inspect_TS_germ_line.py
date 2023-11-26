#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inspect TS_germ line.h5ad downloaded from figshare
https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219
"""


# =============================================================================

# 1 set up Python session

# set working directory
import os
wd = '/home/nikola/Project_Data/Python_data/Spyder/Tabular_Sapiens'
os.chdir(wd)

# import packages
import scanpy as sc
# import time
# import anndata
# import numpy as np
# import pandas as pd
# from scipy.sparse import csr_matrix

# input directories
i_dir = '/home/nikola/Documents/Rsch/resources/HCA/Tabula_Sapiens'

# output directory
o_dir = os.path.join(wd, 'output_figshare')
# os.mkdir(o_dir)

# =============================================================================


# =============================================================================

# 2 inspect the Tabular Sapiens germ line dataset

# 2.1 load data
fl = 'TS_germ line.h5ad'
fn = '/'.join((i_dir, fl))
# adata = sc.read_h5ad(fn)

# # test
# try:
#     adata = sc.read_h5ad(fn, backed='r')
#     adata.isbacked
#     adata.filename
#     # 2.1 inspect the AnnData object
#     adata
# finally:
#     adata.file.close()

# # test 2
# with sc.read_h5ad(fn, backed='r') as adata:
#     adata.isbacked
#     adata.filename
#     adata

# # AttributeError: __enter__

# test 3
adata = sc.read_h5ad(fn, backed='r')
adata.isbacked
adata.filename

# 2.2 inspect the AnnData object
adata
"""
AnnData object with n_obs × n_vars = 11 × 58870
    obs: 'organ_tissue', 'method', 'donor', 'anatomical_information', 'n_counts_UMIs', 
    'n_genes', 'cell_ontology_class', 'free_annotation', 'manually_annotated', 'compartment', 
    'gender'
    var: 'gene_symbol', 'feature_type', 'ensemblid', 'highly_variable', 'means', 
    'dispersions', 'dispersions_norm', 'mean', 'std'
    uns: '_scvi', '_training_mode', 'dendrogram_cell_type_tissue', 
    'dendrogram_computational_compartment_assignment', 'dendrogram_consensus_prediction', 
    'dendrogram_tissue_cell_type', 'donor_colors', 'donor_method_colors', 'hvg', 
    'method_colors', 'neighbors', 'organ_tissue_colors', 'sex_colors', 'tissue_colors', 'umap'
    obsm: 'X_pca', 'X_scvi', 'X_scvi_umap', 'X_umap'
    layers: 'decontXcounts', 'raw_counts'
    obsp: 'connectivities', 'distances'
"""

# cell meta
adata.obs.columns
# gene meta
adata.var.columns

# export gene info
df_g = adata.var.loc[:,['gene_symbol', 'ensemblid']]
df_g.to_csv(os.path.join(o_dir, 'genes_TS_figshare.csv'), index=True)



adata.file.close()

# =============================================================================
