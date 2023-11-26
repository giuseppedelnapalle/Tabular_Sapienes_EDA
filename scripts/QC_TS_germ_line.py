#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QC of the Tabular Sapiens germ line dataset
$input TS_germ line.h5ad
# source https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219
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
import numpy as np
import pandas as pd
# from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns

# input directories
i_dir = '/home/nikola/Documents/Rsch/resources/HCA/Tabula_Sapiens'

# output directory
o_dir = os.path.join(wd, 'output_figshare')
# os.mkdir(o_dir)

# =============================================================================


# =============================================================================

# 2 QC the Tabular Sapiens germ line dataset

# 2.1 load data
fl = 'TS_germ line.h5ad'
fn = '/'.join((i_dir, fl))
adata = sc.read_h5ad(fn)

# # test
# adata = sc.read_h5ad(fn, backed='r')
# adata.isbacked
# adata.filename
# adata.file.close()

# 2.2 preprocessing
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

# expression data
# The following numeric layers are available:
#     .layers["raw_counts"]: raw, not normalized counts.
#     .layers["decontXcounts"]: decontX corrected counts.
#     .raw.X: log1p normalized decontX corrected counts.
#     .X: log1p normalized and scaled decontX corrected counts.
# ref https://tabula-sapiens-portal.ds.czbiohub.org/organs

# make var names unique
adata.var_names_make_unique()
adata

# reset gene expression data
adata_new = adata.copy()
adata_new.X = adata_new.layers['raw_counts']

pd.Series(np.ravel(adata_new.X.sum(axis=1))).median() # row sums
pd.Series(np.ravel(adata_new.layers['raw_counts'].sum(axis=1))).median()

pd.Series(np.ravel(adata.X.sum(axis=1))).median()
pd.Series(np.ravel(adata.layers['raw_counts'].sum(axis=1))).median()

# Show those genes that yield the highest fraction of counts in each single cell, across all cells
sc.pl.highest_expr_genes(adata_new, n_top=20, show=False)
plt.savefig(os.path.join(o_dir, 'highest_expr_genes_raw_counts.pdf'))

# basic filtering
sc.pp.filter_cells(adata_new, min_genes=200)
sc.pp.filter_genes(adata_new, min_cells=1)


# percentage of mitochondrial genes
adata_new.var['mt'] = adata_new.var_names.str.startswith('MT-')
# .X used for calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata_new, qc_vars=['mt'], percent_top=None, use_raw=False, log1p=True, inplace=True)
adata_new.obs.pct_counts_mt.median()
# 4.048381
# histogram
plt.hist(adata_new.obs.pct_counts_mt, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_pct_counts_mt_raw_counts.pdf'))
plt.close()
# density plot
sns.kdeplot(adata_new.obs.pct_counts_mt, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_pct_counts_mt_raw_counts.pdf'))
plt.close()

# percentage of ribosomal genes
adata_new.var['rb'] = adata_new.var_names.str.match('RP[SL]')
# .X used for calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata_new, qc_vars=['rb'], percent_top=None, use_raw=False, log1p=True, inplace=True)
adata_new.obs.pct_counts_rb.median()
# 29.785204
# histogram
plt.hist(adata_new.obs.pct_counts_rb, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_pct_counts_rb_raw_counts.pdf'))
plt.close()
# density plot
sns.kdeplot(adata_new.obs.pct_counts_rb, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_pct_counts_rb_raw_counts.pdf'))
plt.close()

# violin plot
sc.pl.violin(adata_new, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb'],
             jitter=0.4, multi_panel=True, show=False)
plt.savefig(os.path.join(o_dir, 'violin_plot_qc_metrics_raw_counts.pdf'))

# scatter plots
# pct_counts_mt vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='pct_counts_mt', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_mt_n_counts_raw_counts.pdf'))
# n_genes_by_counts vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='n_genes_by_counts', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_n_genes_n_counts_raw_counts.pdf'))
# pct_counts_rb vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_rb_n_counts_raw_counts.pdf'))
# pct_counts_rb vs pct_counts_mt
sc.pl.scatter(adata_new, x='pct_counts_mt', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_rb_pct_mt_raw_counts.pdf'))

# filtering
# adata_new[adata_new.obs.n_genes_by_counts < 4500, :]
adata_new = adata_new[adata_new.obs.n_genes_by_counts < 4500, :]
# adata_new[adata_new.obs.pct_counts_mt < 10, :]
adata_new = adata_new[adata_new.obs.pct_counts_mt < 10, :]

adata_new.raw = adata_new

# Total-count normalize (library-size correct) the data matrix X
sc.pp.normalize_total(adata_new, target_sum=1e4)
adata_new.layers['lib_size_norm'] = adata_new.X

# logarithmize the data
sc.pp.log1p(adata_new)
adata_new.layers['log1p_norm'] = adata_new.X

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata_new, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_new, show=False)
plt.savefig(os.path.join(o_dir, 'highly_variable_genes_log1p_norm.pdf'))

# 2.3 PCA
sc.tl.pca(adata_new, svd_solver='arpack')

# scatter plot in PCA coordinates
colours = ['donor', 'method', 'organ_tissue', 'gender', 'compartment']
for c in colours:
    print('colour:', c, sep=' ')
    sc.pl.pca(adata_new, color=c, show=False)
    f = f"scatter_plot_PCA_{c}_log1p_norm.pdf"
    plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
    print('completed.')

# PCA variance ratio
sc.pl.pca_variance_ratio(adata_new, log=True, show=False)
plt.savefig(os.path.join(o_dir, 'PCA_variance_ratio_test_log1p_norm.pdf'))

# # 2.4 clustering
# # computing the neighborhood graph
# sc.pp.neighbors(adata_new, n_neighbors=10, n_pcs=40)
# # ValueError: `X_pca` does not have enough PCs. Rerun `sc.pp.pca` with adjusted `n_comps`.

# # embedding the neighborhood graph using UMAP
# sc.tl.umap(adata_new)
# for c in colours:
#     print('colour:', c, sep=' ')
#     sc.pl.umap(adata, color=c, show=False)
#     f = f"UMAP_{c}_log1p_norm.pdf"
#     plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
#     print('completed.')

# # clustering the neighborhood graph
# sc.tl.leiden(adata_new)
# colours_2 = ['leiden', 'donor', 'method', 'organ_tissue', 'gender', 'compartment']
# for c in colours_2:
#     print('colour:', c, sep=' ')
#     sc.pl.umap(adata, color=c, show=False)
#     f = f"UMAP_Leiden_clustering_{c}_log1p_norm.pdf"
#     plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
#     print('completed.')

# =============================================================================
