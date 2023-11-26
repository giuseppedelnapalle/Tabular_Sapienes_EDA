#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QC of the subset of Tabular Sapiens dataset & export filtered data as csv files
$input TS_{tissue}.h5ad
# source https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219
$output QC plots, filtered subset of specified tissue as csv files
"""


# =============================================================================

# 1 set up Python session

import os
import scanpy as sc
import argparse
import time
# import numpy as np
# import pandas as pd
# from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.use('TkAgg')
# supported values are ['GTK3Agg', 'GTK3Cairo', 'GTK4Agg', 'GTK4Cairo', 'MacOSX', 'nbAgg', 'QtAgg', 
# 'QtCairo', 'Qt5Agg', 'Qt5Cairo', 'TkAgg', 'TkCairo', 'WebAgg', 'WX', 'WXAgg', 'WXCairo', 'agg', 
# 'cairo', 'pdf', 'pgf', 'ps', 'svg', 'template']

# working directory
wd = '/home/nikola/Project_Data/Python_data/terminal/Tabular_Sapiens'
# output directory
o_dir = os.path.join(wd, 'output_figshare')

# argument parser
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file', metavar='h5ad file', type=str, required=False, help="path to h5ad file as input")
parser.add_argument('-t', '--tissue', metavar='tissue or organ', type=str, required=True, 
                    help="tissue to select obs from the full dataset. choose from Lymph_Node, Blood, \
                    Lung, Spleen, Thymus, Muscle, Salivary_Gland, Bladder, Fat, Prostate, Vasculature, \
                        Tongue, Large_Intestine, Pancreas, Small_Intestine, Bone_Marrow, Heart, Mammary. \
                            Eye, Kidney, Trachea, Skin, Uterus, Liver")
parser.add_argument('-o', '--output', metavar='output path',  nargs='?', const=o_dir, default=o_dir, type=str, required=False,help="output path")

# parse command-line arguments
args = parser.parse_args()

# input file
fl = f"TS_{args.tissue}.h5ad"
fn = '/'.join((args.output, fl))

args.file = args.file if args.file else fn

# create output directory if it does not exist
if not os.path.exists(args.output):
    os.mkdir(args.output)

os.chdir(args.output)

# plot directory
p_dir = os.path.join(args.output, f"plots_{args.tissue}")
if not os.path.exists(p_dir):
    os.mkdir(p_dir)
    
# display arguments
# print("args=%s" % args)
print("args.file=%s" % args.file)
print("args.tissue=%s" % args.tissue)
print("args.output=%s" % args.output)


# =============================================================================


# =============================================================================

# 2 QC the Tabular Sapiens germ line dataset

# 2.1 load data
adata = sc.read_h5ad(args.file)

# 2.2 preprocessing
print(adata)

# expression data
# The following numeric layers are available:
#     .layers["raw_counts"]: raw, not normalized counts.
#     .layers["decontXcounts"]: decontX corrected counts.
#     .raw.X: log1p normalized decontX corrected counts.
#     .X: log1p normalized and scaled decontX corrected counts.
# ref https://tabula-sapiens-portal.ds.czbiohub.org/organs

# make var names unique
print('make var names unique.')
adata.var_names_make_unique()
print(adata)

# reset gene expression data
adata_new = adata.copy()
adata_new.X = adata_new.layers['raw_counts']

# pd.Series(np.ravel(adata_new.X.sum(axis=1))).median() # row sums
# pd.Series(np.ravel(adata_new.layers['raw_counts'].sum(axis=1))).median()

# pd.Series(np.ravel(adata.X.sum(axis=1))).median()
# pd.Series(np.ravel(adata.layers['raw_counts'].sum(axis=1))).median()

# Show those genes that yield the highest fraction of counts in each single cell, across all cells
sc.pl.highest_expr_genes(adata_new, n_top=20, show=False)
plt.savefig(os.path.join(p_dir, 'highest_expr_genes_raw_counts.pdf'), bbox_inches='tight')
plt.close()

# basic filtering
sc.pp.filter_cells(adata_new, min_genes=200)
sc.pp.filter_genes(adata_new, min_cells=3)

# percentage of mitochondrial genes
adata_new.var['mt'] = adata_new.var_names.str.startswith('MT-')
# .X used for calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata_new, qc_vars=['mt'], percent_top=None, use_raw=False, log1p=True, inplace=True)
# adata_new.obs.pct_counts_mt.median()
# histogram
plt.hist(adata_new.obs.pct_counts_mt, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(p_dir, 'histogram_pct_counts_mt_raw_counts.pdf'))
plt.close()
# density plot
sns.kdeplot(adata_new.obs.pct_counts_mt, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(p_dir, 'density_plot_pct_counts_mt_raw_counts.pdf'))
plt.close()

# percentage of ribosomal genes
adata_new.var['rb'] = adata_new.var_names.str.match('RP[SL]')
# .X used for calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata_new, qc_vars=['rb'], percent_top=None, use_raw=False, log1p=True, inplace=True)
# adata_new.obs.pct_counts_rb.median()
# histogram
plt.hist(adata_new.obs.pct_counts_rb, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(p_dir, 'histogram_pct_counts_rb_raw_counts.pdf'))
plt.close()
# density plot
sns.kdeplot(adata_new.obs.pct_counts_rb, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(p_dir, 'density_plot_pct_counts_rb_raw_counts.pdf'))
plt.close()

# violin plot
sc.pl.violin(adata_new, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb'],
             jitter=0.4, multi_panel=True, show=False)
plt.savefig(os.path.join(p_dir, 'violin_plot_qc_metrics_raw_counts.pdf'))
plt.close()

# scatter plots
# pct_counts_mt vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='pct_counts_mt', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_pct_mt_n_counts_raw_counts.pdf'))
plt.close()
# n_genes_by_counts vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='n_genes_by_counts', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_n_genes_n_counts_raw_counts.pdf'))
plt.close()
# pct_counts_rb vs total_counts
sc.pl.scatter(adata_new, x='total_counts', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_pct_rb_n_counts_raw_counts.pdf'))
plt.close()
# pct_counts_rb vs pct_counts_mt
sc.pl.scatter(adata_new, x='pct_counts_mt', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_pct_rb_pct_mt_raw_counts.pdf'))
plt.close()
# log1p
# pct_counts_mt vs log1p_total_counts
sc.pl.scatter(adata_new, x='log1p_total_counts', y='pct_counts_mt', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_pct_mt_log1p_n_counts_raw_counts.pdf'))
plt.close()
# log1p_n_genes_by_counts vs log1p_total_counts
sc.pl.scatter(adata_new, x='log1p_total_counts', y='log1p_n_genes_by_counts', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_log1p_n_genes_log1p_n_counts_raw_counts.pdf'))
plt.close()
# pct_counts_rb vs total_counts
sc.pl.scatter(adata_new, x='log1p_total_counts', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(p_dir, 'scatter_plot_pct_rb_log1p_n_counts_raw_counts.pdf'))
plt.close()

# filtering
# total_counts
# adata_new[adata_new.obs.total_counts > 2500, :]
adata_new = adata_new[adata_new.obs.total_counts > 2500, :]
# pct_counts_mt
mad_mt = (adata_new.obs.pct_counts_mt - adata_new.obs.pct_counts_mt.mean()).abs().mean()
thres_mt = adata_new.obs.pct_counts_mt.mean() + mad_mt*3 # 3 median absolute deviations from the mean
# adata_new[adata_new.obs.pct_counts_mt < thres_mt, :]
adata_new = adata_new[adata_new.obs.pct_counts_mt < thres_mt, :]

# adata_new.raw = adata_new

# Total-count normalize (library-size correct) the data matrix X
sc.pp.normalize_total(adata_new, target_sum=1e4)
# adata_new.layers['lib_size_norm'] = adata_new.X

# logarithmize the data
sc.pp.log1p(adata_new)
# adata_new.layers['log1p_norm'] = adata_new.X

# Identify highly-variable genes
sc.pp.highly_variable_genes(adata_new, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata_new, show=False)
plt.savefig(os.path.join(p_dir, 'highly_variable_genes_log1p_norm.pdf'))
plt.close()

# 2.3 PCA
sc.tl.pca(adata_new, svd_solver='arpack')

# scatter plot in PCA coordinates
colours = ['donor', 'method', 'organ_tissue', 'gender', 'compartment']
for c in colours:
    print('colour:', c, sep=' ')
    sc.pl.pca(adata_new, color=c, show=False)
    f = f"scatter_plot_PCA_{c}_log1p_norm.pdf"
    plt.savefig(os.path.join(p_dir, f), bbox_inches='tight')
    plt.close()
    print('completed.')

# PCA variance ratio
sc.pl.pca_variance_ratio(adata_new, log=True, show=False)
plt.savefig(os.path.join(p_dir, 'PCA_variance_ratio_test_log1p_norm.pdf'))
plt.close()

# 2.4 clustering
# computing the neighborhood graph
try:
    sc.pp.neighbors(adata_new, n_neighbors=10, n_pcs=40)
except ValueError:
    print('skip clustering.')
else:
    sc.tl.umap(adata_new)
    for c in colours:
        print('colour:', c, sep=' ')
        sc.pl.umap(adata, color=c, show=False)
        f = f"UMAP_{c}_log1p_norm.pdf"
        plt.savefig(os.path.join(p_dir, f), bbox_inches='tight')
        plt.close()
        print('completed.')

    # clustering the neighborhood graph
    sc.tl.leiden(adata_new)
    for c in colours:
        print('colour:', c, sep=' ')
        sc.pl.umap(adata, color=c, show=False)
        f = f"UMAP_Leiden_clustering_{c}_log1p_norm.pdf"
        plt.savefig(os.path.join(p_dir, f), bbox_inches='tight')
        plt.close()
        print('completed.')


# 2.5 filter out lowly expressed genes
# adata_s = adata[adata_new.obs_names,adata_new.var_names] # view of adata
# filter genes based on dispersions_norm
# print('filter genes based on normalised dispersions.')
# thres_g = adata_s.var.dispersions_norm.median()
# adata_s = adata_s[:,adata_s.var.dispersions_norm >= thres_g]
# print(adata_s)

# 2.6 write h5ad file
f = f"TS_{args.tissue}_filtered.h5ad"
print('tissue:', args.tissue, sep=' ')
adata_s = adata[adata_new.obs_names,adata_new.var_names] # view of adata
print(adata_s)
print('writing filtered data as h5ad file.')
begin = time.time()
adata_s.write(os.path.join(args.output, f), compression='gzip')
end = time.time()
print(f"Total runtime of the program is {end - begin}.")


# =============================================================================
