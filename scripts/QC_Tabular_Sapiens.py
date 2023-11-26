# -*- coding: utf-8 -*-
"""
quality control of the Tabular Sapiens dataset
$input all_cells.h5ad
# source https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
$output QC metrics & plots
# QC metrics
#  1. the number of genes expressed in the count matrix
#  2. the total counts per cell
#  3. the percentage of counts in mitochondrial genes
#  4. the percentage of counts in ribosomal genes
"""

# =============================================================================

# 1 set up Python session

# set working directory
import os
wd = '/home/nikola/Project_Data/Python_data/Spyder/Tabular_Sapiens'
os.chdir(wd)

# import packages
import scanpy as sc
# import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
# from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import seaborn as sns

# set plot options
sc.set_figure_params(dpi=100, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
plt.style.use('seaborn-white') # 'classic'

# input directories
i_dir = '/home/nikola/Documents/Rsch/resources/HCA/Tabula_Sapiens'
ds_dir = '/home/nikola/Project_Data/R_data/Tabula_Sapiens/output'

# output directory
o_dir = os.path.join(wd, 'output')
# os.mkdir(o_dir)

# =============================================================================


# =============================================================================

# 2 QC of the Tabular Sapiens dataset

# load data
fl = 'all_cells.h5ad'
fn = '/'.join((i_dir, fl))
adata = sc.read_h5ad(fn)

# 2.1 inspect the AnnData object
adata
"""
AnnData object with n_obs × n_vars = 483152 × 58559
    obs: 'tissue_in_publication', 'assay_ontology_term_id', 'donor_id', 
    'anatomical_information', 'n_counts_UMIs', 'n_genes', 'cell_ontology_class', 
    'free_annotation', 'manually_annotated', 'compartment', 'sex_ontology_term_id', 
    'is_primary_data', 'organism_ontology_term_id', 'disease_ontology_term_id', 
    'self_reported_ethnicity_ontology_term_id', 'development_stage_ontology_term_id', 
    'cell_type_ontology_term_id', 'tissue_ontology_term_id', 'suspension_type', 
    'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 
    'development_stage'
    var: 'feature_type', 'ensemblid', 'highly_variable', 'means', 'dispersions', 
    'dispersions_norm', 'mean', 'std', 'feature_is_filtered', 'feature_name', 
    'feature_reference', 'feature_biotype'
    uns: '_scvi', '_training_mode', 'compartment_colors', 'default_embedding', 
    'dendrogram_cell_type_tissue', 'dendrogram_computational_compartment_assignment', 
    'dendrogram_consensus_prediction', 'dendrogram_tissue_cell_type', 'donor_id_colors', 
    'donor_method_colors', 'hvg', 'method_colors', 'neighbors', 'schema_version', 
    'sex_colors', 'tissue_in_publication_colors', 'title', 'umap'
    obsm: 'X_pca', 'X_scvi', 'X_scvi_umap', 'X_umap'
    layers: 'decontXcounts'
    obsp: 'connectivities', 'distances'
"""
# n_vars = 58870 in dataset obtained from figshare
# raw_counts layer missing

# expression data
# The following numeric layers are available:
#     .layers["raw_counts"]: raw, not normalized counts.
#     .layers["decontXcounts"]: decontX corrected counts.
#     .raw.X: log1p normalized decontX corrected counts.
#     .X: log1p normalized and scaled decontX corrected counts.
# ref https://tabula-sapiens-portal.ds.czbiohub.org/organs
adata.layers['decontXcounts']
adata.X
adata.raw.X

# sample meta
adata.obs[:8] # pd DataFrame
adata.obs.columns # alternative adata.obs.keys()
# adata.obs.iloc[:6, :6]
adata.obs_names[:8]

# feature meta
adata.var[:8]
adata.var.columns
adata.var_names[:8]

# # highest expressing genes
# sc.pl.highest_expr_genes(adata, n_top=20, )
# # insufficient memory

# 2.2 quality control
# 2.2.1 make var names unique
adata.var_names_make_unique()
adata

# save meta as csv files
adata.write_csvs(os.path.join(o_dir, 'meta'), skip_data=True)

# save gene info as csv file
adata.var.columns
df_g = adata.var.loc[:,['ensemblid', 'feature_name']]
df_g.to_csv(os.path.join(o_dir, 'genes_TS_CELLxGENE.csv'), index=True)

# save data subset as csv file
# subset 1
adata.X[1000:1006,40000:40006]
df = pd.DataFrame(csr_matrix.todense(adata.X[1000:1006,40000:40006]))
df.index = adata.obs.index[1000:1006]
df.columns = adata.var.index[40000:40006]
df
df.to_csv(os.path.join(o_dir, 'X_c1k_g40k.csv'), index=True, header=True)
# save data subset & associated meta data
adata[1000:1006,40000:40006].write_csvs(os.path.join(o_dir, 'c1k_g40k'), skip_data=False)

# subset 2
# adata.X[8*10**4:(8*10**4+20),5*10**4:(5*10**4+20)]
df_2 = pd.DataFrame(csr_matrix.todense(adata.X[8*10**4:(8*10**4+20),5*10**4:(5*10**4+20)]))
df_2.index = adata.obs.index[8*10**4:(8*10**4+20)]
df_2.columns = adata.var.index[5*10**4:(5*10**4+20)]
df_2
df_2.to_csv(os.path.join(o_dir, 'X_c80k_g50k.csv'), index=True, header=True)
# save data subset & associated meta data
adata[8*10**4:(8*10**4+20),5*10**4:(5*10**4+20)].write_csvs(os.path.join(o_dir, 'c80k_g50k'), skip_data=False)

# # save gene expression data & meta as csv files
# adata.write_csvs(os.path.join(o_dir, 'all'), skip_data=False)
# # MemoryError: Unable to allocate 105. GiB for an array with shape (483152, 58559) and data type float32

# 2.2.2 QC metrics

# 2.2.2.1 inspect # expressed genes
# saved # expressed genes
ng = adata.obs.loc[:,'n_genes']
type(ng)
ng.shape
ng.index
ng
ng.median()
# 2158.0
ng.mean()
# 2526.0704333211906
ng.min()
# 200
ng.max()
# 34543
ng.std()
# 1414.8366970581403
# histogram
plt.hist(ng, bins=50)
plt.savefig(os.path.join(o_dir, 'histogram_n_genes_meta.pdf'))

# 2.2.2.2 inspect # gene counts
# saved # gene counts
nc = adata.obs.loc[:,'n_counts_UMIs']
nc.median()
# 7913.5
nc.mean()
# 52043.086
nc.min()
# 2500.0
nc.max()
# 22275940.0
nc.std()
# 304864.9

# histogram
plt.hist(nc, bins=50)
plt.savefig(os.path.join(o_dir, 'histogram_n_counts_meta.pdf'))
# density plot
sns.kdeplot(nc, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_n_counts_meta.pdf'))

# calculated # gene counts
# decontXcounts layer
nc_decontX = adata.layers["decontXcounts"].sum(axis=1) # row sums
type(nc_decontX)
# convert matrix to 1-dimensional ndarray
nc_decontX = np.ravel(nc_decontX)
nc_decontX = pd.Series(nc_decontX)
nc_decontX.shape
nc_decontX.index = ng.index # modify indices for index alignment
nc_decontX
# statistics
nc_decontX.median()
# 7911.0
nc_decontX.mean()
# 52033.168
nc_decontX.min()
# 2490.0
nc_decontX.max()
# 22275582.0
nc_decontX.std()
# 304822.56
# histogram
plt.hist(nc_decontX, bins=50)
plt.savefig(os.path.join(o_dir, 'histogram_n_counts_decontX.pdf'))
# density plot
sns.kdeplot(nc_decontX, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_n_counts_decontX.pdf'))
# calculate the difference
diff_nc_decontX = nc - nc_decontX
# statistics
diff_nc_decontX.median()
# 1.0
diff_nc_decontX.mean()
# 9.90175
diff_nc_decontX.min()
# 0
diff_nc_decontX.max()
# 10138.0
diff_nc_decontX.std()
# 94.92633
# histogram
plt.hist(diff_nc_decontX, bins=100, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_diff_n_counts_decontX.pdf'))
# density plot
sns.kdeplot(diff_nc_decontX, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_diff_n_counts_decontX.pdf'))

# raw.X
nc_raw_X = adata.raw.X.sum(axis=1) # row sums
type(nc_raw_X)
# convert matrix to 1-dimensional ndarray
nc_raw_X = np.ravel(nc_raw_X)
nc_raw_X = pd.Series(nc_raw_X)
nc_raw_X.shape
nc_raw_X.index = ng.index # modify indices for index alignment
nc_raw_X
# statistics
nc_raw_X.median()
# 8448.0
nc_raw_X.mean()
# 54954.05
nc_raw_X.min()
# 2511.0
nc_raw_X.max()
# 22305076.0
nc_raw_X.std()
# 313048.22
# histogram
plt.hist(nc_raw_X, bins=50)
plt.savefig(os.path.join(o_dir, 'histogram_n_counts_raw_X.pdf'))
# density plot
sns.kdeplot(nc_raw_X, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_n_counts_raw_X.pdf'))
# calculate the difference
diff_nc_raw_X = nc - nc_raw_X
# statistics
diff_nc_raw_X.median()
# -112.0
diff_nc_raw_X.mean()
# -2910.9788
diff_nc_raw_X.min()
# -6104775.0
diff_nc_raw_X.max()
# 10064.0
diff_nc_raw_X.std()
# 33662.793
# histogram
plt.hist(diff_nc_raw_X, bins=100, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_diff_n_counts_raw_X.pdf'))
# density plot
sns.kdeplot(diff_nc_raw_X, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_diff_n_counts_raw_X.pdf'))
# conclusion statistics identical to those of the Seurat object

# X
nc_X = adata.X.sum(axis=1) # row sums
type(nc_X)
# convert matrix to 1-dimensional ndarray
nc_X = np.ravel(nc_X)
nc_X = pd.Series(nc_X)
nc_X.shape
nc_X.index = ng.index # modify indices for index alignment
nc_X
# nc_X the same as total_counts column of meta after calling calculate_qc_metrics (see below)
nc_X.median()
# 6032.158
nc_X.mean()
# 6083.623
nc_X.min()
# 263.94144
nc_X.max()
# 135744.52
nc_X.std()
# 2047.9714
# histogram
plt.hist(nc_X, bins=50)
plt.savefig(os.path.join(o_dir, 'histogram_n_counts_X.pdf'))
# density plot
sns.kdeplot(nc_X, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_n_counts_X.pdf'))
# calculate the difference
diff_nc_X = nc - nc_X
# statistics
diff_nc_X.median()
# 1804.3965
diff_nc_X.mean()
# 45959.45
diff_nc_X.min()
# -54115.836
diff_nc_X.max()
# 22268050.0
diff_nc_X.std()
# 304895.47
# histogram
plt.hist(diff_nc_X, bins=100, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_diff_n_counts_X.pdf'))
# density plot
sns.kdeplot(diff_nc_X, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_diff_n_counts_X.pdf'))

# 2.2.2.3 plot QC metrics
# annotate the group of mitochondrial genes as 'mt'
adata.var.loc[:,'feature_name'].head()
# adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['mt'] = adata.var.feature_name.str.startswith('MT-')
# raw.X used for calculate_qc_metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, use_raw=True, log1p=False, inplace=True)
# if use_raw=False(default), adata.X, i.e. log1p normalized and scaled decontX corrected counts, 
# used for calculate_qc_metrics
adata.obs.pct_counts_mt
# pct_counts_mt calculated as total_counts_mt / total_counts * 100
adata.obs.pct_counts_mt.median()
# 7.9542966
adata.obs.pct_counts_mt.mean()
# 10.831316
adata.obs.pct_counts_mt.min()
# 0.0
adata.obs.pct_counts_mt.max()
# 10.849361
adata.obs.pct_counts_mt.std()
# 10.849361
# histogram
plt.hist(adata.obs.pct_counts_mt, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_pct_counts_mt.pdf'))
# density plot
sns.kdeplot(adata.obs.pct_counts_mt, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_pct_counts_mt.pdf'))

# annotate the group of ribosomal genes as 'rb'
adata.var['rb'] = adata.var.feature_name.str.match('RP[SL]') # alternative .contains('^RP[SL]')
sc.pp.calculate_qc_metrics(adata, qc_vars=['rb'], percent_top=None, use_raw=True, log1p=False, inplace=True)
adata.obs.pct_counts_rb
adata.obs.pct_counts_rb.median()
# 16.809708
adata.obs.pct_counts_rb.mean()
# 18.29946
adata.obs.pct_counts_rb.min()
# 0.0
adata.obs.pct_counts_rb.max()
# 78.21024
adata.obs.pct_counts_rb.std()
# 11.147203
# histogram
plt.hist(adata.obs.pct_counts_rb, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_pct_counts_rb.pdf'))
# density plot
sns.kdeplot(adata.obs.pct_counts_rb, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_pct_counts_rb.pdf'))

# violin plots of QC metrics
# n_genes
sc.pl.violin(adata, 
              ['n_genes', 'n_genes_by_counts'], 
              jitter=.4, multi_panel=True, show=False)
plt.savefig(os.path.join(o_dir, 'violin_plot_n_genes.png'))
# n_counts
sc.pl.violin(adata, 
              ['n_counts_UMIs', 'total_counts'], 
              jitter=.4, multi_panel=True, show=False)
plt.savefig(os.path.join(o_dir, 'violin_plot_n_counts.png'))
# pct_counts
sc.pl.violin(adata, 
              ['pct_counts_mt', 'pct_counts_rb'], 
              jitter=.4, multi_panel=True, show=False)
plt.savefig(os.path.join(o_dir, 'violin_plot_pct_counts.png'))

# scatter plots
# pct_counts_mt vs total_counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_mt_n_counts.png'))
# n_genes_by_counts vs total_counts
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_n_genes_n_counts.png'))
# pct_counts_rb vs total_counts
sc.pl.scatter(adata, x='total_counts', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_rb_n_counts.png'))
# pct_counts_rb vs pct_counts_mt
sc.pl.scatter(adata, x='pct_counts_mt', y='pct_counts_rb', show=False)
plt.savefig(os.path.join(o_dir, 'scatter_plot_pct_rb_pct_mt.png'))

# scatter plots with marginal histograms
# hexagonal bins
# n_genes_by_counts vs total_counts
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="hex"
)
plt.savefig(os.path.join(o_dir, 'scatter_plot_hist_hex_n_genes_n_counts.pdf'))
# add a linear regression fit and univariate KDE curves
# n_genes_by_counts vs total_counts
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="n_genes_by_counts",
    kind="reg"
)
plt.savefig(os.path.join(o_dir, 'scatter_plot_hist_reg_n_genes_n_counts.png'))
# pct_counts_mt vs total_counts
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="pct_counts_mt",
    kind="reg"
)
plt.savefig(os.path.join(o_dir, 'scatter_plot_hist_reg_pct_mt_n_counts.png'))
# pct_counts_rb vs total_counts
sns.jointplot(
    data=adata.obs,
    x="total_counts",
    y="pct_counts_rb",
    kind="reg"
)
plt.savefig(os.path.join(o_dir, 'scatter_plot_hist_reg_pct_rb_n_counts.png'))

# calculated # expressed genes
# raw.X
adata.obs.n_genes_by_counts.median()
# 2220.0
adata.obs.n_genes_by_counts.mean()
#  2620.3306785442264
adata.obs.n_genes_by_counts.min()
# 201
adata.obs.n_genes_by_counts.max()
# 34908
adata.obs.n_genes_by_counts.std()
# 1470.7687638441207
# histogram
plt.hist(adata.obs.n_genes_by_counts, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_n_genes_raw_X.pdf'))
# density plot
sns.kdeplot(adata.obs.n_genes_by_counts, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_diff_n_genes_raw_X.pdf'))
# calculate the difference
diff_ng_raw_X = ng - adata.obs.n_genes_by_counts
diff_ng_raw_X.median()
# -11.0
diff_ng_raw_X.mean()
# -94.2602452230354
diff_ng_raw_X.min()
# -9959
diff_ng_raw_X.max()
# 30
diff_ng_raw_X.std()
# 378.51121784656664
# histogram
plt.hist(diff_ng_raw_X, bins=50, density=False, alpha=.5, histtype='stepfilled', 
          color='steelblue', edgecolor='none')
plt.savefig(os.path.join(o_dir, 'histogram_diff_n_genes_raw_X.pdf'))
# density plot
sns.kdeplot(diff_ng_raw_X, fill=True, alpha=.5, linewidth=0)
plt.savefig(os.path.join(o_dir, 'density_plot_diff_n_genes_raw_X.pdf'))
# conclusion statistics slightly different from those of the those of the Seurat object

# 2.2.3 comparison with data of Seurat object
# load data subset of the Seurat object
X_s_sub = pd.read_csv(os.path.join(ds_dir, 'X_sub.csv'))
X_s_sub.shape
# calculate # gene counts
X_s_sub.sum(axis=1)
X_s_sub.sum(axis=1).min()
# 2511
X_s_sub.sum(axis=1).median()
# 25112.0
X_s_sub.sum(axis=1).max()
# 22305076
# compare # gene counts
X_s_sub.sum(axis=1)
"""
GGATCTACAAAGGCAC_TSP7_Blood_NA_10X_1_1                       12511
TAAGCCATCCGTCAAA_TSP7_Blood_NA_10X_2_1                       12511
ATGCATGAGATACATG_TSP7_SalivaryGland_Parotid_10X_1_2          25117
ATTTCACTCCGCGAGT_TSP7_Spleen_NA_10X_1_1                      25115
CTCAGGGGTGTCTTGA_TSP7_Tongue_Anterior_10X_1_2                25114
  
B107919_D14_S243.homo.gencode.v30.ERCC.chrM                1284484
TSP2_SI_distal_SS2_B114584_B133323_Epithelial_D16_S40       284480
CTGGACGCACTGTCCT_TSP2_Trachea_NA_10X_1_1                      8448
ACATCCCAGTACAACA_TSP2_Vasculature_Aorta_10X_1_2               8448
TSP10_Fat_MAT_SS2_B134181_B115064_Endothelial_G23_L003    22305076
Length: 118, dtype: int64
"""
select_obs = X_s_sub.index
adata.obs.total_counts.loc[select_obs]
"""
GGATCTACAAAGGCAC_TSP7_Blood_NA_10X_1_1                       12511.0
TAAGCCATCCGTCAAA_TSP7_Blood_NA_10X_2_1                       12511.0
ATGCATGAGATACATG_TSP7_SalivaryGland_Parotid_10X_1_2          25117.0
ATTTCACTCCGCGAGT_TSP7_Spleen_NA_10X_1_1                      25115.0
CTCAGGGGTGTCTTGA_TSP7_Tongue_Anterior_10X_1_2                25114.0
   
B107919_D14_S243.homo.gencode.v30.ERCC.chrM                1284484.0
TSP2_SI_distal_SS2_B114584_B133323_Epithelial_D16_S40       284480.0
CTGGACGCACTGTCCT_TSP2_Trachea_NA_10X_1_1                      8448.0
ACATCCCAGTACAACA_TSP2_Vasculature_Aorta_10X_1_2               8448.0
TSP10_Fat_MAT_SS2_B134181_B115064_Endothelial_G23_L003    22305076.0
Name: total_counts, Length: 118, dtype: float32
"""
# conclusion data from two souces are identical

# 2.2.4 preprocessing
# # set the .raw attribute
# adata.raw = adata

# total-count normalise the data to 10,000 reads per cell
sc.pp.normalize_total(adata, target_sum=1e4)
# layer : Optional[str] (default: None)
# Layer to normalize instead of X. If None, X is normalized.

# logarithmise the data
sc.pp.log1p(adata)
# adata.X used for log1p

# identify highly-variable genes
# layer : Optional[str] (default: None)
#     If provided, use adata.layers[layer] for expression values instead of adata.X.
# adata.X used for highly_variable_genes
sc.pp.highly_variable_genes(adata, min_mean=.0125, max_mean=3, min_disp=.5)
# insufficient memory if calling highly_variable_genes with adata.raw.to_adata()
# plot
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(os.path.join(o_dir, 'filter_genes_dispersion.png'))

# # Set the .raw attribute of the AnnData object to the normalized and logarithmized gene expression
# adata.raw = adata
# # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed
# sc.pp.regress_out(adata, ['n_counts_UMIs', 'pct_counts_mt'])
# # MemoryError: Unable to allocate 105. GiB for an array with shape (483152, 58559) and data type float32
# # Scale each gene to unit variance. Clip values exceeding standard deviation 10
# sc.pp.scale(adata, max_value=10)

# 2.2.5 PCA
sc.tl.pca(adata, svd_solver='arpack')

# scatter plot in PCA coordinates
colours = ['donor_id', 'assay', 'tissue_in_publication', 'sex', 'compartment']
for c in colours:
    print('colour:', c, sep=' ')
    sc.pl.pca(adata, color=c, show=False)
    f = f"scatter_plot_PCA_{c}.pdf"
    plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
    print('completed.')
    
# PCA variance ratio
sc.pl.pca_variance_ratio(adata, log=True, show=False)
plt.savefig(os.path.join(o_dir, 'PCA_variance_ratio_test.pdf'))

# 2.2.6 clustering
# computing the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# embedding the neighborhood graph using UMAP
sc.tl.umap(adata)
for c in colours:
    print('colour:', c, sep=' ')
    sc.pl.umap(adata, color=c, show=False)
    f = f"UMAP_{c}.pdf"
    plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
    print('completed.')

# clustering the neighborhood graph
sc.tl.leiden(adata)
colours_2 = ['leiden', 'donor_id', 'assay', 'tissue_in_publication', 'sex', 'compartment']
for c in colours_2:
    print('colour:', c, sep=' ')
    sc.pl.umap(adata, color=c, show=False)
    f = f"UMAP_Leiden_clustering_{c}.pdf"
    plt.savefig(os.path.join(o_dir, f), bbox_inches='tight')
    print('completed.')

# # In some occasions, you might still observe disconnected clusters and similar connectivity violations
# # They can usually be remedied by running
# sc.tl.paga(adata)
# sc.pl.paga(adata, plot=False) # remove `plot=False` if you want to see the coarse-grained graph
# sc.tl.umap(adata, init_pos='paga')

# # basic filtering
# adata = adata
# sc.pp.filter_cells(adata, min_genes=200)
# sc.pp.filter_genes(adata, min_cells=3)
# # insufficient memory

# end of the session

# =============================================================================
