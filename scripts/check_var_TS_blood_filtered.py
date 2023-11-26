#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
check var of TS_Blood_filtered.h5ad
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
import pandas as pd
# # from scipy.sparse import csr_matrix
# import matplotlib.pyplot as plt
# import seaborn as sns

# set plot options
sc.set_figure_params(dpi=100, color_map = 'viridis_r')
sc.settings.verbosity = 1
sc.logging.print_header()
# plt.style.use('seaborn-white') # 'classic'

# input directories
i_dir = '/home/nikola/Project_Data/Python_data/terminal/Tabular_Sapiens/output_figshare'
v_dir = '/home/nikola/Project_Data/Python_data/Spyder/Tabular_Sapiens/output_figshare/meta'

# output directory
o_dir = os.path.join(wd, 'output_figshare')
# os.mkdir(o_dir)

# =============================================================================


# =============================================================================

fl = 'TS_Blood_filtered.h5ad'
fn = '/'.join((i_dir, fl))
adata = sc.read_h5ad(fn)
adata.shape

# var
f_v = 'var.csv'
fn_v = '/'.join((v_dir, f_v))
var_info = pd.read_csv(fn_v, index_col=0)
var_info.shape
# (49274, 47706)

# count how many of adata.var_names included in var_info.gene_symbol
lst = [x in var_info.gene_symbol for x in adata.var_names]
pd.Series(lst).sum()
# 47706

adata.var_names.unique().size

# save meta as csv files
adata.write_csvs(os.path.join(o_dir, 'blood_filtered'), skip_data=True)
