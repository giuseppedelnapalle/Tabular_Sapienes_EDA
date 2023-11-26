#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
partitioning of the Tabular Sapiens dataset
$input TabulaSapiens.h5ad
# source https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219
$output gene expression data subset of cells from distinct organs
# format csv
"""

# =============================================================================

# 1 set up Python session

# set working directory
import os
wd = '/home/nikola/Project_Data/Python_data/Spyder/Tabular_Sapiens'
os.chdir(wd)

# import packages
import scanpy as sc
import time
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

# 2 partition the Tabular Sapiens dataset

# load data
fl = 'TabulaSapiens.h5ad'
fn = '/'.join((i_dir, fl))
# adata = sc.read_h5ad(fn)
begin = time.time()
adata = sc.read_h5ad(fn, backed='r')
end = time.time()
print(f"Total runtime of the program is {end - begin}")
# Total runtime of the program is 273.7775900363922
adata.isbacked
adata.filename

# inspect annotation data
adata.obs.columns
# gender
adata.obs.gender.value_counts()
"""
male      245068
female    238084
Name: gender, dtype: int64
"""
# donor
adata.obs.donor.value_counts()
"""
TSP14    171719
TSP2     115557
TSP7      57127
TSP1      39431
TSP10     23463
TSP4      20930
TSP8      15067
TSP12     11505
TSP6       7978
TSP9       6699
TSP15      5611
TSP5       2592
TSP11      2488
TSP3       2447
TSP13       538
Name: donor, dtype: int64
"""
# method
adata.obs.method.value_counts()
"""
10X          456101
smartseq2     27051
Name: method, dtype: int64
"""
# compartment
adata.obs.compartment.value_counts()
"""
immune         264824
epithelial     104148
stromal         82478
endothelial     31691
germ line          11
Name: compartment, dtype: int6
"""
# organ_tissue
adata.obs.organ_tissue.value_counts()
"""
Lymph_Node         53275
Blood              50115
Lung               35682
Spleen             34004
Thymus             33664
Muscle             30746
Salivary_Gland     27199
Bladder            24583
Fat                20263
Prostate           16375
Vasculature        16037
Tongue             15020
Large_Intestine    13680
Pancreas           13497
Small_Intestine    12467
Bone_Marrow        12297
Heart              11505
Mammary            11375
Eye                10650
Kidney              9641
Trachea             9522
Skin                9424
Uterus              7124
Liver               5007
Name: organ_tissue, dtype: int64
"""
adata.obs.organ_tissue.value_counts().size
# 24
# cell_ontology_class
adata.obs.cell_ontology_class.value_counts()
"""
macrophage                               35204
fibroblast                               31999
b cell                                   19807
memory b cell                            15676
mesenchymal stem cell                    15459
 
respiratory mucous cell                      8
retina horizontal cell                       6
double-positive, alpha-beta thymocyte        6
retinal ganglion cell                        5
myeloid dendritic cell                       3
Name: cell_ontology_class, Length: 177, dtype: int64
"""
adata.obs.cell_ontology_class.value_counts().size
# 177

try:
    # export data subsets as csv files
    
    # drop = ['Salivary_Gland', 'Vasculature', 'Tongue', 'Heart']
    # select_tissues = adata.obs.organ_tissue.value_counts().drop(labels=drop).index
    # for t in select_tissues:
    #     print('tissue:', t, sep=' ')
    #     adata_s = adata[adata.obs.organ_tissue == t]
    #     begin = time.time()
    #     adata_s.write_csvs(os.path.join(o_dir, t), skip_data=False)
    #     end = time.time()
    #     print(f"Total runtime of the program is {end - begin}")
        
    # X
    # t = 'Blood'
    # adata_s = adata[adata.obs.organ_tissue == t]
    # begin = time.time()
    # adata_s.write_csvs(os.path.join(o_dir, t), skip_data=False)
    # end = time.time()
    # print(f"Total runtime of the program is {end - begin}")
    
    # h5ad
    t = 'Blood'
    f = f"TS_{t}.h5ad"
    adata_s = adata[adata.obs.organ_tissue == t]
    begin = time.time()
    adata_s.write(os.path.join(o_dir, f), compression='gzip')
    end = time.time()
    print(f"Total runtime of the program is {end - begin}")
    
finally:
    adata.file.close()
# Total runtime of the program is 128.21598839759827
# end of the session

# =============================================================================
