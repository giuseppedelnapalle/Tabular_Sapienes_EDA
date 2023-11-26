#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exporting subset of the Tabular Sapiens dataset
$input TabulaSapiens.h5ad from figshare
# source https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219
$output anndata object of specified subset
# format h5ad
"""

import os
import scanpy as sc
import argparse
import time

# input directory
i_dir = '/home/nikola/Documents/Rsch/resources/HCA/Tabula_Sapiens'

# working directory
wd = '/home/nikola/Project_Data/Python_data/terminal/Tabular_Sapiens'
# output directory
o_dir = os.path.join(wd, 'output_figshare')

# input file
fl = 'TabulaSapiens.h5ad'
fn = '/'.join((i_dir, fl))

# argument parser
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file', metavar='h5ad file', nargs='?', const=fn, default=fn, type=str, required=False, help="path to h5ad file as input")
parser.add_argument('-t', '--tissue', metavar='tissue or organ', type=str, required=True, 
                    help="tissue to select obs from the full dataset. one of Lymph_Node, Blood, \
                    Lung, Spleen, Thymus, Muscle, Salivary_Gland, Bladder, Fat, Prostate, Vasculature, \
                        Tongue, Large_Intestine, Pancreas, Small_Intestine, Bone_Marrow, Heart, Mammary. \
                            Eye, Kidney, Trachea, Skin, Uterus, Liver")
parser.add_argument('-o', '--output', metavar='output path',  nargs='?',const=o_dir, default=o_dir, type=str, required=False,help="output path")

# parse command-line arguments
args = parser.parse_args()

# create output directory if it does not exist
if not os.path.exists(args.output):
    os.mkdir(args.output)

os.chdir(args.output)

# display arguments
# print("args=%s" % args)
print("args.file=%s" % args.file)
print("args.tissue=%s" % args.tissue)
print("args.output=%s" % args.output)

# import dataset
adata = sc.read_h5ad(fn, backed='r')

# select obs from the dataset by specified tissue
try:
    f = f"TS_{args.tissue}.h5ad"
    print('tissue:', args.tissue, sep=' ')
    adata_s = adata[adata.obs.organ_tissue == args.tissue]
    begin = time.time()
    adata_s.write(os.path.join(args.output, f), compression='gzip')
    end = time.time()
    print(f"Total runtime of the program is {end - begin}.")
    
finally:
    adata.file.close()
