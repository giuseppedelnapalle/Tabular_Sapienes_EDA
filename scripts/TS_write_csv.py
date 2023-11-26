#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
exporting subset of the Tabular Sapiens dataset as csv files
$input TS_tissue.h5ad created by export_TS_subset.py
$output csv files for the Tabular Sapiens dataset
# format csv
"""

import os
import scanpy as sc
import argparse
import time

# set working directory
wd = '/home/nikola/Project_Data/Python_data/terminal/Tabular_Sapiens'
os.chdir(wd)

# output directory
o_dir = os.path.join(wd, 'output_figshare')

# argument parser
parser = argparse.ArgumentParser()

parser.add_argument('-f', '--file', metavar='h5ad file', type=str, required=False, help="path to h5ad file as input")
parser.add_argument('-t', '--tissue', metavar='tissue or organ', type=str, required=True, 
                    help="tissue to select obs from the full dataset. one of Lymph_Node, Blood, \
                    Lung, Spleen, Thymus, Muscle, Salivary_Gland, Bladder, Fat, Prostate, Vasculature, \
                        Tongue, Large_Intestine, Pancreas, Small_Intestine, Bone_Marrow, Heart, Mammary. \
                            Eye, Kidney, Trachea, Skin, Uterus, Liver")
parser.add_argument('-g', '--gene_expression', metavar='gene expression values', nargs='?', const='raw.X', default='raw.X', type=str, required=False, 
                   help="gene expression values to export, either raw.X or X")
parser.add_argument('-o', '--output', metavar='output path',  nargs='?', const=o_dir, default=o_dir, type=str, required=False,help="output path")

# parse command-line arguments
args = parser.parse_args()

# input file
fl = f"TS_{args.tissue}.h5ad"
fn = '/'.join((o_dir, fl))

args.file = args.file if args.file else fn

# display arguments
# print("args=%s" % args)
print("args.file=%s" % args.file)
print("args.tissue=%s" % args.tissue)
print("args.gene_expression=%s" % args.gene_expression)
print("args.output=%s" % args.output)

# import dataset
adata = sc.read_h5ad(args.file)

# write csv files
print('tissue:', args.tissue, sep=' ')

if args.gene_expression=='raw.X':
    f = f"{args.tissue}_raw_X"
    begin = time.time()
    adata.raw.to_adata().write_csvs(os.path.join(args.output, f), skip_data=False)
    end = time.time()
    print(f"Total runtime of the program is {end - begin}")
else:
    f = f"{args.tissue}_X"
    begin = time.time()
    adata.write_csvs(os.path.join(args.output, f), skip_data=False)
    end = time.time()
    print(f"Total runtime of the program is {end - begin}")
    

