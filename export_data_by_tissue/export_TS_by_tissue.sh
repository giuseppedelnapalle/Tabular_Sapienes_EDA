#!/bin/bash
# export the Tabular Sapiens dataset by tissue
#@param file path to h5ad file
#@param tissue tissue or organ in the TS dataset, required
#@param output output path
#@param gene_expression gene expression values from the anndata object
#$author giuseppe
#$date Dec 2022

# py310bio
py3=~/miniconda3/envs/py310bio/bin/python

# path_to_scripts
path_s=~/Project_Data/Python_data/terminal/Tabular_Sapiens/scripts
# py_script_1
s_1=${path_s}/export_TS_subset.py
# py_script_2
s_2=${path_s}/QC_TS_subset.py

while [ ! -z "$1" ]; do
  case "$1" in
     --file|-f)
         shift
         echo "You entered file as: $1"
         f=$1
         ;;
      --tissue|-t)
         shift
         echo "You entered tissue as: $1"
         t=$1
         ;;
      --output|-o)
         shift
         echo "You entered output as: $1"
         o=$1
         ;;
      *)
         echo "invalid argument: $1"
         ;;
  esac
  shift
done

# tissue
if [ -z $t ]; then
  echo "tissue argument not given. program halted."
  exit 1
fi

echo "tissue: $t"

# file
if [ -z $f ]; then
  f=~/Documents/Rsch/resources/HCA/Tabula_Sapiens/TabulaSapiens.h5ad
fi

echo "file: $f"

# output
if [ -z $o ]; then
  o=~/Project_Data/Python_data/terminal/Tabular_Sapiens/output_figshare
fi

echo "output: $o"

# export TS subset as h5ad file
$py3 $s_1 -f $f -t $t -o $o
echo "writing h5ad file for TS subset completed."

# QC & write h5ad file for filtered data
fl=TS_${t}.h5ad
fn=$o/$fl
$py3 $s_2 -f $fn -t $t -o $o
echo "QC & writing h5ad file for filtered TS subset completed."

rm $fn

echo "program finished."
