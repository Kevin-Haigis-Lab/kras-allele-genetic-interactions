#!/bin/bash

## Write the libraries used in R and python. ##

module load gcc R/4.0.1 conda2/4.2.13

source ~/.bashrc

# Output files.
R_LIB_FILE="config/R-libraries.txt"
PY_LIB_FILE="config/python-env.yaml"

# Take a snapshot of the library using 'renv'.
Rscript -e "renv::snapshot()"

# Export conda environment.
source activate rctest
conda env export > $PY_LIB_FILE
conda deactivate
