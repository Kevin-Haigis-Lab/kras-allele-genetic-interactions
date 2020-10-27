#!/bin/bash

## Write the libraries used in R and python. ##

module load gcc R/4.0.1 conda2/4.2.13

# Output files.
R_LIB_FILE="config/R-libraries.txt"
PY_LIB_FILE="config/python-env.yaml"

# Export conda environment.
source activate rctest
conda env export > $PY_LIB_FILE
conda deactivate 

# Make a text file listing all R libraries.
Rscript -e "source('lib/make-project-helpers.R');write_package_list('$R_LIB_FILE')"
