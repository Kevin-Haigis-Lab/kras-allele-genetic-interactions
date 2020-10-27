#!/bin/bash

## Install the libraries used in R and create the python virtual env. ##

module load gcc R/4.0.1 conda2/4.2.13

# Input files.
R_LIB_FILE="config/R-libraries.txt"
PY_LIB_FILE="config/python-env.yaml"

# Install the R pakcages.
Rscript -e "source('lib/make-project-helpers.R');install_pacakges($R_LIB_FILE)"

# Create the conda virtual environment.
conda env create -f $PY_LIB_FILE