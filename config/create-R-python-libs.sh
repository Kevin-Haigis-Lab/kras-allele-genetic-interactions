#!/bin/bash

## Install the libraries used in R and create the python virtual env. ##

module load gcc R/4.0.1 conda2/4.2.13

bash ~/.bashrc

# Input files.
COMUTATION_ENV_FILE="config/comutation_environment.yaml"
RCTEST_ENV_FILE="config/rctest_environment.yaml"

# Install the R packages using 'renv'.
Rscript -e "renv::restore()"

# Create the conda virtual environment.
conda env create -f "$COMUTATION_ENV_FILE"
conda env create -f "$RCTEST_ENV_FILE"
