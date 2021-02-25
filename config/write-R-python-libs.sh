#!/bin/bash

## Write the libraries used in R and python. ##

module load gcc R/4.0.1 conda2/4.2.13

bash /home/jc604/.bashrc

# Output files.
COMUTATION_ENV_FILE="config/comutation_environment.yaml"
RCTEST_ENV_FILE="config/rctest_environment.yaml"

# Take a snapshot of the library using 'renv'.
Rscript -e "renv::snapshot(force=TRUE)"

# Export comutation conda environment.
conda env export --name comutation --no-builds > "$COMUTATION_ENV_FILE"

# Export RC-test conda environment.
conda env export --name rctest --no-builds > "$RCTEST_ENV_FILE"
