#!/bin/bash
#SBATCH -c 1
#SBATCH --mem 2G
#SBATCH -t 0-00:10
#SBATCH -p short
#SBATCH -o logs/comut-fdr/%A_%a.log
#SBATCH -e logs/comut-fdr/%A_%a.log

module load gcc R/4.0.1
Rscript src/20_82_comutation-fdr-simulation-run.R $1
exit
