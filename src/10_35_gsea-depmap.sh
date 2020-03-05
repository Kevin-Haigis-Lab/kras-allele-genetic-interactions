#!/bin/bash

#SBATCH -p short
#SBATCH --mem 6G
#SBATCH -c 2
#SBATCH -t 0-4:30
#SBATCH --mail-type NONE
#SBATCH -o logs/gsea-slurm/gsea_%A_%a.log
#SBATCH -e logs/gsea-slurm/gsea_%A_%a.log
#SBATCH --x11=batch

# Run GSEA for DepMap data
# This script runs a single input file as a part of a job array.

module load gcc java

# Paths to shared items: GSEA jar, gene sets, output directory.
GSEA_PATH=/home/jc604/mysoftware/gsea-3.0.jar
GENE_SETS='data/gsea/genesets/h.all.v7.0.symbols.gmt,data/gsea/genesets/c2.all.v7.0.symbols.gmt'
OUT_DIR=data/gsea/output

# The expression file based on the array index.
EXPRS_FILE=$(ls data/gsea/input/*txt | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Get the CLS file and result name from the expression file.
CLS_FILE=data/gsea/input/$(basename $EXPRS_FILE .txt).cls
RESULTS_NAME=$(basename $EXPRS_FILE .txt)

# Print out to log.
echo $EXPRS_FILE
echo $CLS_FILE
echo $RESULTS_NAME

# Run GSEA
java -cp $GSEA_PATH \
	-Xmx5000m \
	xtools.gsea.Gsea \
	-gmx $GENE_SETS \
	-res $EXPRS_FILE \
	-cls $CLS_FILE \
	-collapse false \
	-norm meandiv \
	-nperm 10000 \
	-permute gene \
	-scoring_scheme weighted \
	-rpt_label $RESULTS_NAME \
	-create_svgs true \
	-make_sets true \
	-plot_top_x 100 \
	-rnd_seed 23 \
	-set_max 500 \
	-set_min 15 \
	-zip_report false \
	-out $OUT_DIR \
	-gui false \
	-plot_top_x 50

exit
