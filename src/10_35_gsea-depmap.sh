#!/bin/bash

#SBATCH -p priority
#SBATCH --mem 6G
#SBATCH -c 2
#SBATCH -t 0-5:00
#SBATCH --mail-type NONE
#SBATCH -o logs/gsea-slurm/gsea_%A.log
#SBATCH -e logs/gsea-slurm/gsea_%A.log
#SBATCH --x11=batch

# Run GSEA for DepMap data
# This script runs a single input file as a part of a job array.

module load gcc java

rm -r data/gsea/output/*

GSEA_PATH=/home/jc604/mysoftware/gsea-3.0.jar
GENE_SETS=data/gsea/genesets/h.all.v7.0.symbols.gmt
OUT_DIR=data/gsea/output

for EXPRS_FILE in $(ls data/gsea/input/*txt); do

	CLS_FILE=data/gsea/input/$(basename $EXPRS_FILE .txt).cls
	RESULTS_NAME=$(basename $EXPRS_FILE .txt)
	
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
		-plot_top_x 20 \
		-rnd_seed timestamp \
		-set_max 500 \
		-set_min 15 \
		-zip_report false \
		-out $OUT_DIR \
		-gui false
done

exit
