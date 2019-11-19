#!/bin/bash

#SBATCH -p priority
#SBATCH -c 1
#SBATCH -t 0-6:00
#SBATCH --mem 8G
#SBATCH --mail-type NONE
#SBATCH -o logs/annovar/24_annovar-annotation-mutations.log
#SBATCH -e annovar_run.log

## INPUT (positional)
#   1: path to input file
#   2: path for output file

module load gcc annovar/20170601

ANNOVAR_PATH=/n/app/annovar/20170601

$ANNOVAR_PATH/table_annovar.pl \
	$1 \
	$ANNOVAR_PATH/humandb \
	--protocol refGene,dbnsfp30a,clinvar_20160302,cosmic68wgs,icgc21 \
	--operation g,f,f,f,f \
	--verbose \
	--buildver hg19 \
	--remove \
	--nastring NA \
	--otherinfo \
	-out $2

echo "ANNOVAR FINISHED"
