#!/bin/bash
#SBATCH -p priority
#SBATCH -c 1
#SBATCH -t 5
#SBATCH --mem 1G
#SBATCH -o /dev/null
#SBATCH -e /dev/null

# A simple script for converting the SVG figures into PDF files.

module load gcc python/3.7.4

source ~/base-env/bin/activate

ORIGINAL_DIR=$(pwd)

cd paper/figures

for img in *svg
do
	echo "Converting $img"
	PDF_NAME="pdfs/$(basename $img svg)pdf"
	cairosvg $img -o $PDF_NAME -s 1.333
done

deactivate

cd $ORIGINAL_DIR

exit
