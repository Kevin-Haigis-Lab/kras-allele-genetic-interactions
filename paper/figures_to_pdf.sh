#!/bin/bash

## A simple script for converting the SVG figures into PDF files.

ORIGINAL_DIR=$(pwd)

cd paper/figures

for img in *svg
do
	echo "Converting $img"
	PDF_NAME="$(basename $img svg)pdf"
	cairosvg $img -o $PDF_NAME -s 1.333
done

cd $ORIGINAL_DIR
exit
