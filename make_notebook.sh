#!/bin/bash

## Copy the graphs from the analysis to the "reports/static/img/"
## to be used in the notebook.

REPORTS_IMG_DIR="reports/static/img"
[ ! -d $REPORTS_IMG_DIR ] && mkdir $REPORTS_IMG_DIR

cp -r graphs $REPORTS_IMG_DIR


## Copy the subversioned figures to "reports/static/img/figures/" 
FIGURE_DIR="paper/figures"
DEST_FIG_DIR="${REPORTS_IMG_DIR}/figures"
[ ! -d $DEST_FIG_DIR ] && mkdir $DEST_FIG_DIR

for f in $( find $FIGURE_DIR -name "*Figure*svg" ); do
	cp $f $DEST_FIG_DIR
done


## Copy the full figures to "reports/content/home/gallery/gallery/"
GALLERY_DIR="reports/content/home/gallery/gallery"
for f in $( find $FIGURE_DIR -d 1 -name "*Figure*svg" ); do
	cp $f $GALLERY_DIR
done


## Build hugo website

cd reports/
hugo -d ../docs
cd ..


## Commit to repo

git add -A
git commit -m "build notebook: $(date)"
git push
