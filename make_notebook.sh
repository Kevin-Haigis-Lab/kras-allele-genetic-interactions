#!/bin/bash

## Copy the graphs from the analysis to the "reports/static/img/"
## to be used in the notebook.

cp -r graphs reports/static/img


## Build hugo website

cd reports/

hugo -d ../docs

cd ..


## Commit to repo

git add -A

git commit -m "build notebook: $(date)"

git push

