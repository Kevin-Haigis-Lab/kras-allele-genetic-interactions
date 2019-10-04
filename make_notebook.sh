#!/bin/bash

# Run to build, commit, and push notebook

cd reports/

hugo ../docs

cd ..

git add -A

git commit -m "build notebook: $(date)"

git push

