#!/bin/bash

# Run to build, commit, and push notebook

cd reports/

hugo -d ../docs

cd ..

git add -A

git commit -m "build notebook: $(date)"

git push

