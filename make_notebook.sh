#!/bin/bash

# Run to build, commit, and push notebook

hugo -s reports -d docs --cleanDestinationDir

git add -A

git commit -m "build notebook: $(date)"

git push

