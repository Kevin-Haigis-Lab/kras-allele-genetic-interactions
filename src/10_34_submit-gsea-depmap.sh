#!/bin/bash

# Run the GSEA for all files as a job-array

sbatch \
  --array=1-$(ls data/gsea/input/*txt | wc -l)%100 \
  src/10_35_gsea-depmap.sh

exit