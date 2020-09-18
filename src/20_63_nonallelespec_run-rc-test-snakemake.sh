#!/bin/bash

#SBATCH -c 1
#SBATCH -p priority
#SBATCH -t 2-00:00
#SBATCH --mem 4G
#SBATCH -o logs/rc-nonallelespec/snakemake_%A.log
#SBATCH -e logs/rc-nonallelespec/snakemake_%A.log

module unload python
module load gcc conda2 slurm-drmaa/1.1.1
conda activate rctest


snakemake \
  --snakefile src/20_62_nonallelespec_rc-test-Snakefile.py \
  --jobs 9980 \
  --restart-times 0 \
  --cluster-config config/rc-test-snakemake-cluster.json \
  --latency-wait 120 \
  --drmaa " -c {cluster.cores} -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -o {cluster.out} -e {cluster.err} -J {cluster.J}"

conda deactivate
