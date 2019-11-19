#!/bin/bash

#SBATCH -c 1
#SBATCH -p priority
#SBATCH -t 1-00:00
#SBATCH --mem 4G
#SBATCH -o logs/rc-test_slurm_logs/snakemake_%A.log
#SBATCH -e logs/rc-test_slurm_logs/snakemake_%A.log


module load gcc conda2 slurm-drmaa/1.1.0
source activate rctest

snakemake \
  --snakefile src/20_20_rc-test-Snakefile \
  --jobs 9950 \
  --restart-times 0 \
  --cluster-config config/rc-test-snakemake-cluster.json \
  --latency-wait 120 \
  --drmaa " -c {cluster.cores} -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -o {cluster.out} -e {cluster.err} -J {cluster.J}"


# to make a dag
# snakemake \
#   --snakefile src/20_20_rc-test-Snakefile \
#   --dag |  \
#   dot -Tpdf > graphs/20_20_rc-test-Snakefile/dag.pdf
