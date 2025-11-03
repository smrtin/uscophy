#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -j n
#$ -N uscophy_main
#$ -q fast.q
#$ -pe smp 56

module load miniforge/24.7.1-0
conda activate ./test_environment

# if compute nodes are not connected to internet, download busco lineage on the headnode
# uscophy run --genomic input_test --output output --best_duplicated --genetree --lineage metazoa_odb10 --frag --snakemake "--until busco_download"

uscophy run    \
   --genomic input_test \
   --output output\
   --best_duplicated \
   --genetree \
   --lineage metazoa_odb10    \
   --frag    \
   --alignment_software mafft    \
   --modeltesting TEST    \
   --threads ${NSLOTS} \
   --min_genes 0.8 \
   --snakemake "--latency-wait 30  --keep-incomplete"

