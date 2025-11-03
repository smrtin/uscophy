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
# uscophy run --genomic input --output output --best_duplicated --genetree --lineage metazoa_odb10 --frag --snakemake "--until busco_download"

uscophy run    \
   --genomic input \
   --output output \
   --lineage metazoa_odb10    \
   --frag    \
   --best_duplicated \
   --alignment_software mafft    \
   --modeltesting TEST    \
   --threads 56 \
   --min_genes 0.5 \
   --category-csv input/samples_metadata.csv

