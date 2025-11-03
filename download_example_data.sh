#!/bin/bash

# download reference genomes from NCBI 
# in this example 9 genomes from the Balanomorpha, which are an order of Barnacles

echo "download genomes"

#datasets \
#	download \
#	genome \
#	taxon Balanomorpha \
#	--filename Example_Dataset.zip \
#	--reference

# we can use the Taxonomy name to download genomes or provide a list with accession numbers
#datasets \
#	download \
#	genome \
#	taxon  \
#	--inputfile outgroup_list.txt \
#	--filename Outgroup_Dataset.zip \
#	--reference


# extract the genomes form the datasets zip

echo "extract genomes"

python uscophy/workflow/scripts/extract_dataset.py --generate-category --output-dir input Example_Dataset.zip
