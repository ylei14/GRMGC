#!/bin/bash

# Kill script if any commands fail

set -e
echo "Job Start at `date`"



############################04.Taxonomy annotion####################################

# Build the diamond index for the nr database.

mkdir -p ~/metagenome/06_taxonomy/NR_2022_10_12

cd ~/metagenome/06_taxonomy

diamond makedb --threads 180 --in ~/Database/nr.fa --db ~/Database/NR_2022_10_12

# Aligning protein sequence of gene catalog to the NCBI-nr database

diamond blastp \
    -q ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep \
    -d ~/Database/NR_2022_10_12.dmnd \
    --threads 100 --evalue 1e-5 -k 50 --id 30 --sensitive --outfmt 6\
    --out ~/metagenome/06_taxonomy/allsample_nr_geneset_diamond2nr

    # -Note: Match taxonomic information with annotated seqids using taxonkit.


