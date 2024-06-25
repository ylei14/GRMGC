#!/bin/bash

# Kill script if any commands fail

set -e
echo "Job Start at `date`"

############################03.Construction of the gene catalog####################################

mkdir -p ~/metagenome/05_genecatalog

cd ~/metagenome/05_genecatalog && mkdir before_cdhit cdhit_cluster

# Gene prediction

cd ~/metagenome/05_genecatalog/before_cdhit

prodigal \
    -a ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs.pep \
    -d ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs.fa \
    -i ~/metagenome/04_assembly/Allcontigs/allsample_contigs.fasta \
    -f gff \
    -o ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs.gff \
    -q \
    -p meta

# Remove ORFs shorter than 100bp using a custom script

python del_less_than_100bp.py ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs.fa ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs_filtered.fa

# Construct a non-redundant gene catalog

cd-hit-est \
    -i ~/metagenome/05_genecatalog/before_cdhit/allsample_contigs_filtered.fa \
    -o ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.fa \
    -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 1 -T 80

# Translate the non-redundant gene catalog into a non-redundant protein catalog

transeq -sequence \
    ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.fa \
    -table 11 \
    -trim \
    -outseq ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep 1>transeq.o 2>transeq.e

