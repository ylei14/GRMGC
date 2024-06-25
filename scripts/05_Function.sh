#!/bin/bash

# Kill script if any commands fail

set -e

echo "Job Start at `date`"


############################05. Functional annotion#####################################
mkdir ~/metagenome/07_function

cd ~/metagenome/07_function

mkdir KEGG eggNOG CAZy CARD VFDB

# KEGG

source ~/mambaforge/bin/activate ~/mambaforge/envs/kofamnead

~/biosoft/kofam_scan-1.3.0/exec_annotation \
    -p ~/Database/KofamKOALA/db/profiles \
    -k ~/Database/KofamKOALA/db/ko_list \
    --cpu 20 \
    -E 1e-5 \
    -f mapper \
    -o ~/metagenome/07_function/KEGG/allsample_nr_geneset_KEGG.txt \
    ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep

# eggNOG

emapper.py \
    -i ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep \
    --temp_dir ~/metagenome/07_function/eggNOG/temp \
    --dmnd_db eggnog_proteins.dmnd \
    --data_dir ~/Database/eggnogdb/ \
    --matrix BLOSUM62 \
    --keep_mapping_files \
    -m diamond --cpu 50 --evalue 1e-5 \
    -o ~/metagenome/07_function/eggNOG/allsample_nr_geneset_eggNOG.txt

# CAZy

hmmscan -o ~/metagenome/07_function/CAZy/allsample_nr_geneset.dbcan \
    --tblout ~/metagenome/07_function/CAZy/allsample_nr_geneset.dbcan.tab \
    -E 1e-5 --cpu 20 \
    ~/Database/dbCAN/dbCAN-HMMdb-V8.txt \
    ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep \

# CARD 

diamond makedb --threads 180 --in protein_fasta_protein_homolog_model.fasta --db ~/Database/CARD

diamond blastp --db ~/Database/CARD.dmnd --query ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep \
    --evalue 1e-5 --query-cover 90 --subject-cover 90 --range-cover 90 --id 60 --out ~/metagenome/CARD/args_anno.txt

# VFDB

diamond makedb --in ~/Database/VFDB/VFDB_setA_pro.fas --db ~/Database/VFDB

diamond blastp --db  ~/metagenome/VFDB.dmnd --query ~/metagenome/05_genecatalog/cdhit_cluster/allsample_nr_geneset.pep --out ~/metagenome/VFDB/vf_anno.txt \
    --threads 100 --evalue 1e-5 --outfmt 6






