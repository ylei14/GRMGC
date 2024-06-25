#!/bin/bash

# Kill script if any commands fail

set -e
echo "Job Start at `date`"

##############06 Binning: Contigs to MAGs###########################

mkdir ~/metagenome/Binning

cd ~/metagenome/Binning && mkdir 01_Binning 02_dRep_99 02_dRep_95 03_Taxonomy 04_Genome_annotation 05_MAG_quant 06_Phylophlan

cd ~/metagenome/Binning/01_Binning

mkdir INITIAL_BINNING BIN_REFINEMENT CHECKM

# Perform single-sample binning using three different methods

mamba activate metawrap-env

metawrap binning -t 6 \
    -a ~/metagenome/04_assembly/Allcontigs/DG1.contigs.fa \
    -o ~/metagenome/Binning/INITIAL_BINNING/DG1/ \
    -l 500 \
    --universal \
    --metabat2 \
    --maxbin2 \
    --concoct \
    --interleaved \
    ~/metagenome/03_rm_dna_contamination/DG1_1.hq.fq \
    ~/metagenome/03_rm_dna_contamination/DG1_2.hq.fq  

# Consolidate bin sets with the Bin_refinement module

metawrap bin_refinement -t 40 -m 500 \
    -o ~/metagenome/Binning/BIN_REFINEMENT/DG1/ \
    -A ~/metagenome/Binning/INITIAL_BINNING/DG1/metabat2_bins/ \
    -B ~/metagenome/Binning/INITIAL_BINNING/DG1/maxbin2_bins/ \
    -C ~/metagenome/Binning/INITIAL_BINNING/DG1/concoct_bins/ \
    -c 50 -x 10 \
    --quick

# Perform quality assessment using CheckM

checkm lineage_wf \
    -t 80 \
    -x fa \
    --tab_table \
    -f ~/metagenome/Binning/CHECKM/all_mags_checkm.tsv \
    ~/metagenome/Binning/BIN_REFINEMENT/*/*.fa \
    ~/metagenome/Binning/CHECKM/

# Collect all MAGs

cd ~/metagenome/Binning/01_Binning

mkdir Allmags

mv ~/metagenome/Binning/BIN_REFINEMENT/*/*.fa ~/metagenome/Binning/01_Binning/Allmags

# Generate strain-level MAGs

source ~/mambaforge/bin/activate /ifs1/User/leiyu/mambaforge/envs/drep-env 

dRep dereplicate \
    ~/metagenome/Binning/02_dRep_99 \
    -g ~/metagenome/Binning/01_Binning/Allmags/*.fa \
    -p 128 \
    -d \
    --ignoreGenomeQuality \
    -pa 0.95 \
    -sa 0.99 \
    -cm larger

# Generate species-level genome bins (SGBs)

dRep dereplicate \
    ~/metagenome/Binning/02_dRep_95 \
    -g ~/metagenome/Binning/01_Binning/Allmags/*.fa \
    -p 128 \
    -d \
    --ignoreGenomeQuality \
    -pa 0.9 \
    -sa 0.95 \
    -cm larger

# Taxonomic classification of strain-level MAGs and SGBs

source ~/mambaforge/bin/activate /ifs1/User/leiyu/mambaforge/envs/GTDB-tk

cd ~/metagenome/Binning/03_Taxonomy

mkdir GTDB_99 GTDB_95

gtdbtk classify_wf \
    --cpus 40 \
    --out_dir ~/metagenome/Binning/03_Taxonomy/GTDB_99 \
    --genome_dir ~/metagenome/Binning/02_dRep_99/dereplicated_genomes \
    --extension fa \
    --skip_ani_screen

gtdbtk classify_wf \
    --cpus 40 \
    --out_dir ~/metagenome/Binning/03_Taxonomy/GTDB_95 \
    --genome_dir ~/metagenome/Binning/02_dRep_95/dereplicated_genomes \
    --extension fa \
    --skip_ani_screen

# Calling ORFs if MAGs

    # -Note: use a MAG as an example

prodigal -i ~/metagenome/Binning/01_Binning/Allmags/DG1_bin.13.fa \
    -a ~/metagenome/Binning/04_Genome_annotation/DG1_bin.13.pep \
    -d ~/metagenome/Binning/04_Genome_annotation/DG1_bin.13.cds \
    -f gff -g 11 \
    -o ~/metagenome/Binning/04_Genome_annotation/DG1_bin.13.gff \
    -p single

    # -Note: Select the result files to be retained and delete the unnecessary result files.

# Phylogenetic analysis of MAGs

source ~/mambaforge/bin/activate /ifs1/User/leiyu/mambaforge/envs/phylophlan

phylophlan \
    -i ~/metagenome/Binning/04_Genome_annotation/ \  # Save the folder containing MAG .pep files
    -t a \
    -d phylophlan \
    --diversity medium \
    -f supermatrix_aa.cfg \
    --nproc 60 \
    -o ~/metagenome/Binning/06_Phylophlan \
    --min_num_markers 50 \
    --proteome_extension .pep

# Calculate the relative abundance (TPM) of MAGs in each sample

source ~/mambaforge/bin/activate /ifs1/User/leiyu/mambaforge/envs/metawrap-env

metawrap quant_bins \
    -b ~/metagenome/Binning/02_dRep_99/dereplicated_genomes \
    -o ~/metagenome/Binning/05_MAG_quant \
    -t 64 \
    ~/metagenome/03_rm_dna_contamination/*.fq








