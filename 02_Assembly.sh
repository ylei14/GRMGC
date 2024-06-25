
#!/bin/bash

# Kill script if any commands fail

set -e
echo "Job Start at `date`"

############################02.Assembly: High-quelity data to contigs####################################
# Note: Take DG1 as an example; Please add the directories of the required software to the environment variables in advance.

hqdata=~/metagenome/03_rm_dna_contamination
mkdir -p ~/metagenome/04_assembly/DG1_spades

# Single sample assembly

cd ~/metagenome/04_assembly/DG1_spades

metaspades.py \
    -1 ${hqdata}/DG1_1.hq.fq\
    -2 ${hqdata}/DG1_2.hq.fq\
    --only-assembler \
    -m 400 \
    -t 40 \
    -o ~/metagenome/04_assembly/DG1_spades

    # -Note: remove contigs shorter than 500bp using the custom script.
    
## change the contig file name

mv ~/metagenome/04_assembly/DG1_spades/contigs.fasta ~/metagenome/04_assembly/DG1_spades/DG1.contigs.fa

# Acquire the unassembled reads

mkdir -p ~/metagenome/04_assembly/mix_unassemble/

## Build index

bowtie2-build --threads 20 ~/metagenome/04_assembly/DG1_spades/DG1.contigs.fa ~/metagenome/04_assembly/mix_unassemble/DG1.contigs

## Generate .sam file

bowtie2 --threads 20 -x \
    ~/metagenome/04_assembly/mix_unassemble/DG1.contigs \
    -1 ${hqdata}/DG1_1.hq.fq \
    -2 ${hqdata}/DG1_2.hq.fq \
    -S ~/metagenome/04_assembly/mix_unassemble/DG1.sam

## Extract unmapped sequence and sorted the unmapped sequence .bam file

samtools view -@ 20 -b -f 12 -F 256 ~/metagenome/04_assembly/mix_unassemble/DG1.sam >~/metagenome/04_assembly/mix_unassemble/DG1_unmapped.bam

samtools sort -@ 20 -n ~/metagenome/04_assembly/mix_unassemble/DG1_unmapper.bam -o ~/metagenome/04_assembly/mix_unassemble/DG1_unmapped_sorted.bam

## Transfer .bam to fastq 

bedtools bamtofastq -i \
    ~/metagenome/04_assembly/mix_unassemble/DG1_unmapped_sorted.bam \
    -fq ~/metagenome/04_assembly/mix_unassemble/DG1_unass_r1.fastq \
    -fq2 ~/metagenome/04_assembly/mix_unassemble/DG1_unass_r2.fastq

    # -Note: Delete the .sam and .bam files generated during the analysis if necessary. 

# Merge all unmapped reads

cat ~/metagenome/04_assembly/mix_unassemble/*_unass_r1.fastq >~/metagenome/04_assembly/mix_unassemble/unmapped_r1.fastq
cat ~/metagenome/04_assembly/mix_unassemble/*_unass_r2.fastq >~/metagenome/04_assembly/mix_unassemble/unmapped_r2.fastq

# All unmapped reads were co-assembled

megahit -1 ~/metagenome/04_assembly/mix_unassemble/unmapped_r1.fastq \
        -2 ~/metagenome/04_assembly/mix_unassemble/unmapped_r2.fastq \
        --min-count 2 \
        --k-min 21 \
        --k-max 141 \
        --k-step 10 \
        -t 40 \
        -m 500 \
        --min-contig-len 500 \
        -o ~/metagenome/04_assembly/mix_unassemble/

mv ~/metagenome/04_assembly/mix_unassemble/final.contigs.fa ~/metagenome/04_assembly/mix_unassemble/unmapped.contigs.fa

# Merge all contigs from single sample assembled and co-assembled from unmapped reads

mkdir -p ~/metagenome/04_assembly/Allcontigs

mv ~/metagenome/04_assembly/*_spades/*.contigs.fa ~/metagenome/04_assembly/Allcontigs

mv ~/metagenome/04_assembly/mix_unassemble/unmapped.contigs.fa ~/metagenome/04_assembly/Allcontigs

cat ~/metagenome/04_assembly/Allcontigs/*.contigs.fa >allsample_contigs.fasta

