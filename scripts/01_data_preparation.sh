#!/bin/bash

# Kill script if any commands fail

set -e
echo "Job Start at `date`"

############################01.QC Generate high-quality data####################################
 # -Note: Take DG1 as an example; Please add the directories of the required software to the environment variables in advance.

mkdir -p ~/reference_genome/Goat_bwa_index
mkdir -p ~/metagenome
cd metagenome
mkdir 00_rawdata 01_cleandata

rawdata=~/metagenome/00_rawdata
cleandata=~/metagenome/01_cleandata
ref=~/reference_genome/Goat_bwa_index/

## QC Step1: filter low quality sequence by fastp

cd ${cleandata}

fastp -i ${rawdata}/DG1_1.raw.fq.gz \
    -o ${cleandata}/DG1_clean_R1.fastq.gz \
    -I ${rawdata}/DG1_2.raw.fq.gz \
    -O ${cleandata}/DG1_clean_R2.fastq.gz \
    --cut_by_quality3 \
    -W 4 -M 20 -n 5 -c -l 50 -w 3

## QC Step2: filt host sequence

### index of reference genome

bwa index -a bwtsw ${ref}/GCF_001704415.2_ARS1.2_genomic.fna

### Mapping to reference

#### Generate .sam file

mkdir -p ~/metagenome/02_QC

bwa mem -t 6 -M -R '@RG\tID:DG1\tLB:DG1\tPL:ILLUMINA\tSM:DG1' \
    ${ref}/GCF_001704415.2_ARS1.2_genomic.fna \
    ${cleandata}/DG1_clean_R1.fastq.gz \
    ${cleandata}/DG1_clean_R2.fastq.gz >~/metagenome/02_QC/DG1.sam

#### Transfer .sam to .bam

samtools view \
    -@ 20 -bS \
    ~/metagenome/02_QC/DG1.sam >~/metagenome/02_QC/DG1.bam

#### Generate DNA contamination rate file

samtools flagstat ~/metagenome/02_QC/DG1.bam >DG1.map.rate

printf DG1 >>~/metagenome/02_QC/DG1.map.rate

#### Extract unmapped sequence and sorted the unmapped sequence .bam file

samtools view -@ 20 -u -f 12 -F 256 -o  ~/metagenome/02_QC/DG1.unmapped.bam ~/metagenome/02_QC/DG1.bam

samtools sort -o ~/metagenome/02_QC/DG1.unmapped.sorted.bam -@ 20 -O bam ~/metagenome/02_QC/DG1.unmapped.bam

#### Change .bam to .fastq

mkdir -p ~/metagenome/03_rm_dna_contamination

bedtools bamtofastq -i ~/metagenome/02_QC/DG1.unmapped.sorted.bam -fq ~/metagenome/03_rm_dna_contamination/DG1_1.hq.fq -fq2 ~/metagenome/03_rm_dna_contamination/DG1_2.hq.fq

#### Return to home directory

cd ~




