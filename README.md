# GRMGC
Construction of the microbial reference gene catalog and metagenomic assembled genome database for goat rumen microbiome.
# INSTALLATION
Before running, please ensure that the appropriate software and databases are installed.
## Software installation
|Software   |Availability   |
|--|--|
|fastp(v0.20.0)   |https://github.com/OpenGene/fastp   |
|bwa(v0.7.17)   |https://github.com/lh3/bwa   |
|Samtools(v1.17)   |https://github.com/samtools/samtools/releases/   |
|metaSPAdes(v3.15.5)   |https://github.com/ablab/spades   |
|Bowtie2(v2.5.1)   |https://github.com/BenLangmead/bowtie2   |
|MEGAHIT(v1.2.9)   |https://github.com/voutcn/megahit  |
|Prodigal(v2.6.3)   |https://github.com/hyattpd/Prodigal   |
|CD-HIT(v4.8.1)   |https://github.com/weizhongli/cdhit   |
|EMBOSS Transeq(v6.6.0.0)   |http://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html  |
|SOAPaligner(v2.21)   |https://github.com/ShujiaHuang/SOAPaligner   |
|DIAMOND(v2.1.8.162)   |https://github.com/bbuchfink/diamond/   |
|eggnog-mapper(v2.1.11)   |http://eggnog5.embl.de/   |
|KOfamScan(v1.1.0)   |https://github.com/takaram/kofam_scan   |
|eggnog-mapper(v2.1.11)   |http://eggnog5.embl.de/   |
|metawrap(v1.3.2)   |https://github.com/bxlab/metaWRAP   |
|dRep(v3.2.2)   |https://github.com/MrOlm/drep   |
|GTDB-tk(2.3.0)   |http://gtdb.ecogenomic.org/   |
|PhyloPhlAn(v3.0.67)   |https://github.com/biobakery/phylophlan   |
|BLASTP(v2.14.1+)   |https://blast.ncbi.nlm.nih.gov/Blast.cgi   |

Note: Some software is easier to use when configured in a virtual environment, such as MetaWRAP and KofamScan.
## Database installation
|Database   |Version   |Description   |Availability   |
|--|--|--|--|
|goat   |Genome assembly ARS1.2   |Goat reference genome   |https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001704415.2/|
|NCBI-nr   | 20221012  |Protein database   |https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/   |
|EggNOG   |EggNOG5.0   |EggNOG annotation	|http://eggnog5.embl.de/#/app/downloads  |
|CARD   |v3.0.9   |CARD annotion   |https://card.mcmaster.ca   |   
|VFDB   |20200703   |VFDB annotion   |http://www.mgc.ac.cn/VFs/main.htm   |
|GTDB-tk   |release214    |Taxonomic assignments of MAGs   |https://gtdb.ecogenomic.org/downloads   |

Note: Update the database versions if necessary.

# Overview of the scripts
## 01_data_preparation.sh
Metagenomic raw data processing: reads trimming and DNA contamination removal, generating high-quality data.
## 02_Assembly.sh
Assembly: assemble short reads into contigs.
## 03_Genecatalog.sh
Construct non-redundant gene catalog: predict ORFs and remove redundancy to build the non-redundant gene catalog.
## 04_Taxonomy.sh
Taxonomic annotation: align protein sequences from the non-redundant gene catalog to the NCBI-nr database to obtain taxonomic annotation information.
## 05_Function.sh
Functional annotation: align protein sequences from the non-redundant gene catalog to the eggNOG, KEGG, CAZy, CARD, and VFDB databases to obtain functional annotation 
information.
## 06_Genome_reconstraction.sh
The steps related to the reconstruction of metagenome-assembled genomes (MAGs) are outlined in this section. Binning, refining, reassembling, annotating genomes, and estimating MAG abundances were carried out using the modules of the metaWRAP pipeline. MAG dereplication was performed using dRep software. Taxonomic classification and phylogenetic analysis were conducted with GTDB-tk and PhyloPhlAn, respectively.

Note: many custom scripts were involved in the analysis process. Please contact me for access.
