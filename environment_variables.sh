#######################
# Project directories #
#######################
COUNT_DIR=output/counts
CLEAN_FASTQ_DIR=clean_fastq
DECOY_COUNT_DIR=output/decoys
JXN_DIR=output/junctions
INSERT_METRIC_DIR=output/insert_metrics
MERGE_BAM_DIR=output/merge_bam
RNASEQ_METRIC_DIR=output/rnaseq_metrics
SJ_TAB_DIR=output/sj_tabs
SORT_BAM_DIR=output/sort_bam
STAR_DIR=output/star
STRINGTIE_OUTPUT=output/stringtie
UNSPLICED_COUNT_DIR=output/unspliced_counts

##############
# Extensions #
##############
FASTQ_EXT=fastq.gz
BAM_EXT=bam

###############
# Sample info #
###############
FASTQ_BASE_NAMES=data/fastq_base_names.txt
RGSM_NAMES=data/rgsm_list.txt
RG_INFO_FILE=data/ssc_rnaseq_readgroup_info.csv

#############
# Junctions #
#############
NEW_SJ_TABS=output/novel_sj_tabs.tsv
SJ_TAB_FILE_LIST=output/sj_tab_files.txt

###################
# Reference files #
###################
DECOY_SAF_NAME=decoys_saf.txt
GENOME_FASTA_NAME=Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
GENOME_FILE_FOR_COVER_BED=Homo_sapiens.GRCh38.dna_sm.primary_assembly.genome_file
GENOME_GTF_NAME=Homo_sapiens.GRCh38.85.gtf
GENOME_SAF_NAME=GRCh38_gene_saf.txt

###################
# Reference paths #
###################
GENOME_DIR=/ssd1/resources/human/GRCh38_r85

DECOY_BED_FILE=${GENOME_DIR}/decoys.bed
DECOY_GENOME_FILE=${GENOME_DIR}/fasta/${GENOME_FILE_FOR_COVER_BED}
DECOY_SAF=${GENOME_DIR}/saf/${DECOY_SAF_NAME}
FASTA=${GENOME_DIR}/fasta/${GENOME_FASTA_NAME}
FASTA_FAI=${FASTA}.fai
GTF=${GENOME_DIR}/gtf/${GENOME_GTF_NAME}
#LEAFCUTTER_DIR=${GENOME_DIR}/leafcutter
SAF=${GENOME_DIR}/saf/${GENOME_SAF_NAME}

##################
# RNASeq metrics #
##################
REF_FLAT=${GENOME_DIR}/refFlat/refFlat.txt.gz
RIBOSOMAL_INTERVALS=${GENOME_DIR}/intervals/ribosomal.interval_list

########################
# Gene annotation file #
########################
GENE_ANNOTATION=biomart_gene_transcript_map.txt.gz

##########
# Salmon #
##########
SALMON_REF=/ssd1/resources/human/GRCh38_r85/salmon
SALMON_GENEMAP=/ssd1/resources/human/GRCh38_r85/cdna/salmon_gene_map.tsv
SALMON_CORES=4
SALMON_BOOTSTRAPS=50

##############
# Leafcutter #
##############
LEAFCUTTER_EXONS=data/leafcutter_exons_GRCh38_r85.txt.gz

#######################
# RNA-STAR parameters #
#######################
# star ref is 2x150.
STAR_REF=/scratch/eli/resource/GRCh38_r85_sm_decoy/star
STAR_OUT_FILTER=0.45
STAR_MULTIMAP_MAX=10

##########
# PICARD #
##########
PICARD=~/src/picard/picard.jar
RECORDS_PER_GB=250000

############################
# Feature Count parameters #
############################
MINMAPQUAL=10

#####################################################################################
# wgcna threads                                                                     #
# this is important, as if you run this on a cluster, make sure your PPN or similar #
# argument is *at least* this value to avoid the job being rejected                 #
#####################################################################################
WGCNA_THREADS=12
