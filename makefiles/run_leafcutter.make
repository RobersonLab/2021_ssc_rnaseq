#########
# Setup #
#########
SHELL := /bin/bash
BASE := $(shell cat data/fastq_base_names.txt)
RGSM := $(shell cat data/rgsm_list.txt)
TISSUES := PBMC Skin

.SECONDARY:
	#pass
	
###############################
# pull from environment files #
###############################
LEAFCUTTER_JUNCTION_DIR := $(shell cat environment_variables.sh | grep JXN_DIR | awk -F= '{print $$2}')
LEAFCUTTER_CORES := $(shell cat environment_variables.sh | grep SALMON_CORES | awk -F= '{print $$2}')
LEAFCUTTER_EXONS := $(shell cat environment_variables.sh | grep LEAFCUTTER_EXONS | awk -F= '{print $$2}')

###############
# directories #
###############
CLEAN_FASTQ := clean_fastq
DATA := data
LOGS := logs
OUTPUT := output
PBS := pbs
RESULTS := results

TISSUE_INTRON_CLUSTER_DIR := $(OUTPUT)/leafcutter_intron_clusters
TISSUE_DIFF_SPLICING_DIR := $(OUTPUT)/leafcutter_diff_splicing

DIRS := $(LOGS) $(OUTPUT) $(LEAFCUTTER_JUNCTION_DIR) $(TISSUE_INTRON_CLUSTER_DIR) $(TISSUE_DIFF_SPLICING_DIR)

##############
# leafcutter #
##############
SAMPLE_JXNS := $(addsuffix .junc, $(addprefix $(LEAFCUTTER_JUNCTION_DIR)/, $(RGSM)))
CLUSTER_LISTS := $(addsuffix _perind_numers.counts.gz, $(addprefix $(TISSUE_INTRON_CLUSTER_DIR)/, $(TISSUES)))
TISSUE_DIFF_SPLICING := $(addsuffix _significance.txt, $(addprefix $(TISSUE_DIFF_SPLICING_DIR)/, $(TISSUES)))

#######
# run #
#######
all: $(DIRS) $(SAMPLE_JXNS) $(CLUSTER_LISTS) $(TISSUE_DIFF_SPLICING)

################################
# make any missing directories #
################################
$(DIRS):
	mkdir -p $@

########################
# leafcutter junctions #
########################
$(LEAFCUTTER_JUNCTION_DIR)/%.junc: $(OUTPUT)/sort_bam/%_rgid.bam
	bash src/leafcutter_scripts/scripts/bam2junc.sh $^ $@ 1>$(LOGS)/$*_junc.log

##############################
# leafcutter intron clusters #
##############################
$(TISSUE_INTRON_CLUSTER_DIR)/%_perind_numers.counts.gz: $(SAMPLE_JXNS)
	python src/leafcutter_scripts/scripts/leafcutter_cluster.py \
	-j $(DATA)/leafcutter_$*_0m_intron_cluster.txt \
	-m 30 \
	-l 500000 \
	-o $(TISSUE_INTRON_CLUSTER_DIR)/$* 1>$(LOGS)/$*_leafcutter_intron_clustering.log 2>&1 && \
	mv *_$*_*.junc.*sorted.gz $(TISSUE_INTRON_CLUSTER_DIR)
	
####################################
# leafcutter differential splicing #
####################################
$(TISSUE_DIFF_SPLICING_DIR)/%_significance.txt: $(TISSUE_INTRON_CLUSTER_DIR)/%_perind_numers.counts.gz $(DATA)/leafcutter_%_0m_contrast.txt
	Rscript src/leafcutter_scripts/scripts/leafcutter_ds.R \
	--num_threads $(LEAFCUTTER_CORES) \
	--exon_file $(LEAFCUTTER_EXONS) \
	--min_samples_per_group 4 \
	--min_coverage 20 \
	-o $(TISSUE_DIFF_SPLICING_DIR)/$* \
	$^
	