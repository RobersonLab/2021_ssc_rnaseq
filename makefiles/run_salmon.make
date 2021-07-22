#########
# Setup #
#########
SHELL := /bin/bash
BASE := $(shell cat data/fastq_base_names.txt)
RGSM := $(shell cat data/rgsm_list.txt)

.SECONDARY:
	#pass
	
###############################
# pull from environment files #
###############################
SALMON_REF := $(shell cat environment_variables.sh | grep SALMON_REF | awk -F= '{print $$2}')
SALMON_GENEMAP := $(shell cat environment_variables.sh | grep SALMON_GENEMAP | awk -F= '{print $$2}')
SALMON_CORES := $(shell cat environment_variables.sh | grep SALMON_CORES | awk -F= '{print $$2}')
SALMON_BOOTSTRAPS := $(shell cat environment_variables.sh | grep SALMON_BOOTSTRAPS | awk -F= '{print $$2}')

###############
# directories #
###############
CLEAN_FASTQ := clean_fastq
DATA := data
LOGS := logs
OUTPUT := output
PBS := pbs
RESULTS := results

OUTPUT_SUBDIRS := count_gene count_gene_region multiqc order rnaseq_qc splice_juncs salmon star
OUTPUT_SUBDIRS := $(addprefix $(OUTPUT)/, $(OUTPUT_SUBDIRS))

DIRS := $(RESULTS) $(CLEAN_FASTQ) $(DATA) $(LOGS) $(OUTPUT) $(OUTPUT_SUBDIRS) $(PBS)
	
#########################
# confirm order correct #
#########################
ORDER := $(addprefix $(OUTPUT)/order/, $(addsuffix _order.txt, $(BASE)))

#########################
# Salmon quantification #
#########################
SALMON_QUANT := $(addsuffix /quant.sf, $(addprefix $(OUTPUT)/salmon/, $(RGSM)))
H5_ABUNDANCE := $(addsuffix /abundance.h5, $(addprefix $(OUTPUT)/salmon/, $(RGSM)))
TISSUES := PBMC Skin
DET := $(addprefix results/sleuth_, $(addsuffix _SSc_vs_Control_transcript_id.csv, $(TISSUES)))
DEG_DET := $(addprefix results/sleuth_, $(addsuffix _SSc_vs_Control_onlyDeseq2Genes.csv, $(TISSUES)))

#######
# run #
#######
all: $(DIRS) $(ORDER) $(SALMON_QUANT) $(H5_ABUNDANCE) $(DET) $(DEG_DET)

################################
# make any missing directories #
################################
$(DIRS):
	mkdir -p $@

#########################
# compile compare order #
#########################
bin/compare_order:
	g++ -m64 -O3 -o bin/compare_order src/compare_order/*.cpp 

#####################
# check FASTQ order #
#####################
$(OUTPUT)/order/%_order.txt: bin/compare_order $(CLEAN_FASTQ)/%_1_clean.fastq.gz $(CLEAN_FASTQ)/%_2_clean.fastq.gz
	bin/compare_order --format fastq --file1 <(gunzip -c $(CLEAN_FASTQ)/$*_1_clean.fastq.gz) --file2 <(gunzip -c $(CLEAN_FASTQ)/$*_2_clean.fastq.gz) 1>$@ 2>&1
	
################
# salmon quant #
################
$(OUTPUT)/salmon/%/quant.sf:
	$(eval LOCAL_IN_STRING=$(shell python3 bin/get_salmon_input_string.py --fastqpath $(CLEAN_FASTQ) $* data/ssc_rnaseq_readgroup_info.csv))

	salmon quant \
	--index $(SALMON_REF) \
	--libType ISF \
	--seqBias \
	--gcBias \
	--validateMappings \
	--geneMap $(SALMON_GENEMAP) \
	$(LOCAL_IN_STRING) \
	--threads $(SALMON_CORES) \
	--numBootstraps $(SALMON_BOOTSTRAPS) \
	--output $(OUTPUT)/salmon/$* \
	1>logs/$*_salmon_run.log \
	2>&1

#####################
# wasabi for sleuth #
#####################
$(OUTPUT)/salmon/%/abundance.h5: $(OUTPUT)/salmon/%/quant.sf
	R --vanilla < src/prep_for_sleuth.R $* 1>logs/$*_wasabi.log 2>&1

#################
# transcript DE #
#################
results/sleuth_%_SSc_vs_Control_transcript_id.csv: src/sleuth_%_det.R
	R --vanilla < $^ 1>logs/sleuth_$*_transcript_de.log 2>&1

###################################
# transcripts only diff exp genes #
###################################
results/sleuth_%_SSc_vs_Control_onlyDeseq2Genes.csv: src/sleuth_%_degfilter_det.R
	R --vanilla < $^ 1>logs/sleuth_$*_deg_transcript_de.log 2>&1
