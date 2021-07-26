# 2021 systemic sclerosis PBMC and skin punch biopsy RNA-Seq

Preprint: [TBA](https://www.biorxiv.org)
Downloadable data: [FigShare Project](https://figshare.com/projects/2021_Roberson_lab_systemic_sclerosis_transcriptome_data/118698)

## Project description
This project used samples collected from controls without autoimmune disease and inviduals with rheumatologist-confirmed systemic sclerosis (**SSc**). Biological samples included isolated peripheral blood mononuclear cells (**PBMCS**) and skin punch biopsies stabilized to reduce RNA degradation.

We used Takara / Clontech stranded, ribosomal depletion, total RNA library preprations. We sequenced samples with Illumina paired-end sequencing.

## Code description
This repository contains the code used to process data and to generate statistical analysis / paper figures.

Beware all ye who enter here: this project was started before we switched to [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow management. This analysis is instead coordinated using [Make](https://www.gnu.org/software/make/manual/make.html). Analyses for highly parallelizable operations were performed on a cluster at Washington University using PBS scripts. Other operations were performed on a local lab server. Much of the figure generation and analysis was performed using a local laptop. This all means that the overall workflow is less easy to follow and reimplement than a Snakemake workflow, but more reproducible than having no workflow management at all.

## Software Requirements

### Software - general
* R
* Python
* Perl
* salmon
* featureCounts
* RNA-STAR
* pigz

### Annotation and related data
* GRCh38 human genome
* GRCh38 GTF file
* Additional contigs for FASTA file (see preprint or published paper for details; not strictly necessary)
* Additional contig SAF file (see featureCounts documentation)
* Gene start, end, strand SAF file for feature counts (export from Biomart)
* Salmon index
* RNA-STAR index of modified FASTA file

### Software - R packages
* cowplot
* data.table
* DESeq2
* doMC (depending on platform)
* doSNOW (depending on platform)
* dynamicTreeCut
* ggplotify
* ggrepel
* grid
* here
* knitr
* leafcutter
* optparse
* reshape2
* sleuth
* sva
* tidyverse
* UpSetR
* wasabi
* WGCNA

## Recreating analysis
* Clone repository to local environment.

* Install necessary softare, associated packages, and human genome annotation files to your local environment. and packages to your computing environment.

* If you want to include the additional microbial genomes to your FASTA file as we did (detailed in manuscript), concatenate them onto the GRCh38 human reference.

* Generate all necessary (RNA-STAR and salmon) indexes.

* Setup the environmental variables file so the paths work on your system.

* Run PBS scripts 1-8 to do basic alignments and QC. If you aren't running in a PBS / Torque high-performance computing cluster these files could be ported to shell file loops for local execution on a single system.

* Run Makefiles for salmon and leafcutter analysis.

* Knit the R Markdown analysis files in order.

**Critically important note** All command-line analyses should be run from the top-level of the repository, i.e. "make -f makefiles/run_salmon.make" to ensure proper use of relative directories. R markdown files use the "here" package and can be knitted as is.

## Recreating analysis - shortcut method (recommended)
The most important parts of this analysis are performed using R Markdown. If you want to perform the analyses (DESeq2, WGCNA, etc) without running the alignments and such, you can do the R Markdown portion of the analysis. The most critical files will be downloaded automatically from FigShare by running the first script.
