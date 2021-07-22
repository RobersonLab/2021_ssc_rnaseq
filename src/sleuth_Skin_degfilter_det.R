#############
# libraries #
#############
library( tidyverse )
library( here )
library( sleuth )

#############
# constants #
#############
tissue <- "Skin"
case_name <- "SSc"
control_name <- "Control"

#############
# variables #
#############
env_variables <- read_tsv( file = here::here( 'environment_variables.sh' ), col_names = FALSE, col_types = 'c', comment = '#' ) %>%
	separate( data = ., col = 'X1', sep = '=', into = c( 'name', 'value' ) )

salmon_bootstraps <- filter( .data = env_variables, name == 'SALMON_BOOTSTRAPS' ) %>%
	pull( .data = ., value ) %>%
	as.integer( . )
	
salmon_gene_transcript_map_file <- filter( .data = env_variables, name == 'SALMON_GENEMAP' ) %>%
	pull( .data = ., value ) %>%
	as.character( . )
	
gene_annotation_file <- filter( .data = env_variables, name == "GENE_ANNOTATION" ) %>%
	pull( .data = ., value ) %>%
	here::here( 'data', . )

###############
# load design #
###############
read_groups <- read_csv( file = here::here( "data","ssc_rnaseq_readgroup_info.csv" ), col_names = TRUE )

#######################
# load annotation map #
#######################
salmon_gene_transcript_map <- read_tsv( file = salmon_gene_transcript_map_file, col_names = c( 'transcript_id', 'gene_id' ) ) %>%
	mutate( .data = ., gene_id_with_version = gene_id ) %>%
	mutate( .data = ., gene_id = str_replace( string = gene_id, pattern = "\\.[0-9]+$", replacement = "" ) )
	
transcript_ids_to_test <- read_tsv( file = here::here( 'results', paste0( 'SSc_0m_DiffExp_', tissue, '_gene.tsv' ) ), col_names = TRUE ) %>%
	filter( qval < 0.05 ) %>%
	select( gene_id ) %>%
	merge( x = salmon_gene_transcript_map, y = ., by = 'gene_id' ) %>%
	dplyr::pull( transcript_id )

full_annotation <- read_tsv( file = gene_annotation_file, col_names = TRUE, col_types = cols( .default = 'c' ) ) %>%
	rename( .data = ., gene_id = `Gene stable ID` ) %>%
	rename( .data = ., transcript_id = `Transcript stable ID` ) %>%
	rename( .data = ., symbol = `Gene name` ) %>%
	rename( .data = ., gene_biotype = `Gene type` ) %>%
	rename( .data = ., transcript_biotype = `Transcript type` )
	
gene_annotation <- full_annotation %>%
	select( .data = ., gene_id, symbol ) %>%
	as.data.frame( . ) %>%
	unique( . )
	
transcript_annotation <- full_annotation %>%
	select( .data = ., transcript_id, gene_id, symbol ) %>%
	as.data.frame( . ) %>%
	unique( . )
	
sleuth_genemap <- transcript_annotation %>%
	select( .data = ., transcript_id, gene_id, symbol ) %>%
	rename( .data = ., target_id = transcript_id )

###################
# Run PBMC design #
###################

case_level <- paste0( tissue, "_", case_name )
control_level <- paste0( tissue, "_", control_name )

model_name <- paste0( "Status", case_level )

output_file <- here::here( 'results', paste0( 'sleuth_', tissue, '_', case_name, "_vs_", control_name, "_", "onlyDeseq2Genes.csv" ) )

design <- filter( .data = read_groups, Tissue == tissue ) %>%
	mutate( .data = ., Timepoint = case_when(
		is.na( Timepoint ) ~ '0m',
		TRUE ~ Timepoint ) ) %>%
	filter( .data = ., Timepoint == '0m' ) %>%
	select( .data = ., RGSM, Status ) %>%
	as.data.frame( . ) %>%
	unique( . ) %>%
	mutate( .data =., path = here::here( 'output', 'salmon', RGSM, 'abundance.h5' ) ) %>%
	mutate( .data = ., Status = factor( Status, levels = c( control_level, case_level ) ) ) %>%
	select( .data =., RGSM, Status, path )
names( design )[1] = 'sample'

#genes	
#so <- sleuth_prep( design, ~Status, target_mapping = sleuth_genemap, aggregation_column = 'gene_id', max_bootstrap = salmon_bootstraps, extra_bootstrap_summary = TRUE )

so <- sleuth_prep( design, ~Status, max_bootstrap = salmon_bootstraps, extra_bootstrap_summary = TRUE, filter_target_id = transcript_ids_to_test )

so <- sleuth_fit( so )

models( so )

so <- sleuth_wt( so, which_beta = model_name )

results <- sleuth_results( so, model_name, show_all = TRUE ) %>%
	filter( .data = ., !is.na( qval ) ) %>%
	rename( .data = ., transcript_id = target_id ) %>%
	mutate( .data = ., transcript_id = str_replace( string = transcript_id, pattern = "\\.[0-9]+$", replacement = "" ) ) %>%
	merge( x = transcript_annotation, y = ., by = 'transcript_id', all.y = TRUE ) %>%
	mutate( .data = ., FoldChange = exp( b ) ) %>%
	mutate( .data = ., FoldChange = case_when(
		FoldChange >= 1.0 ~ FoldChange,
		TRUE ~ (-1) / FoldChange)) %>%
	arrange( .data = ., qval, pval )

write_csv( x = results, path = output_file )

# session info

sessionInfo()

