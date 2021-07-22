#############
# libraries #
#############
library( tidyverse )
library( here )
library( sleuth )

#############
# variables #
#############
env_variables <- read_tsv( file = here::here( 'environment_variables.sh' ), col_names = FALSE, col_types = 'c', comment = '#' ) %>%
	separate( data = ., col = 'X1', sep = '=', into = c( 'name', 'value' ) )

salmon_bootstraps <- filter( .data = env_variables, name == 'SALMON_BOOTSTRAPS' ) %>%
	pull( .data = ., value ) %>%
	as.integer( . )
	
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

########
# Skin #
########
tissue <- "Skin"
case_name <- "SSc"
control_name <- "Control"

case_level <- paste0( tissue, "_", case_name )
control_level <- paste0( tissue, "_", control_name )

model_name <- paste0( "Status", case_level )

output_file <- here::here( 'results', paste0( 'sleuth_', tissue, '_', case_name, "_vs_", control_name, "_", "transcript_id.csv" ) )

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

so <- sleuth_prep( design, ~Status, max_bootstrap = salmon_bootstraps, extra_bootstrap_summary = TRUE )

so <- sleuth_fit( so )

models( so )

so <- sleuth_wt( so, which_beta = model_name )

results <- sleuth_results( so, model_name, show_all = TRUE ) %>%
	filter( .data = ., !is.na( qval ) ) %>%
	rename( .data = ., transcript_id = target_id ) %>%
	mutate( .data = ., transcript_id = str_replace( string = transcript_id, pattern = "\\.[0-9]+$", replacement = "" ) ) %>%
	merge( x = transcript_annotation, y = ., by = 'transcript_id', all.y = TRUE )  %>%
	mutate( .data = ., FoldChange = exp( b ) ) %>%
	mutate( .data = ., FoldChange = case_when(
		FoldChange >= 1.0 ~ FoldChange,
		TRUE ~ (-1) / FoldChange)) %>%
	arrange( .data = ., qval, pval )

write_csv( x = results, path = output_file )

# session info

sessionInfo()

