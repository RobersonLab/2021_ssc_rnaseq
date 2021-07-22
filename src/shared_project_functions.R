###########################
# libraries required here #
###########################
library( here )
library( tidyverse )
library( ggrepel )

####################
# shared constants #
####################
min_reads_detection_threshold <- 5
min_detected_samples_for_de <- 5

##################
# directory vars #
##################
results_dir <- here::here( 'results' )
wgcna_output_dir <- here::here( 'results', 'wgcna' )
r_object_dir <- here::here( 'results', 'r_objects' )
figure_dir <- here::here( 'results', 'figures' )
data_dir <- here::here( 'data' )
salmon_dir <- here::here( 'output', 'salmon' )

##############################
# make directories as needed #
##############################
dir.create( path = results_dir, showWarnings = FALSE )
dir.create( path = figure_dir, showWarnings = FALSE )
dir.create( path = wgcna_output_dir, showWarnings = FALSE )
dir.create( path = r_object_dir, showWarnings = FALSE )
dir.create( path = file.path( figure_dir, 'wgcna' ), showWarnings = FALSE )

#########################################
# Some sprintf magic for summary tables #
#########################################
n_plus_percent_string <- function( query_value, total_value ) {
  percent_string <- ( ( query_value / total_value ) * 100.0 ) %>%
    round( x = ., digits = 1 ) %>%
    sprintf( "%.1f", . )

  paste0( query_value, " (", percent_string, ")" ) %>%
    return( . )
}

#########################################
# more sprintf magic for summary tables #
#########################################
mean_sd_string <- function( mean_value, sd_value ) {
  mean_string <- round( x = mean_value, digits = 1 ) %>%
    sprintf( "%.1f", . )
  sd_string <- round( x = sd_value, digits = 1 ) %>%
    sprintf( "%.1f", . )

  paste0( mean_string, " (", sd_string, ")" ) %>%
    return( . )
}

###########################
# Detect values in matrix #
###########################
detect <- function( input_vector, cutoff, upper = TRUE )
{
  if ( upper == TRUE )
  {
    return( length( which( input_vector >= cutoff ) ) )
  } else
  {
    return( length( which( input_vector <= cutoff ) ) )
  }
}

########################################
# non-parametric correlation for purrr #
########################################
kendall_tau_matrix_corr <- function( x_vector, y_vector ) {
  if ( !all( names( x_vector ) == names( y_vector ) ) ) {
    stop( paste0( "Variables out of order for cor.test" ) )
  }

  kendall_out <- cor.test( x = x_vector, y = y_vector, alternative = 'two.sided', method = 'kendall', exact = FALSE )

  c( Tau = unname( kendall_out$estimate ), P_kendall = unname( kendall_out$p.value ) ) %>%
    return( . )
}

###########################
# merge gtf featurecounts #
###########################
readFeatureCounts <- function( path_vector, rgsm_name ) {
  tb <- tibble()

  for ( idx in 1:length( path_vector ) ) {
    curr_rgid <- path_vector[ idx ]

    temp_tb <- read_tsv( file = curr_rgid, skip = 2, col_names = c( 'gene_id', 'chrom', 'start', 'end', 'strand', 'length', 'tempname' ), col_types = 'cccccii' ) %>%
      select( .data = ., gene_id, tempname )
    colnames( temp_tb )[2] <- paste0( 'tempname', idx )

    if ( sum( dim( tb ) ) > 0 ) {
      tb <- merge( x = tb, y = temp_tb, by = 'gene_id', all = TRUE )
    } else {
      tb = temp_tb
    }
  }

  count_mat <- tb %>%
    as.data.frame( . ) %>%
    column_to_rownames( ., 'gene_id' ) %>%
    as.matrix( . )

  count_tibble <- tibble( gene_id = rownames( count_mat ), temp = rowSums( count_mat ) )
  colnames( count_tibble )[2] = rgsm_name

  return( count_tibble )
}

################################
# pretty format for DE results #
################################
format_deseq_results <- function( deseq_results, annotation_df ) {
  as.data.frame( deseq_results ) %>%
    rownames_to_column( "gene_id" ) %>%
    filter( !is.na( padj ) ) %>%
    filter( !is.na( pvalue ) ) %>%
    mutate( FoldChange = 2^log2FoldChange ) %>%
    mutate( FoldChange = case_when(
      FoldChange < 1.0 ~ -1 / FoldChange,
      TRUE ~ FoldChange )
    ) %>%
    dplyr::rename( pval = pvalue ) %>%
    dplyr::rename( qval = padj ) %>%
    select( gene_id, FoldChange, pval, qval, baseMean, log2FoldChange, lfcSE, stat ) %>%
    merge( annotation_df, ., by="gene_id" ) %>%
    mutate_at( .tbl = ., .vars = c( 'FoldChange', 'pval', 'qval', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat' ), .funs = as.double ) %>%
    arrange( qval, pval, gene_id ) %>%
    return( . )
}

##################################
# Report differential expression #
##################################
summarize_pretty_de_results <- function( pretty_de, fc_cutoff ) {
  down <- pretty_de %>%
    filter( qval < 0.05 & FoldChange < 1.0 ) %>%
    pull( gene_id ) %>%
    length( . )

  down_fc <- pretty_de %>%
    filter( qval < 0.05 & FoldChange < ( -1 * fc_cutoff ) ) %>%
    pull( gene_id ) %>%
    length( . )

  up <- pretty_de %>%
    filter( qval < 0.05 & FoldChange > 1.0 ) %>%
    pull( gene_id ) %>%
    length( . )

  up_fc <- pretty_de %>%
    filter( qval < 0.05 & FoldChange > fc_cutoff ) %>%
    pull( gene_id ) %>%
    length( . )

  data.frame( Change = c( 'Down', paste0( 'Down - min ', sprintf( fmt="%.2f", fc_cutoff ) ), 'Up', paste0( 'Up - min ', sprintf( fmt="%.2f", fc_cutoff ) ) ), n = c( down, down_fc, up, up_fc ) ) %>%
    return( . )
}

#################
# volcano plots #
#################
make_ggplot_volcano <- function( deg_dataframe, case_name, control_name, axis_steps = 2, fold_change_cutoff = 1.5, qvalue_cutoff = 0.05, max_label = 30 )
{
  ##############################
  # set significance threshold #
  ##############################
  deg_dataframe <- deg_dataframe %>%
    mutate( Significant = case_when(
      qval < qvalue_cutoff & abs( FoldChange ) >= fold_change_cutoff ~ "Large",
      qval < qvalue_cutoff ~ "Modest",
      TRUE ~ "Not" ) ) %>%
    mutate( Significant = factor( Significant, levels=c( "Not", "Modest", "Large" ) ) )

  ################################
  # set values for square x axis #
  ################################
  x_volcano_value <- ( abs( deg_dataframe$log2FoldChange[ is.finite( deg_dataframe$log2FoldChange ) ] ) + 0.051 ) %>%
    max( . ) %>%
    round( ., 1 )

  if ( x_volcano_value < 1.0 ) {
    x_volcano_value = 1.0
  }

  x_num_for_limits <- round( x_volcano_value, 0 )

  x_volcano_low <- x_volcano_value * -1
  x_volcano_high <- x_volcano_value

  x_break_list <- seq( -1 * x_num_for_limits, x_num_for_limits, by = axis_steps )

  ##############
  # plot lines #
  ##############
  horizontal_line <- log10( qvalue_cutoff ) * -1
  vertical_line_1 <- log2( fold_change_cutoff )
  vertical_line_2 <- vertical_line_1 * -1

  ###################################
  # actually make the volcano plots #
  ###################################
  plot_volcano <- ggplot( deg_dataframe, aes( x=log2FoldChange, y=-log10( qval ), colour=Significant ) ) +
    scale_colour_manual( values = c( "darkgray", "blue", "red" ) ) +
    scale_x_continuous( limits = c( x_volcano_low, x_volcano_high ), breaks = x_break_list ) +
    theme_bw() +
    gg_bigger_texts +
    gg_no_legend +
    gg_no_grid +
    gg_center_title +
    geom_point( size=1.2 ) +
    geom_hline( yintercept = horizontal_line, linetype=2 ) +
    geom_vline( xintercept=c( vertical_line_1, vertical_line_2 ), linetype=2 ) +
    geom_text_repel( data=subset( deg_dataframe, Significant == "Large" )[c(1:max_label),], colour="black", aes( label=symbol ), size=3 ) +
    xlab( parse( text=paste0( "log[2]~(", case_name, "/", control_name, ")" ) ) ) +
    ylab( parse( text = paste0( "-log[10]~(Adj.~p-value)" ) ) )

  return( plot_volcano )
}

#################################################
# process raw gprofile into an updated gprofile #
#################################################
convertRefseqList2SymbolList <- function( input_string, name_frame, string_sep = ",", out_sep=";" ) {
  str_split( string = input_string, pattern = string_sep ) %>%
    unlist( "," ) %>%
    map_chr( .x = ., ~ name_frame[ .x, 'symbol' ] ) %>%
    sort( . ) %>%
    paste( ., collapse = out_sep ) %>%
    return( . )
}

parse_raw_gprofile <- function( gprof_file_path, symbol_df, title.truncate = 65, title.wrap = 35 ) {
  new_filename_string <- str_replace( gprof_file_path, pattern = ".csv", "_processed.tsv" )

  this_tibble <- read_csv( file = gprof_file_path, comment="#" ) %>%
    mutate( name = str_replace( term_name, "^ +", "" ) ) %>%
    mutate( symbols = map( intersections, convertRefseqList2SymbolList, name_frame = symbol_df ) ) %>%
    unnest( cols = c( symbols ) ) %>%
    mutate( title = paste0( term_id, ' - ', term_name ) ) %>%
    mutate( title = str_replace_all( string = title, pattern = "\\; +match class\\: +[0-9]+", replacement = "" ) ) %>%
    mutate( title = str_trunc( title, width = title.truncate ) ) %>%
    mutate( title = str_wrap( title, width = title.wrap ) ) %>%
    select( adjusted_p_value, source, term_id, term_name, symbols, intersections, title ) %>%
    mutate( intersections = str_replace( string = intersections, pattern = ",", replacement = ";" ) ) %>%
    arrange( adjusted_p_value )

  select( this_tibble, -title ) %>%
    write_tsv( path = new_filename_string )

  return( this_tibble )
}

############################
# expression to PCA matrix #
############################
run_prcomp <- function( expression_matrix, sample_list, gene_list ) {
  input_matrix <- expression_matrix[ gene_list, sample_list ] %>%
    t( . ) %>%
    prcomp( x = ., scale. = TRUE ) %>%
    return( . )
}

###############################
# color blind friendly colors #
###############################
colorBlindPalette <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )
colorBlindPalette2 <- c( "#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )
blues5Palette <- c( '#ffffcc', '#a1dab4', '#41b6c4', '#2c7fb8', '#253494' )
greens5Palette <- c( '#ffffcc', '#c2e699', '#78c679', '#31a354', '#006837' )
purples5Palette <- c( '#feebe2', '#fbb4b9', '#f768a1', '#c51b8a', '#7a0177' )
reds5Palette <- c( '#ffffb2', '#f3cc5c', '#fd8d3c', '#f03b20', '#bd0026' )
reds4Palette <- c( '#fef0d9', '#fdcc8a', '#fc8d59', '#d7301f' )
reds3Palette <- c( '#fee0d2', '#fc9272', '#de2d26' )

####################
# ggplot modifiers #
####################
gg_bigger_texts = theme(
  axis.title = element_text( size = 22 ),
  axis.text = element_text( size = 20 ),
  legend.text = element_text( size = 14 ),
  legend.title = element_text( size = 15 ),
  plot.title = element_text( size = 22 ),
  strip.text.x = element_text( size = 17, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_multiplot_texts = theme(
  axis.title = element_text( size = 20 ),
  axis.text = element_text( size = 18 ),
  legend.text = element_text( size = 12 ),
  legend.title = element_text( size = 13 ),
  plot.title = element_text( size = 20 ),
  strip.text.x = element_text( size = 16, margin = margin( b = 5, t = 5 ) ),
  strip.text.y = element_text( size = 15 )
)

gg_quadplot_smaller_text = theme(
  axis.title = element_text( size=14 ),
  axis.text = element_text( size=9 ),
  plot.title = element_text( size=15 )
)

gg_reduce_pathway_text = theme(
    axis.title = element_text( size=14 ),
    axis.text.y = element_text( size=8 ),
    axis.text.x = element_text( size=10 ),
    plot.title = element_text( size=15 )
)

gg_no_legend = theme(
  legend.position='none'
)

gg_no_grid = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

gg_no_x_grid = theme(
  panel.grid.major.x = element_blank() )

gg_no_y_grid = theme(
  panel.grid.major.y = element_blank() )

gg_center_title = theme(
  plot.title = element_text( hjust = 0.5 )
)

gg_no_x_label = theme(
  axis.title.x = element_blank()
)

gg_no_y_label = theme(
  axis.title.y = element_blank()
)

gg_angled_x_text = theme (
  axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1, color = 'black' )
)

#################################
# Convert UpSet plot to cowplot #
#################################
UpsetPlot2Cowplot <- function( upset_object, title = NULL ) {
  # I tried using this to make a reasonable figure but it was pretty cruddy
  # but figured that it was worth saving in case I needed to work on it later

  library( UpSetR )
  library( ggplotify )
  library( grid )

  top_row_ratio <- c( 0.25, 0.75 )
  bottom_row_ratio <- c( 0.185, 0.815 )
  vertical_ratios <- c( 0.7, 0.3 )

  blank_plot <- grid.rect( gp = gpar( col = "white" ) )

  main_upset_plot <- as_grob( upset_object$Main_bar )

  top_row <- plot_grid( blank_plot,
                        main_upset_plot,
                        rel_widths = top_row_ratio )

  if ( !is.null( title ) ) {
    title_gg <- ggplot() + labs( title = title ) + gg_bigger_texts + gg_center_title

    top_row <- plot_grid( top_row, title_gg, ncol = 1, rel_heights = c( 0.10, 1 ) )
  }

  bottom_row <- plot_grid( as_grob( upset_object$Sizes ),
                           as_grob( upset_object$Matrix ),
                           rel_widths = bottom_row_ratio )

  plot_grid( top_row,
                           bottom_row,
                           rel_heights = vertical_ratios,
                           ncol = 1 ) %>%
    return( . )
}
