library( wasabi )
library( here )

cargs <- commandArgs()
base_name <- cargs[3]

salmon_dir <- here::here( 'output', 'salmon', base_name )

prepare_fish_for_sleuth( salmon_dir )

sessionInfo()

