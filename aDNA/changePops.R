#####################################################################
####
#### R Script to make each aDNA sequence belong to its own population
####
#####################################################################

cat("\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%% Changing Populations %%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")

## Set-up
library(genio, quietly = TRUE)

## Get individuals' data
individuals <- read_ind("Data/v53.1_1240k_pub_original.ind")
individuals$label <- individuals$id

## Write new ind file
write_ind("Data/v53.1_1240k_pub.ind", individuals)