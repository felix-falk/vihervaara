## ChIA-PET Filtering script

# This R script returns the left-hand and right-hand ChIA-PET connection regions
# from a BED file as two BED files. 

# Run in terminal as: 
# Rscript ~/Documents/Vihervaara/vihervaara_scripts/chia_pet_script.R ~/Documents/Vihervaara/hg19/chiapet/ChIAPET_dTREd_to_Pol2blocks_allConn_iu_indivReducedToUniqueConnections_TCDinRed_d2dBlac_d2nTinBrOrange_nTnTgreenish.sorted.bed ~/Documents/Vihervaara/hg19/chiapet/ChIAPETLeft.bed ~/Documents/Vihervaara/hg19/chiapet/ChIAPETRight.bed

args <- commandArgs(trailingOnly = TRUE) # Get command-line arguments

if (length(args) < 3) {
  stop("Please provide a filename as an argument.")
} # Check if the correct number of arguments is provided

chiapetfilepath <- args[1] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPET_dTREd_to_Pol2blocks_allConn_iu_indivReducedToUniqueConnections_TCDinRed_d2dBlac_d2nTinBrOrange_nTnTgreenish.sorted.bed"

lefthandfilepath <- args[2] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPETLeft.bed"

righthandfilepath <- args[3] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPETRight.bed"

library(dplyr)
library(tidyverse)
library(stringr) # Load required packages

ChIAPET <- read.delim(chiapetfilepath, header=FALSE) # Import the ChIA-PET file

# Create a list of unique names, as long as the ChIA-PET file
uniqueNames <- paste("Name", seq(1, nrow(ChIAPET)), sep = "_")

ChIAPET$uniqueNames <- uniqueNames # Add the unique names as a new column in the data frame

# Column 2 is the start of the leftmost connection
# Column 3 is the end of the rightmost connection
# Column 4 contains the start and end coordinates of both connections

# Segregate out column 4

# Split up column 4 by the "-" character

ChIAPET <- ChIAPET %>%
  mutate(
    V4beforeDash = str_split(V4, "-", simplify = TRUE)[, 1],
    V4afterDash = str_split(V4, "-", simplify = TRUE)[, 2]
    )

ChIAPET <- ChIAPET %>% filter(sub(":.*", "", V4beforeDash) == sub(":.*", "", V4afterDash)) # Filter out all interchromosomal connections

ChIAPET$V4beforeDash <- sub("^[^:]+:", "", ChIAPET$V4beforeDash) # Remove beginning of V4beforeDash string

ChIAPET$V4afterDash <- sub("^[^:]+:", "", ChIAPET$V4afterDash) # Remove beginning of V4afterDash string

ChIAPET$V4afterDash <- sub(",.*$", "", ChIAPET$V4afterDash) # Remove the end of V4afterDash string

# Split up the V4beforeDasha and V4afterDash columns by the ".." character
ChIAPET <- ChIAPET %>%
  mutate(
    leftStart = str_split(V4beforeDash, "\\.\\.", simplify = TRUE)[, 1],
    leftEnd = str_split(V4beforeDash, "\\.\\.", simplify = TRUE)[, 2],
    rightStart = str_split(V4afterDash, "\\.\\.", simplify = TRUE)[, 1],
    rightEnd = str_split(V4afterDash, "\\.\\.", simplify = TRUE)[, 2]
  )

ChIAPETLeft <- ChIAPET[, c("V1", "leftStart", "leftEnd", "uniqueNames")]

ChIAPETRight <- ChIAPET[, c("V1", "rightStart", "rightEnd", "uniqueNames")]

write.table(ChIAPETLeft, lefthandfilepath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(ChIAPETRight, righthandfilepath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
