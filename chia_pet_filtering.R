## ChIA-PET Filtering script

# This R script returns the left-hand and right-hand ChIA-PET anchor regions
# from a .bed file as two .bed files. 

# Run in terminal as: 
# Rscript ~/Documents/Vihervaara/vihervaara_scripts/chia_pet_script.R ~/Documents/Vihervaara/hg19/chiapet/ChIAPET_dTREd_to_Pol2blocks_allConn_iu_indivReducedToUniqueConnections_TCDinRed_d2dBlac_d2nTinBrOrange_nTnTgreenish.sorted.bed ~/Documents/Vihervaara/hg19/chiapet/ChIAPETLeft.bed ~/Documents/Vihervaara/hg19/chiapet/ChIAPETRight.bed

args <- commandArgs(trailingOnly = TRUE) # Get command-line arguments

if (length(args) < 3) {
  stop("Please provide a filename as an argument.")
} # Check if the correct number of arguments is provided

chiapetfilepath <- args[1] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPET_dTREd_to_Pol2blocks_allConn_iu_indivReducedToUniqueConnections_TCDinRed_d2dBlac_d2nTinBrOrange_nTnTgreenish.sorted.bed"

lefthandfilepath <- args[2] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPETLeft.bed"

righthandfilepath <- args[3] # "~/Documents/Vihervaara/hg19/chiapet/ChIAPETRight.bed"

# Load required packages
library(dplyr)
library(tidyverse)
library(stringr) 

# Import the ChIA-PET file
ChIAPET <- read.delim(chiapetfilepath, header=FALSE) 

# Create a list of unique names, as long as the ChIA-PET file
uniqueNames <- paste("Name", seq(1, nrow(ChIAPET)), sep = "_")

# Add the unique names as a new column in the data frame
ChIAPET$uniqueNames <- uniqueNames 

# Column 2 is the start of the leftmost anchor
# Column 3 is the end of the rightmost anchor
# Column 4 contains the start and end coordinates of both anchors

# Segregate out column 4

# Split up column 4 by the "-" character

ChIAPET <- ChIAPET %>%
  mutate(
    V4beforeDash = str_split(V4, "-", simplify = TRUE)[, 1],
    V4afterDash = str_split(V4, "-", simplify = TRUE)[, 2]
    )

# Filter out all interchromosomal connections
ChIAPET <- ChIAPET %>% filter(sub(":.*", "", V4beforeDash) == sub(":.*", "", V4afterDash)) 

# Remove beginning of V4beforeDash string
ChIAPET$V4beforeDash <- sub("^[^:]+:", "", ChIAPET$V4beforeDash) 

# Remove beginning of V4afterDash string
ChIAPET$V4afterDash <- sub("^[^:]+:", "", ChIAPET$V4afterDash) 

# Remove the end of V4afterDash string
ChIAPET$V4afterDash <- sub(",.*$", "", ChIAPET$V4afterDash)

# Split up the V4beforeDasha and V4afterDash columns by the ".." character
ChIAPET <- ChIAPET %>%
  mutate(
    leftStart = str_split(V4beforeDash, "\\.\\.", simplify = TRUE)[, 1],
    leftEnd = str_split(V4beforeDash, "\\.\\.", simplify = TRUE)[, 2],
    rightStart = str_split(V4afterDash, "\\.\\.", simplify = TRUE)[, 1],
    rightEnd = str_split(V4afterDash, "\\.\\.", simplify = TRUE)[, 2]
  )

# Keep relevant rows
ChIAPETLeft <- ChIAPET[, c("V1", "leftStart", "leftEnd", "uniqueNames")]
ChIAPETRight <- ChIAPET[, c("V1", "rightStart", "rightEnd", "uniqueNames")]

# Export left hand ChIA-PET anchors .bed file
write.table(ChIAPETLeft, lefthandfilepath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

# Export right hand ChIA-PET anchors .bed file
write.table(ChIAPETRight, righthandfilepath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
