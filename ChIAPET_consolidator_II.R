#!/usr/bin/env Rscript

# CHIAPET Consolidator II
# Consolidates two chiapet connection files, the pro-en and en-pro files into one

# Run in BASH like so: 
# Rscript ~/Documents/Vihervaara/vihervaara_scripts/ChIAPET_consolidator_II.R ~/Documents/Vihervaara/hg19/chiapet/lEn_rPro_ChIAPET.bed ~/Documents/Vihervaara/hg19/chiapet/lPro_rEn_ChIAPET.bed ~/Documents/Vihervaara/hg19/chiapet/PRO_HS_connections.bed

# Load required library
library(dplyr)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 3) {
  stop("Please provide 3 arguments: input file 1, input file 2, and output file path")
}

# File paths from arguments
lEnhancer_filepath <- args[1]
rEnhancer_filepath <- args[2]
output_path <- args[3]


# Read the files
lEnhancer <- read.delim(lEnhancer_filepath, header = FALSE)
rEnhancer <- read.delim(rEnhancer_filepath, header = FALSE)

# Process the lEnhancer file
lEnhancer <- lEnhancer %>%
  mutate(orientation = "enpro") %>%
  mutate(distance = V10 - V4) %>%
  group_by(V1) %>%
  filter(distance == max(distance)) %>%
  slice(1) %>%
  ungroup()

# Process the rEnhancer file
rEnhancer <- rEnhancer %>%
  mutate(orientation = "proen") %>%
  mutate(distance = V3 - V11) %>%
  group_by(V1) %>%
  filter(distance == max(distance)) %>%
  slice(1) %>%
  ungroup()

# Combine the two datasets
merged <- rbind(lEnhancer, rEnhancer)

# Write the merged result to the output file
write.table(merged, output_path, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

cat("Merge complete! Output written to: ", output_path, "\n")
