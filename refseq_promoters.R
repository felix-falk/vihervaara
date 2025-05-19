## RefSeq Promoters Script

# This R script returns the calculated promoters as a BED file, as calculated
# from a given refseq.txt file. 

# Run in terminal as: Rscript refseq_promoters_script.R refseq.txt refseqpromoters.bed

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE) 

# Check if the correct number of arguments is provided
if (length(args) < 1) {
  stop("Please provide a filename as an argument.")
} 

# Import file paths
RefSeqFilepath <- args[1]
RefSeqPromotersFilepath <- args[2]

#RefSeqFilepath = "~/Documents/Vihervaara/hg19/ncbiRefSeqCuratedhg19.txt"
#RefSeqPromotersFilepath = "~/Documents/Vihervaara/hg19/refSeqPromoters.bed"

# Load required libraries
library(dplyr)
library(tidyverse)
library(stringr)

# Import the RefSeq select genome file
refSeq <- read.delim(RefSeqFilepath, header=FALSE) 

# Keep relevant columns
refSeq <- refSeq[, c(3, 4, 5, 6, 13)] 

# Create "+" strand refseq data frame
refSeqPlus <- refSeq %>% filter(V4 == "+") 

# Create "-" strand refseq data frame
refSeqMinus <- refSeq %>% filter(V4 == "-") 

# Calculate "+" strand promoter coordinates
refSeqPlus <- refSeqPlus %>% mutate(promoterStart = V5 - 100, promoterEnd = V5) 

# Calculate "-" strand promoter coordinates
refSeqMinus <- refSeqMinus %>% mutate(promoterStart = V6 + 1, promoterEnd = V6 + 101) 

# Join back into one data frame
refSeqPromoters <- rbind(refSeqPlus, refSeqMinus) 

# Create list of 0s for .bed file formatting
zeros <- rep(0, nrow(refSeqPromoters))
refSeqPromoters$zeros <- zeros 

# Create refSeqPromoters dataframe in BED-format
refSeqPromoters <- refSeqPromoters[, c("V3", "promoterStart", "promoterEnd", 
                                       "V13", "zeros", "V4")] 

# Remove RefSeq "fix" items 
refSeqPromoters <- refSeqPromoters[!grepl("fix", refSeqPromoters$V3), ] 

# Remove RefSeq "alt" items 
refSeqPromoters <- refSeqPromoters[!grepl("alt", refSeqPromoters$V3), ] 

# Export refseq promoters .bed file
write.table(refSeqPromoters, RefSeqPromotersFilepath, quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = "\t") 
