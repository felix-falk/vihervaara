## RefSeq Promoters Script

# This R script returns the calculated promoters as a BED file, as calculated
# from a given refseq.txt file. 

# Run in terminal as: Rscript refseq_promoters_script.R refseq.txt refseqpromoters.bed

# In the refseq file: 
# Column 3 represents the chromosome
# Column 4 represents strandness
# Column 5 represents the first nucleotide - 1 postion for + strand genes
# Column 6 represents the first nucleotide of - strand genes
# Column 13 represents the gene name

args <- commandArgs(trailingOnly = TRUE) # Get command-line arguments

if (length(args) < 1) {
  stop("Please provide a filename as an argument.")
} # Check if the correct number of arguments is provided

RefSeqFilepath <- args[1] # refseq curated filepath ("~/Documents/Vihervaara/hg19/ncbiRefSeqCuratedhg19.txt")

RefSeqPromotersFilepath <- args[2] # promoters calculated from refseq filepath ("~/Documents/Vihervaara/hg19/refSeqPromoters.bed")

#RefSeqFilepath = "~/Documents/Vihervaara/hg19/ncbiRefSeqCuratedhg19.txt"
#RefSeqPromotersFilepath = "~/Documents/Vihervaara/hg19/refSeqPromoters.bed"

library(dplyr)
library(tidyverse)
library(stringr)

refSeq <- read.delim(RefSeqFilepath, 
                     header=FALSE) # Import the hg19 RefSeq select genome file

refSeq <- refSeq[, c(3, 4, 5, 6, 13)] # Pick out relevant columns

# refSeq <- refSeq %>% select(, c(1, 2, 3, 4, 5))

refSeqPlus <- refSeq %>% filter(V4 == "+") # Create "+" strand refseq data frame

refSeqMinus <- refSeq %>% filter(V4 == "-") # Create "-" strand refseq data frame

refSeqPlus <- refSeqPlus %>% 
  mutate(promoterStart = V5 - 100, 
         promoterEnd = V5) # Calculate "+" strand promoter coordinates

refSeqMinus <- refSeqMinus %>% 
  mutate(promoterStart = V6 + 1, 
         promoterEnd = V6 + 101) # Calculate "-" strand promoter coordinates

refSeqPromoters <- rbind(refSeqPlus, 
                         refSeqMinus) # Join back into one data frame

zeros <- rep(0, nrow(refSeqPromoters))

refSeqPromoters$zeros <- zeros # Create list of 0s for BED-file formatting

refSeqPromoters <- refSeqPromoters[
  , c("V3", 
      "promoterStart", 
      "promoterEnd", 
      "V13", 
      "zeros", 
      "V4")] # Create refSeqPromoters dataframe in BED-format

refSeqPromoters <- refSeqPromoters[
  !grepl("fix", refSeqPromoters$V3), ] # Remove RefSeq "fix" items 

refSeqPromoters <- refSeqPromoters[
  !grepl("alt", refSeqPromoters$V3), ] # Remove RefSeq "alt" items 

write.table(refSeqPromoters, 
            RefSeqPromotersFilepath, 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t") # Export refseq promoters BED-file
