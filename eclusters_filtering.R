## Enhancer clusters filtering script

# Load required library
library(dplyr)

# Assign input and output filepaths
input <- "~/Documents/Vihervaara/hg19/enhancers/K562_hg19_eClusters_and_clusteredEnhancers_noHeader.bed"
output <- "~/Documents/Vihervaara/hg19/enhancers/eClusters_filtered.bed"

# Import enhancer clusters file
eClusters <- read.delim(input, header=FALSE)

# Keep relevant columns
eClusters <- eClusters %>% select(c("V1", "V2", "V3", "V4"))

# Keep treatment or cluster name in column 4
eClusters[[4]] <- sub(".*_(.*)$", "\\1", eClusters[[4]])

# Keep clusters
eClusters <- eClusters[!grepl("Nhs|HS|both", eClusters[[4]]), ]

# Create .bed file of filtered enhancer cluster dataframe
write.table(eClusters, output, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
