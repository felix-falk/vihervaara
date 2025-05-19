## PRO-seq heat shock filtering script

# Load required library
library(dplyr)

# Assign filepaths
input <- "~/Documents/Vihervaara/hg19/enhancers/hg19_K562_enhancers_NatCom_dTREs.bed"
output <- "~/Documents/Vihervaara/hg19/enhancers/PRO_HS_enhancers.bed"

# Import heat shock enhancers file
PRO_HS_enhancers <- read.delim(input, header=FALSE)

# Keep columns of interest
PRO_HS_enhancers <- PRO_HS_enhancers %>% select(c("V1", "V2", "V3", "V7"))

# Export enhancers data frame as .bed file
write.table(PRO_HS_enhancers, output, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
