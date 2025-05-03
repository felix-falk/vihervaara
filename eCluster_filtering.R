## eCluster filtering script

library(dplyr)

# Import hemin enhancers file
eClusters <- K562_hg19_eClusters_and_clusteredEnhancers_noHeader <- read.delim("~/Documents/Vihervaara/hg19/enhancers/K562_hg19_eClusters_and_clusteredEnhancers_noHeader.bed", header=FALSE)

eClusters <- eClusters %>% select(c("V1", "V2", "V3", "V4"))

eClusters[[4]] <- sub(".*_(.*)$", "\\1", eClusters[[4]])

eClusters <- eClusters[!grepl("Nhs|HS|both", eClusters[[4]]), ]

write.table(eClusters, "~/Documents/Vihervaara/hg19/enhancers/eClusters_filtered.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
