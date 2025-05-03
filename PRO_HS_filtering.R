## PRO-seq filtering script

library(dplyr)

# Import hemin enhancers file
PRO_HS_enhancers <- read.delim("~/Documents/Vihervaara/hg19/enhancers/hg19_K562_enhancers_NatCom_dTREs.bed", header=FALSE)

PRO_HS_enhancers <- PRO_HS_enhancers %>% select(c("V1", "V2", "V3", "V7"))

write.table(PRO_HS_enhancers, "~/Documents/Vihervaara/hg19/enhancers/PRO_HS_enhancers.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
