## PRO-seq filtering script

library(dplyr)

# Import hemin enhancers file
PRO_hemin_enhancers <- read.delim("~/Documents/Vihervaara/hg19/proseq/K562_hg19_hemin_enhancers.txt")

# Keep columns of interest, those with eRPK data
PRO_hemin_enhancers <- PRO_hemin_enhancers %>% select(-c("He15_FCe", "He15_log2FCe", "He15_dRPKe", "He15_Reg_pl_mn", "He30_FCe", "He30_log2FCe", "He30_Reg_pl_mn", "He60_FCe", "He60_log2FCe", "He60_dRPKe", "He60_Reg_pl_mn", "He24h_FCe", "He24h_log2FCe", "He24h_dRPKe", "He24h_Reg_pl_mn", "He60R48h_FCe", "He60R48h_log2FCe", "He60R48h_dRPKe", "He60R48h_Reg_pl_mn", "He24hR48h_FCe", "He24hR48h_log2FCe", "He24hR48h_dRPKe", "He24hR48h_Reg_pl_mn", "He15_eReg", "He30_eReg", "He60_eReg", "He24h_eReg", "He60R48h_eReg", "He24hR48h_eReg"))

# Create a mean of the two 48h eRPK columns
PRO_hemin_enhancers <- PRO_hemin_enhancers %>%
  mutate(mean_He48h_eRPK = rowMeans(select(., "He60R48h_eRPK", "He24hR48h_eRPK"), na.rm = TRUE)) %>%
  select(-"He60R48h_eRPK", -"He24hR48h_eRPK")

# Set each eRPK below 5 to 0, as a filter, remove rows for which all are 0
PRO_hemin_enhancers <- PRO_hemin_enhancers %>%
  mutate(across(c(He0_eRPK, He15_eRPK, He30_eRPK, He60_eRPK, He24h_eRPK, mean_He48h_eRPK), ~ ifelse(. < 5, 0, .)))
PRO_hemin_enhancers <- PRO_hemin_enhancers %>%
  filter(!(He0_eRPK == 0 & He15_eRPK == 0 & He30_eRPK == 0 & He24h_eRPK == 0 & mean_He48h_eRPK == 0))

He0_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, He0_eRPK)) %>% filter(!(He0_eRPK == 0))

He15_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, He15_eRPK)) %>% filter(!(He15_eRPK == 0))

He30_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, He30_eRPK)) %>% filter(!(He30_eRPK == 0))

He60_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, He60_eRPK)) %>% filter(!(He60_eRPK == 0))

He24h_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, He24h_eRPK)) %>% filter(!(He24h_eRPK == 0))

He48h_eRPK <- PRO_hemin_enhancers %>% select(c(chr, eStart, eEnd, mean_He48h_eRPK)) %>% filter(!(mean_He48h_eRPK == 0))

write.table(He0_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He0_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(He15_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He15_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(He30_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He30_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(He60_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He60_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(He24h_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He24h_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

write.table(He48h_eRPK, "~/Documents/Vihervaara/hg19/enhancers/He48h_eRPK.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
