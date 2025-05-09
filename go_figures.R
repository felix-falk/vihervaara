# GO-enrichment results

library(ggplot2)
library(dplyr)

# Load GO data
GO_results <- read.csv("~/Documents/GitHub/vihervaara/thesis_datasets/go_per_treatment_results.csv", sep=";")

colnames(GO_results) <- c("name", "GO_biological_process_complete", 
                                     "A", "B", "expected", "over_under", 
                                     "foldenrichment", "raw_p_value", "FDR", 
                                     "hemin", "sequencing_method")
GO_results <- GO_results %>% filter(!grepl("heatshock", hemin))

# Plot PRO-seq Hemin 
tiff("~/Documents/vihervaara/goterms_pro_figure.tiff", units="in", width=13, height=5, res=600)
GO_results %>%
  filter(grepl("proseq", sequencing_method)) %>%
  mutate(hemin = factor(hemin, levels = c("0 minutes", "15 minutes", "30 minutes", "60 minutes", "24 hours", "48 hours"))) %>%
  ggplot(aes(x = foldenrichment, 
             y = reorder(GO_biological_process_complete, foldenrichment), 
             color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO terms enrichment of PRO-seq enhancer-associated genes, per hemin treatment duration",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(~ hemin, nrow = 1) +
  ylab(NULL)
dev.off()

# Plot CAGE-seq Hemin
tiff("~/Documents/vihervaara/goterms_cage_figure.tiff", units="in", width=13, height=8, res=600)
GO_results %>%
  filter(grepl("cageseq", sequencing_method)) %>%
  mutate(hemin = factor(hemin, levels = c("0 minutes", "15 minutes", "30 minutes", "60 minutes", "24 hours", "48 hours"))) %>%
  ggplot(aes(x = foldenrichment, 
             y = reorder(GO_biological_process_complete, foldenrichment), 
             color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO terms enrichment of CAGE-seq enhancer-associated genes, per hemin treatment duration",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(~ hemin, nrow = 1) +
  ylab(NULL)
dev.off()

# Load REX absence / presence GO data
GO_results_absence_presence <- read.csv("~/Documents/GitHub/vihervaara/thesis_datasets/go_rex_presence_absence_results.csv", sep=";")

tiff("~/Documents/vihervaara/goterms_pro_absence_figure.tiff", units="in", width=10, height=5, res=600)
ggplot(GO_results_absence_presence %>% filter(grepl("proseq", seq)), aes(x = foldenrichment, 
                                                                      y = reorder(GObiologicalprocesscomplete, foldenrichment), 
                                                                      color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO terms enrichment, PRO-seq, REX presence or absence",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(motif ~ REX, nrow = 1) +
  ylab(NULL)
dev.off()

tiff("~/Documents/vihervaara/goterms_cage_absence_figure.tiff", units="in", width=10, height=5, res=600)
ggplot(GO_results_absence_presence %>% filter(grepl("cageseq", seq)), aes(x = foldenrichment, 
                                                                         y = reorder(GObiologicalprocesscomplete, foldenrichment), 
                                                                         color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO terms enrichment, CAGE-seq, REX presence or absence",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(motif ~ REX, nrow = 1) +
  ylab(NULL)
dev.off()

