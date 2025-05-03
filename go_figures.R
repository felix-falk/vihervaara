# GO-enrichment results

library(ggplot2)
library(dplyr)

GO_results <- read.csv("~/Desktop/GO-enrichment-results.csv", sep=";")

colnames(GO_results) <- c("name", "GO_biological_process_complete", 
                                     "A", "B", "expected", "over_under", 
                                     "foldenrichment", "raw_p_value", "FDR", 
                                     "hemin", "sequencing_method")
GO_results <- GO_results %>% filter(!grepl("heatshock", hemin))

# Plot
ggplot(GO_results %>% filter(grepl("proseq", sequencing_method)), aes(x = foldenrichment, 
                       y = reorder(GO_biological_process_complete, foldenrichment), 
                       color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Hemin-treatmet, PRO-seq Enhancers, GO Term Enrichment",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(~ hemin, nrow = 1) +
  ylab(NULL)

ggplot(GO_results %>% filter(grepl("cageseq", sequencing_method)), aes(x = foldenrichment, 
                                                                      y = reorder(GO_biological_process_complete, foldenrichment), 
                                                                      color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Hemin-treatmet, CAGE-seq Enhancers, GO Term Enrichment",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(~ hemin, nrow = 1) +
  ylab(NULL)

GO_results_absence_presence <- read.csv("~/Desktop/GO_absence_presence_results.csv", sep=";")

ggplot(GO_results_absence_presence %>% filter(grepl("proseq", seq)), aes(x = foldenrichment, 
                                                                      y = reorder(GObiologicalprocesscomplete, foldenrichment), 
                                                                      color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "PRO-seq, GO Terms Enrichment",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(motif ~ REX, nrow = 1) +
  ylab(NULL)

ggplot(GO_results_absence_presence %>% filter(grepl("cageseq", seq)), aes(x = foldenrichment, 
                                                                         y = reorder(GObiologicalprocesscomplete, foldenrichment), 
                                                                         color = -log10(FDR))) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "CAGE-seq, GO Terms Enrichment",
       x = "Fold Enrichment",
       color = "-log10(FDR)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_wrap(motif ~ REX, nrow = 1) +
  ylab(NULL)


