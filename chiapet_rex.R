# Find the REX motif density before, in promoter, between promoter and enhancer, 
# in enhancer and after, based on the filtered ChIA-PET file. 

# Load required packages
library(dplyr); library(ggplot2); library(tidyr); 
library(stringr); library(data.table); library(gridExtra); library(grid)

# Open REX BED file
rex <- read.delim("~/Documents/Vihervaara/hg19/hg19_chr/rex_new.bed", header=FALSE)
colnames(rex) <- c("chr", "rexStart", "rexEnd")
rex$rexMid <- rowMeans(rex[, c("rexStart", "rexEnd")], na.rm = TRUE)

# Import consolidated connections files
PRO_HS_connections <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_HS_connections.bed", header=FALSE)
CAGE_hemin_connections_0min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_0min.bed", header=FALSE)
CAGE_hemin_connections_15min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_15min.bed", header=FALSE)
CAGE_hemin_connections_30min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_30min.bed", header=FALSE)
CAGE_hemin_connections_60min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_60min.bed", header=FALSE)
CAGE_hemin_connections_24h <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_24h.bed", header=FALSE)
CAGE_hemin_connections_48h <- read.delim("~/Documents/Vihervaara/hg19/chiapet/CAGE_hemin_connections_48h.bed", header=FALSE)
eClusters_connections <- read.delim("~/Documents/Vihervaara/hg19/chiapet/eClusters_connections.bed", header=FALSE)
PRO_hemin_connections_0min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_0min.bed", header=FALSE)
PRO_hemin_connections_15min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_15min.bed", header=FALSE)
PRO_hemin_connections_30min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_30min.bed", header=FALSE)
PRO_hemin_connections_60min <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_60min.bed", header=FALSE)
PRO_hemin_connections_24h <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_24h.bed", header=FALSE)
PRO_hemin_connections_48h <- read.delim("~/Documents/Vihervaara/hg19/chiapet/PRO_hemin_connections_48h.bed", header=FALSE)

# Rename consolidated columns
colnames(PRO_HS_connections) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(eClusters_connections) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_0min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_15min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_30min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_60min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_24h) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(CAGE_hemin_connections_48h) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_0min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_15min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_30min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_60min) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_24h) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")
colnames(PRO_hemin_connections_48h) <- c("name", "chr", "eStart", "eEnd", "feature", "chrII", "chiaStart", "chiaEnd", "chrIII", "proStart", "proEnd", "gene", "score", "strand", "chrIV", "chiaStartII", "chiaEndII", "orientation", "distance")

# Add HEMIN column
PRO_HS_connections$hemin <- "heatshock"
eClusters_connections$hemin <- "eClusters"
CAGE_hemin_connections_0min$hemin <- "0min"
CAGE_hemin_connections_15min$hemin <- "15min"
CAGE_hemin_connections_30min$hemin <- "30min"
CAGE_hemin_connections_60min$hemin <- "60min"
CAGE_hemin_connections_24h$hemin <- "24h"
CAGE_hemin_connections_48h$hemin <- "48h"
PRO_hemin_connections_0min$hemin <- "0min"
PRO_hemin_connections_15min$hemin <- "15min"
PRO_hemin_connections_30min$hemin <- "30min"
PRO_hemin_connections_60min$hemin <- "60min"
PRO_hemin_connections_24h$hemin <- "24h"
PRO_hemin_connections_48h$hemin <- "48h"

PRO_HS_connections$seq <- "proseq"
eClusters_connections$seq <- "proseq"
CAGE_hemin_connections_0min$seq <- "cageseq"
CAGE_hemin_connections_15min$seq <- "cageseq"
CAGE_hemin_connections_30min$seq <- "cageseq"
CAGE_hemin_connections_60min$seq <- "cageseq"
CAGE_hemin_connections_24h$seq <- "cageseq"
CAGE_hemin_connections_48h$seq <- "cageseq"
PRO_hemin_connections_0min$seq <- "proseq"
PRO_hemin_connections_15min$seq <- "proseq"
PRO_hemin_connections_30min$seq <- "proseq"
PRO_hemin_connections_60min$seq <- "proseq"
PRO_hemin_connections_24h$seq <- "proseq"
PRO_hemin_connections_48h$seq <- "proseq"

# Remove duplicate eClusters
eClusters_connections <- eClusters_connections %>%
  group_by(feature) %>%
  arrange(desc(distance)) %>%  # Sort by distance in descending order
  slice_head(n = 1) %>%  # Keep the first row (the one with the highest distance)
  ungroup()  # Ungroup to return to normal dataframe

chiapet_connections <- rbind(PRO_HS_connections, eClusters_connections, 
                             CAGE_hemin_connections_0min, CAGE_hemin_connections_15min, 
                             CAGE_hemin_connections_30min, CAGE_hemin_connections_60min,
                             CAGE_hemin_connections_24h, CAGE_hemin_connections_48h, 
                             PRO_hemin_connections_0min, PRO_hemin_connections_15min, 
                             PRO_hemin_connections_30min,PRO_hemin_connections_60min, 
                             PRO_hemin_connections_24h, PRO_hemin_connections_48h)

# HEMIN Histogram Plot
chiapet_connections %>%
  filter(hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")) %>%
  mutate(hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h"))) %>%
  ggplot(aes(x = log10(distance))) +
  geom_histogram(binwidth = 0.1) +
  facet_grid(seq ~ hemin) +
  geom_text(data = ~ .x %>% count(seq, hemin, name = "n"),
            aes(x = Inf, y = Inf, label = paste("n =", n)),
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5) +
  labs(x = "log10(Promoter-Enhancer Distance)", y = "Count") +
  theme_minimal() +
  theme(legend.position = "none")

# eClusters and HS histogram plot
data <- chiapet_connections %>%
  filter(grepl("eClusters|heatshock", hemin)) %>%
  mutate(features_category = factor(case_when(
    grepl("Nhs", feature) ~ "Nhs",
    grepl("HS", feature) ~ "HS",
    grepl("both", feature) ~ "both",
    TRUE ~ "enhancer clusters"
  ), levels = c("Nhs", "HS", "both", "enhancer clusters")))

ggplot(data, aes(x = log10(distance))) +
  geom_histogram(binwidth = 0.1) +
  facet_wrap(~features_category, nrow = 1) +
  geom_text(data = count(data, features_category) %>% 
              mutate(label = paste0("n = ", n)), 
            aes(x = Inf, y = Inf, label = label), 
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.1) +
  labs(x = "log10(Promoter-Enhancer Distance)", y = "Count") +
  theme_minimal(base_size = 10) + theme(legend.position = "none")

# Convert to data.table
chiapet_connections <- as.data.table(chiapet_connections)
rex <- as.data.table(rex)

# Precompute REX as list by chromosome for fast access
rex_by_chr <- split(rex$rexMid, rex$chr)

# Function to count and compute density efficiently
count_density <- function(chr, start, end) {
  if (is.na(chr) || is.na(start) || is.na(end) || end <= start) {
    return(c(0L, 0))
  }
  
  rex_pos <- rex_by_chr[[chr]]
  if (is.null(rex_pos)) return(c(0L, 0))  # no REX on this chromosome
  
  len <- end - start
  count <- sum(rex_pos >= start & rex_pos <= end, na.rm = TRUE)
  return(c(count, count / len))
}

# Vectorized apply for enhancer, promoter, cluster
chiapet_connections[, c("enhancer_REX", "enhancer_density") := 
     transpose(mapply(count_density, chr, eStart, eEnd, SIMPLIFY = FALSE))]

chiapet_connections[, c("promoter_REX", "promoter_density") := 
     transpose(mapply(count_density, chr, proStart, proEnd, SIMPLIFY = FALSE))]

## STATISTICS

# Mean REX density in regions compared to the value expected by chance
# Expected number of REX motifs by chance = 0.00073242187
# Observed REX density genome-wide = 0.00101932746

mean((chiapet_connections %>% filter(seq == "proseq"))$enhancer_density) / 0.00073242187
mean((chiapet_connections %>% filter(seq == "cageseq"))$enhancer_density) / 0.00073242187
mean((chiapet_connections %>% filter(seq == "proseq"))$promoter_density) / 0.00073242187 
mean((chiapet_connections %>% filter(seq == "cageseq"))$promoter_density) / 0.00073242187 
mean((chiapet_connections %>% filter(hemin == "eClusters"))$enhancer_density) / 0.00073242187 

mean((chiapet_connections %>% filter(seq == "proseq"))$enhancer_density) / 0.00101932746
mean((chiapet_connections %>% filter(seq == "cageseq"))$enhancer_density) / 0.00101932746
mean((chiapet_connections %>% filter(seq == "proseq"))$promoter_density) / 0.00101932746 
mean((chiapet_connections %>% filter(seq == "cageseq"))$promoter_density) / 0.00101932746 
mean((chiapet_connections %>% filter(hemin == "eClusters"))$enhancer_density) / 0.00101932746 

# Correlation between REX density and DISTANCE

# Per HEMIN dataset
# Calculate AOV p-value and R-squared for each combination of seq and hemin
stats <- chiapet_connections %>%
  mutate(hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h"))) %>%
  filter(hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")) %>%
  group_by(seq, hemin) %>%
  do({
    model <- lm(log10(distance) ~ enhancer_density, data = .)
    aov_res <- summary(aov(model))
    p_value <- aov_res[[1]]["enhancer_density", "Pr(>F)"]
    r_squared <- summary(model)$r.squared
    data.frame(p_value = p_value, r_squared = r_squared)
  }) %>%
  ungroup() %>%
  mutate(
    hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h")),  # ensure factor also here
    label = paste0("R² = ", round(r_squared, 3), "\np = ", format.pval(p_value, digits = 3))
  )

# Now plot
ggplot(chiapet_connections %>%
         mutate(hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h"))) %>% 
         filter(hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")), 
       aes(x = enhancer_density, y = log10(distance))) +
  geom_point() +  # Scatter points
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +  # Linear regression line
  labs(x = "Enhancer REX-motif Density",
       y = "log10(Distance)") +
  facet_grid(seq ~ hemin) +  # Facet by seq (rows) and hemin (columns)
  theme_minimal() +
  geom_text(
    data = stats, 
    aes(x = Inf, y = Inf, label = label), 
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE
  )

# REX density in PROMOTERS for HEMIN DATA
# Calculate AOV p-value and R-squared for each combination of seq and hemin
stats <- chiapet_connections %>% 
  mutate(hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h"))) %>% 
  filter(hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")) %>%
  group_by(seq, hemin) %>%
  do({
    model <- lm(log10(distance) ~ promoter_density, data = .)
    aov_res <- summary(aov(model))
    p_value <- aov_res[[1]]["promoter_density", "Pr(>F)"]
    r_squared <- summary(model)$r.squared
    data.frame(p_value = p_value, r_squared = r_squared)
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("R² = ", round(r_squared, 3), "\np = ", format.pval(p_value, digits = 3))
  )

# Plot
ggplot(chiapet_connections %>%
         mutate(hemin = factor(hemin, levels = c("0min", "15min", "30min", "60min", "24h", "48h"))) %>% 
         filter(hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")), 
       aes(x = promoter_density, y = log10(distance))) +
  geom_point(alpha = 0.5) +  # Scatter points
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +  # Linear regression line
  labs(x = "Promoter Density",
       y = "log10(Distance)") +
  facet_grid(seq ~ hemin) +  # Facet by seq (rows) and hemin (columns)
  theme_minimal() +
  geom_text(
    data = stats, 
    aes(x = Inf, y = Inf, label = label), 
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE)

# REX density in ENHANCER CLUSTERS and ENHANCERS for HEAT SHOCK DATA
# Prepare the data
chiapet_connections_combined <- chiapet_connections %>%
  filter(hemin %in% c("eClusters", "heatshock")) %>%
  mutate(group = case_when(
    hemin == "eClusters" ~ "Enhancer Clusters",
    grepl("Nhs", feature) ~ "Nhs",
    grepl("HS", feature) ~ "HS",
    grepl("both", feature) ~ "both",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = c("Nhs", "HS", "both", "Enhancer Clusters")))  # specify order

# Calculate R² and p-value for each group
stats <- chiapet_connections_combined %>%
  group_by(group) %>%
  do({
    model <- lm(log10(distance) ~ enhancer_density, data = .)
    aov_res <- summary(aov(model))
    p_value <- aov_res[[1]]["enhancer_density", "Pr(>F)"]
    r_squared <- summary(model)$r.squared
    data.frame(p_value = p_value, r_squared = r_squared)
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("R² = ", round(r_squared, 3), "\np = ", format.pval(p_value, digits = 3))
  )

# Plot
ggplot(chiapet_connections_combined, aes(x = enhancer_density, y = log10(distance))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_text(
    data = stats, 
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.1) +
  labs(x = "Enhancer REX-motif Density", 
       y = "log10(Promoter-Enhancer Distance)") +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal()

# Correlation between REX DENSITY and eCluster and HEAT SHOCK data PROMOTER REX density
# Prepare the data
chiapet_connections_combined <- chiapet_connections %>%
  filter(hemin %in% c("eClusters", "heatshock")) %>%
  mutate(group = case_when(
    hemin == "eClusters" ~ "Enhancer Clusters",
    grepl("Nhs", feature) ~ "Nhs",
    grepl("HS", feature) ~ "HS",
    grepl("both", feature) ~ "both",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group)) %>%
  mutate(group = factor(group, levels = c("Nhs", "HS", "both", "Enhancer Clusters")))  # specify order

# Calculate R² and p-value for each group
stats <- chiapet_connections_combined %>%
  group_by(group) %>%
  do({
    model <- lm(log10(distance) ~ promoter_density, data = .)
    aov_res <- summary(aov(model))
    p_value <- aov_res[[1]]["promoter_density", "Pr(>F)"]
    r_squared <- summary(model)$r.squared
    data.frame(p_value = p_value, r_squared = r_squared)
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("R² = ", round(r_squared, 3), "\np = ", format.pval(p_value, digits = 3))
  )

# Plot
ggplot(chiapet_connections_combined, aes(x = promoter_density, y = log10(distance))) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_text(
    data = stats, 
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.1
  ) +
  labs(x = "Promoter REX-motif Density", 
       y = "log10(Promoter-Enhancer Distance)"
  ) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal()

## STATISTICAL TESTING

## REX motif ABSENCE or PRESENCE comparison

# Process, reshape, compute p-values, and plot in a streamlined pipeline
plot <- chiapet_connections %>%
  mutate(
    condition = case_when(
      seq == "cageseq" ~ "CAGE-seq",
      seq == "proseq" & !grepl("Nhs|HS|both", feature) & hemin != "eClusters" ~ "PRO-seq Hemin",
      feature %in% c("both", "Nhs", "HS") ~ "PRO-seq HS",
      hemin == "eClusters" ~ "eClusters",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(condition)) %>%
  pivot_longer(c(enhancer_REX, promoter_REX), names_to = "REX_location", values_to = "REX_value") %>%
  mutate(
    REX_presence = ifelse(REX_value == 0, "Absent", "Present"),
    REX_location = recode(REX_location, enhancer_REX = "Enhancer", promoter_REX = "Promoter"),
    log_distance = log10(distance)
  ) -> processed_data

# Compute p-values
p_values_df <- processed_data %>%
  group_by(condition, REX_location) %>%
  summarise(
    p_value = t.test(log_distance ~ REX_presence)$p.value,
    y_pos = max(log_distance, na.rm = TRUE) * 0.95,
    .groups = 'drop'
  )

# Generate plot
ggplot(processed_data, aes(x = REX_presence, y = log_distance)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  facet_grid(REX_location ~ condition) +
  geom_text(data = p_values_df,
            aes(x = 1.5, y = y_pos,
                label = paste0("p = ", format.pval(p_value, digits = 3, eps = .001))),
            inherit.aes = FALSE, size = 5) +
  labs(
    x = "REX Motif in Element",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_minimal(base_size = 14)

#################


### Write datasets to csv files to be used for GO enrichment

write_genes <- function(filter_expr, fname) {
  write.table(
    chiapet_connections %>% filter(!!rlang::enquo(filter_expr)) %>% pull(gene),
    file = fname, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ","
  )
}

# GO enrichment lists
write_genes(hemin == "eClusters", "eClusters_genes.csv")
for (f in c("Nhs", "HS", "both")) write_genes(feature == f, paste0("proseq_", f, "_genes.csv"))
for (s in c("proseq", "cageseq"))
  for (t in c("0min", "15min", "30min", "60min", "24h", "48h"))
    write_genes(hemin == t & seq == s, paste0(s, "_", t, "_genes.csv"))
write_genes(TRUE, "~/Desktop/GO_reference_genes.csv")

# REX presence/absence lists
rex_filters <- list(
  eClusters = quote(hemin == "eClusters"),
  heatshock = quote(feature %in% c("Nhs", "HS", "both")),
  hemin_proseq = quote(seq == "proseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h")),
  hemin_cageseq = quote(seq == "cageseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h"))
)

for (name in names(rex_filters)) {
  filter_base <- rex_filters[[name]]
  for (region in c("enhancer", "promoter")) {
    for (rex in c("absent", "present")) {
      rex_expr <- if (rex == "absent") {
        rlang::expr(!!rlang::sym(paste0(region, "_REX")) == 0)
      } else {
        rlang::expr(!!rlang::sym(paste0(region, "_REX")) != 0)
      }
      full_expr <- rlang::expr((!!filter_base) & (!!rex_expr))
      fname <- paste0("~/Desktop/", name, "_", region, "_REX_", rex, "_genes.csv")
      write_genes(!!full_expr, fname)
    }
  }
}

# Export gene lists for top 5000 RPK connections
write.table(CAGE_hemin_connections_0min %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_0min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(CAGE_hemin_connections_15min %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_15min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(CAGE_hemin_connections_30min %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_30min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(CAGE_hemin_connections_60min %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_60min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(CAGE_hemin_connections_24h %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_24h_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(CAGE_hemin_connections_48h %>% arrange(desc(feature)) %>% slice_head(n = 3000) %>% pull(gene), file = "CAGE_48h_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

write.table(PRO_hemin_connections_0min %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_0min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(PRO_hemin_connections_15min %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_15min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(PRO_hemin_connections_30min %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_30min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(PRO_hemin_connections_60min %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_60min_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(PRO_hemin_connections_24h %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_24h_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(PRO_hemin_connections_48h %>% arrange(desc(feature)) %>% slice_head(n = 10000) %>% pull(gene), file = "PRO_48h_gene_list.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

