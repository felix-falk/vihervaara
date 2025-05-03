# Find the REX motif density before, in promoter, between promoter and enhancer, 
# in enhancer and after, based on the filtered ChIA-PET file. 

# Load required packages
library(dplyr); library(ggplot2); library(tidyr); 
library(stringr); library(data.table); library(gridExtra); library(grid)

# Open REX BED file
rex <- read.delim("~/Documents/Vihervaara/hg19/hg19_rex_bed/rex_new.bed", header=FALSE)
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
ggplot(chiapet_filtered, aes(x = log10(distance))) + 
  geom_histogram(binwidth = 0.1, position = "dodge") +
  labs(title = "Distance by Sequencing Method and Hemin Treatment",
       x = "log10(Promoter-Enhancer Distance)", 
       y = "Count") +
  theme_minimal() +
  facet_grid(seq ~ hemin) +
  theme(legend.position = "none") +
  geom_text(data = chiapet_filtered %>% group_by(seq, hemin) %>% summarise(n = n(), .groups = "drop"), aes(x = Inf, y = Inf, label = paste("n =", n)), 
            inherit.aes = FALSE, hjust = 1.1, vjust = 1.5)

# eClusters and HS histogram plot
# Calculate the total count for each 'features_category'
count_labels <- chiapet_connections %>%
  filter(grepl("eClusters|heatshock", hemin)) %>%
  mutate(features_category = factor(case_when(
    grepl("Nhs", feature) ~ "Nhs",
    grepl("HS", feature) ~ "HS",
    grepl("both", feature) ~ "both",
    TRUE ~ "enhancer clusters"  # Change "Other" to "enhancer clusters"
  ), levels = c("Nhs", "HS", "both", "enhancer clusters"))) %>%
  group_by(features_category) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(label = paste0("n = ", n))  # Prepare the label with the count

# Plot the graph with the 'n' label for each facet
ggplot(chiapet_connections %>% 
         filter(grepl("eClusters|heatshock", hemin)) %>%
         mutate(features_category = factor(case_when(
           grepl("Nhs", feature) ~ "Nhs",
           grepl("HS", feature) ~ "HS",
           grepl("both", feature) ~ "both",
           TRUE ~ "enhancer clusters"  # Change "Other" to "enhancer clusters"
         ), levels = c("Nhs", "HS", "both", "enhancer clusters"))), 
       aes(x = log10(distance))) + 
  geom_histogram(binwidth = 0.1, position = "dodge") +  # Adjusted binwidth for more bins
  labs(
    title = "Distance by Sequencing Method and Data Set",
    x = "log10(Promoter-Enhancer Distance)", 
    y = "Count"
  ) +
  theme_minimal() +
  facet_wrap(~features_category, nrow = 1) +  # Split by the new 'features_category' variable
  theme(legend.position = "none") +  # Remove the legend
  geom_text(
    data = count_labels, 
    aes(x = Inf, y = Inf, label = label), 
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE  # Position the label in the upper-right corner
  )

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

# Vectorized apply for enhancer, promoter, cluster, within
chiapet_connections[, c("enhancer_REX", "enhancer_density") := 
     transpose(mapply(count_density, chr, eStart, eEnd, SIMPLIFY = FALSE))]

chiapet_connections[, c("promoter_REX", "promoter_density") := 
     transpose(mapply(count_density, chr, proStart, proEnd, SIMPLIFY = FALSE))]

## STATISTICS

# Mean REX density in regions compared to the value expected by chance
# Expected number of REX motifs by chance = 0.0009765625

mean(chiapet_connections$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance

mean((chiapet_connections %>% filter(hemin == "0min", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "15min", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "30min", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "60min", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "24h", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "48h", seq == "proseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "0min", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "15min", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "30min", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "60min", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "24h", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "48h", seq == "cageseq"))$enhancer_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "0min", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "15min", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "30min", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "60min", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "24h", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "48h", seq == "proseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "0min", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "15min", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "30min", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "60min", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "24h", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "48h", seq == "cageseq"))$promoter_density) / 0.0009765625 # 67 % of what would be expected by chance
mean(chiapet_connections$promoter_density) / 0.0009765625 # 31 % of what would be expected by chance
mean((chiapet_connections %>% filter(hemin == "eClusters"))$enhancer_density) / 0.0009765625 # 82 % of what would be expected by chance

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
  labs(title = "Distance vs Enhancer REX-motif Density",
       x = "Enhancer REX-motif Density",
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
  labs(title = "Distance vs Promoter Density",
       x = "Promoter Density",
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
    hjust = 1.1, vjust = 1.1
  ) +
  labs(
    title = "Promoter-Enhancer Distance vs REX Enhancer Density",
    x = "Enhancer REX-motif Density",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
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
  labs(
    title = "Promoter-Enhancer Distance vs REX Promoter Density",
    x = "Promoter REX-motif Density",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  facet_wrap(~ group, nrow = 1) +
  theme_minimal()

## STATISTICAL TESTING

## REX motif ABSENCE or PRESENCE comparison

# REX ABSENCE or PRESENCE in ENHANCERS, for HEMIN CAGE SEQ DATA
plot1 <- ggplot(chiapet_connections %>%
         filter(seq == "cageseq") %>%
         mutate(REX_presence = ifelse(enhancer_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Hemin, CAGE-seq",
    x = "REX Motif in Enhancer",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", format.pval(t.test(log10(distance) ~ ifelse(enhancer_REX == 0, 'Absent', 'Present'), 
                                                     data = chiapet_connections %>% filter(seq == "cageseq"))$p.value, 
                                              digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in ENHANCERS, for HEMIN PRO SEQ DATA
plot2 <- ggplot(chiapet_connections %>%
         filter(seq == "proseq") %>%
         filter(!grepl("Nhs|HS|both", feature) & hemin != "eClusters") %>%
         mutate(REX_presence = ifelse(enhancer_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Hemin, PRO-seq",
    x = "REX Motif in Enhancer",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(enhancer_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(seq == "proseq") %>%
                                               filter(!grepl("Nhs|HS|both", feature) & hemin != "eClusters"))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in PROMOTERS, for HEMIN CAGE SEQ DATA
plot3 <- ggplot(chiapet_connections %>%
         filter(seq == "cageseq") %>%
         mutate(REX_presence = ifelse(promoter_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Hemin, CAGE-seq",
    x = "REX Motif in Promoter",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", format.pval(t.test(log10(distance) ~ ifelse(promoter_REX == 0, 'Absent', 'Present'), 
                                                     data = chiapet_connections %>% filter(seq == "cageseq"))$p.value, 
                                              digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in PROMOTERS, for HEMIN PRO SEQ DATA
plot4 <- ggplot(chiapet_connections %>%
         filter(seq == "proseq") %>%
         filter(!grepl("Nhs|HS|both", feature) & hemin != "eClusters") %>%
         mutate(REX_presence = ifelse(promoter_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Hemin, PRO-seq",
    x = "REX Motif in Promoter",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(promoter_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(seq == "proseq") %>%
                                               filter(!grepl("Nhs|HS|both", feature) & hemin != "eClusters"))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in ENHANCERS, for HEAT SHOCK DATA
plot5 <- ggplot(chiapet_connections %>%
         filter(feature %in% c("both", "Nhs", "HS")) %>%
         mutate(REX_presence = ifelse(enhancer_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Heat shock",
    x = "REX Motif in Enhancer",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(enhancer_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(feature %in% c("both", "Nhs", "HS")))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in PROMOTERS, for HEAT SHOCK DATA
plot6 <- ggplot(chiapet_connections %>%
         filter(feature %in% c("both", "Nhs", "HS")) %>%
         mutate(REX_presence = ifelse(promoter_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "Heat shock",
    x = "REX Motif in Promoter",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(promoter_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(feature %in% c("both", "Nhs", "HS")))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in ENHANCERS, for ENHANCER CLUSTER DATA
plot7 <- ggplot(chiapet_connections %>%
         filter(hemin == "eClusters") %>%
         mutate(REX_presence = ifelse(enhancer_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "eClusters",
    x = "REX Motif in Enhancer Cluster",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(enhancer_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(hemin == "eClusters"))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# REX ABSENCE or PRESENCE in PROMOTERS, for ENHANCER CLUSTER DATA
plot8 <- ggplot(chiapet_connections %>%
         filter(hemin == "eClusters") %>%
         mutate(REX_presence = ifelse(promoter_REX == 0, "Absent", "Present")), 
       aes(x = REX_presence, y = log10(distance))) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  labs(
    title = "eClusters",
    x = "REX Motif in Promoter",
    y = "log10(Promoter-Enhancer Distance)"
  ) +
  theme_classic() +
  annotate("text", 
           x = 1.5, y = max(log10(chiapet_connections$distance), na.rm = TRUE) * 0.95, 
           label = paste0("p = ", 
                          format.pval(t.test(log10(distance) ~ ifelse(promoter_REX == 0, "Absent", "Present"),
                                             data = chiapet_connections %>% filter(hemin == "eClusters"))$p.value, 
                                      digits = 3, eps = .001)), size = 5)

# Arrange the two plots side by side
grid.arrange(plot1, plot3, plot2, plot4, plot5, plot6, plot7, plot8, nrow = 2)

# Write to a .txt file as a comma-separated list, to be used for GO enrichment
write.table(chiapet_connections %>% filter(hemin == "eClusters") %>%
              pull(gene), file = "eClusters_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(feature == "Nhs") %>%
              pull(gene), file = "proseq_Nhs_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(feature == "HS") %>%
              pull(gene), file = "proseq_HS_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(feature == "both") %>%
              pull(gene), file = "proseq_both_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "0min", seq == "proseq") %>%
              pull(gene), file = "proseq_0min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "15min", seq == "proseq") %>%
              pull(gene), file = "proseq_15min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "30min", seq == "proseq") %>%
              pull(gene), file = "proseq_30min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "60min", seq == "proseq") %>%
              pull(gene), file = "proseq_60min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "24h", seq == "proseq") %>%
              pull(gene), file = "proseq_24h_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "48h", seq == "proseq") %>%
              pull(gene), file = "proseq_48h_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "0min", seq == "cageseq") %>%
              pull(gene), file = "cageseq_0min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "15min", seq == "cageseq") %>%
              pull(gene), file = "cageseq_15min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "30min", seq == "cageseq") %>%
              pull(gene), file = "cageseq_30min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "60min", seq == "cageseq") %>%
              pull(gene), file = "cageseq_60min_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "24h", seq == "cageseq") %>%
              pull(gene), file = "cageseq_24h_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% filter(hemin == "48h", seq == "cageseq") %>%
              pull(gene), file = "cageseq_48h_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
write.table(chiapet_connections %>% pull(gene), file = "~/Desktop/GO_reference_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

# Write gene names to CSV files to be used for GO analysis, comparing REX presence or absence
write.table(chiapet_connections %>%
              filter(hemin == "eClusters" & enhancer_REX == 0) %>%
              pull(gene), file = "~/Desktop/eClusters_enhancer_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(hemin == "eClusters" & enhancer_REX != 0) %>%
              pull(gene), file = "~/Desktop/eClusters_enhancer_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(hemin == "eClusters" & promoter_REX == 0) %>%
              pull(gene), file = "~/Desktop/eClusters_promoter_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(hemin == "eClusters" & promoter_REX != 0) %>%
              pull(gene), file = "~/Desktop/eClusters_promoter_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(feature %in% c("Nhs", "HS", "both") & enhancer_REX == 0) %>%
              pull(gene), file = "~/Desktop/heatshock_enhancer_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(feature %in% c("Nhs", "HS", "both") & enhancer_REX != 0) %>%
              pull(gene), file = "~/Desktop/heatshock_enhancer_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(feature %in% c("Nhs", "HS", "both") & promoter_REX == 0) %>%
              pull(gene), file = "~/Desktop/heatshock_promoter_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(feature %in% c("Nhs", "HS", "both") & promoter_REX != 0) %>%
              pull(gene), file = "~/Desktop/heatshock_promoter_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(seq == "proseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & enhancer_REX == 0) %>%
              pull(gene), file = "~/Desktop/hemin_proseq_enhancer_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # SIGNIFICANT FOUND
write.table(chiapet_connections %>%
              filter(seq == "proseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & enhancer_REX != 0) %>%
              pull(gene), file = "~/Desktop/hemin_proseq_enhancer_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # SIGNIFICANT FOUND
write.table(chiapet_connections %>%
              filter(seq == "proseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & promoter_REX == 0) %>%
              pull(gene), file = "~/Desktop/hemin_proseq_promoter_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # SIGNIFICANT FOUND
write.table(chiapet_connections %>%
              filter(seq == "proseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & promoter_REX != 0) %>%
              pull(gene), file = "~/Desktop/hemin_proseq_promoter_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(seq == "cageseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & enhancer_REX == 0) %>%
              pull(gene), file = "~/Desktop/hemin_cageseq_enhancer_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")  # SIGNIFICANT FOUND
write.table(chiapet_connections %>%
              filter(seq == "cageseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & enhancer_REX != 0) %>%
              pull(gene), file = "~/Desktop/hemin_cageseq_enhancer_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
write.table(chiapet_connections %>%
              filter(seq == "cageseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & promoter_REX == 0) %>%
              pull(gene), file = "~/Desktop/hemin_cageseq_promoter_REX_absent_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # SIGNIFICANT FOUND
write.table(chiapet_connections %>%
              filter(seq == "cageseq" & hemin %in% c("0min", "15min", "30min", "60min", "24h", "48h") & promoter_REX != 0) %>%
              pull(gene), file = "~/Desktop/hemin_cageseq_promoter_REX_present_genes.csv", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",") # No significant found
