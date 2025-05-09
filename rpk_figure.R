# Mean and SD enhancer length for each hemin condition

library(ggplot2)
library(dplyr)
library(broom)

overlap_0min <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_0_min.bed", header=FALSE)
overlap_15min <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_15_min.bed", header=FALSE)
overlap_30min <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_30_min.bed", header=FALSE)
overlap_60min <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_60_min.bed", header=FALSE)
overlap_24h <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_24h.bed", header=FALSE)
overlap_48h <- read.delim("~/Documents/Vihervaara/hg19/overlapping_enhancers_He_48h.bed", header=FALSE)

overlap_0min <- overlap_0min %>% mutate(V3 - V2) %>% mutate(time = "0 min", seq = "both")
overlap_15min <- overlap_15min %>% mutate(V3 - V2) %>% mutate(time = "15 min", seq = "both")
overlap_30min <- overlap_30min %>% mutate(V3 - V2) %>% mutate(time = "30 min", seq = "both")
overlap_60min <- overlap_60min %>% mutate(V3 - V2) %>% mutate(time = "60 min", seq = "both")
overlap_24h <- overlap_24h %>% mutate(V3 - V2) %>% mutate(time = "24 h", seq = "both")
overlap_48h <- overlap_48h %>% mutate(V3 - V2) %>% mutate(time = "48 h", seq = "both")

colnames(overlap_0min) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")
colnames(overlap_15min) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")
colnames(overlap_30min) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")
colnames(overlap_60min) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")
colnames(overlap_24h) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")
colnames(overlap_48h) <- c("chr_CAGE_seq", "eStart_CAGE_seq", "eEnd_CAGE_seq", 
                            "RPK_CAGE_seq", "chr_PRO_seq", "eStart_PRO_seq", 
                            "eEnd_PRO_seq", "RPK_PRO_seq", 
                            "length_CAGE_seq", "time", "seq")

overlapping_enhancers <- bind_rows(overlap_0min, overlap_15min, overlap_30min, overlap_60min, overlap_24h, overlap_48h)

# Prepare the data
df <- overlapping_enhancers %>%
  mutate(
    time = factor(time, levels = c("0 min", "15 min", "30 min", "60 min", "24 h", "48 h")),
    log_CAGE = log10(RPK_CAGE_seq),
    log_PRO = log10(RPK_PRO_seq)
  )

# Fit models and calculate R-squared AND AOV p-value per time group
stats_labels <- df %>%
  group_by(time) %>%
  do({
    model = lm(log_PRO ~ log_CAGE, data = .)
    r2 = summary(model)$r.squared
    p_value = anova(model)$`Pr(>F)`[1]
    
    # Format p-value
    p_value_formatted <- ifelse(p_value < 0.001, "< 0.001", signif(p_value, 3))
    
    data.frame(r2 = r2, p_value = p_value_formatted)
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("R² = ", round(r2, 3), "\n", "p = ", p_value),  # combine R² and p-value
    x = Inf,
    y = Inf
  )

# Now plot
ggplot(df, aes(x = log_CAGE, y = log_PRO)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  geom_text(
    data = stats_labels,
    aes(x = x, y = y, label = label),
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE
  ) +
  facet_wrap(~ time, nrow = 2) +
  theme_minimal() +
  labs(
    x = "log10(CAGE-seq RPK)",
    y = "log10(PRO-seq RPK)")
