# CAGE Project Venn Diagram

library(eulerr)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# TOTAL ENHANCERS COMPARISON

# Define the datasets
venn_data <- data.frame(
  comparison = c("CAGE vs PRO (Heat shock)", "CAGE vs PRO (Hemin)", "PRO (HS) vs PRO (Hemin)"),
  group1 = c("CAGE-seq", "CAGE-seq", "PRO (HS)"),
  group2 = c("PRO (HS)", "PRO (Hemin)", "PRO (Hemin)"),
  group1_only = c(11697, 10485, 15036),
  group2_only = c(26158, 26503, 17969),
  overlap = c(3484, 4427, 14577),
  stringsAsFactors = FALSE
)

# Assign consistent grayscale colors to groups
group_colors <- c(
  "CAGE-seq" = "white",   # light grey (gray90)
  "PRO (HS)" = "gray60",   # mid grey (gray60)
  "PRO (Hemin)" = "gray30" # dark grey (gray30)
)

# Blend function for two hex colors
blend_colors <- function(col1, col2) {
  rgb_row <- (col2rgb(col1) + col2rgb(col2)) / 2
  rgb(rgb_row[1], rgb_row[2], rgb_row[3], maxColorValue = 255)
}

# Generate Euler diagrams
plots <- lapply(seq_len(nrow(venn_data)), function(i) {
  row <- venn_data[i, ]
  
  group1 <- row$group1
  group2 <- row$group2
  overlap_name <- paste0(group1, "&", group2)
  
  # Blend color for the overlap
  blended <- blend_colors(group_colors[group1], group_colors[group2])
  
  # Define region sizes
  areas <- c(
    setNames(row$group1_only + row$overlap, group1),
    setNames(row$group2_only + row$overlap, group2),
    setNames(row$overlap, overlap_name)
  )
  
  # Assign corresponding fills
  fills <- list(
    fill = c(group_colors[group1], group_colors[group2], blended),
    alpha = 0.9
  )
  
  plot(euler(areas),
       labels = NULL,
       fills = fills,
       quantities = list(cex = 0.5),
       main = row$comparison)
})

# Build legend
legend <- legendGrob(
  labels = names(group_colors),
  pch = 22,
  gp = gpar(fill = unname(group_colors), col = "black"),
  nrow = 1
)

# Arrange diagrams and legend
grid.arrange(do.call(arrangeGrob, c(plots, ncol = length(plots))),
             legend,
             nrow = 2,
             heights = c(10, 1))

### HEMIN COMPARISON

# Create a tidy data frame of set sizes
venn_data <- data.frame(
  time = factor(c("0min", "15min", "30min", "60min", "24h", "48h"),
                levels = c("0min", "15min", "30min", "60min", "24h", "48h")),
  CAGE_seq_only = c(1144, 1578, 1558, 1837, 1514, 1384),
  PRO_seq_only = c(26572, 24521, 25219, 25453, 25765, 26525),
  overlap = c(1478, 1945, 2114, 2234, 1751, 1702)
)

# Generate Euler diagrams for each time point
plots <- lapply(seq_len(nrow(venn_data)), function(i) {
  row <- venn_data[i, ]
  fit <- euler(c("CAGE-seq" = row$CAGE_seq_only + row$overlap,
                 "PRO-seq" = row$PRO_seq_only + row$overlap,
                 "CAGE-seq&PRO-seq" = row$overlap))
  plot(fit, labels = NULL, quantities = list(cex = 0.5), main = paste("Hemin", row$time))
})

# Add a legend
legend <- legendGrob(
  labels = c("CAGE-seq", "PRO-seq"),
  pch = 22,
  gp = gpar(fill = c("white", "grey"), col = "black"),
  nrow = 1
)

# Display legend under the grid
grid.arrange(do.call(arrangeGrob, c(plots, ncol = length(plots))), legend, nrow = 2, heights = c(10, 1))
