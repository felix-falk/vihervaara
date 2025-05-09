#!/usr/bin/env Rscript

# Load required packages
suppressWarnings(suppressMessages(library(dplyr)))

# === Parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: Rscript intersect_bed_by_column_with_name.R file1.bed file2.bed output.bed")
}

file1 <- as.character(args[1])
file2 <- as.character(args[2])
output_file <- as.character(args[3])

# === Helper: Read and clean a BED file
read_bed <- function(file) {
  lines <- readLines(file)
  split <- strsplit(lines, "\t")
  clean <- split[sapply(split, length) > 0]
  
  if (length(clean) == 0) {
    stop(paste("Error: File", file, "has no readable rows."))
  }
  
  max_cols <- max(sapply(clean, length))
  df <- as.data.frame(do.call(rbind, lapply(clean, function(x) {
    length(x) <- max_cols
    return(x)
  })), stringsAsFactors = FALSE)
  
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  return(df)
}

# === Find the column index containing "Name_" in any row
find_name_column <- function(df) {
  match_col <- which(sapply(df, function(col) any(grepl("Name_", col))))
  if (length(match_col) != 1) {
    stop("❌ Error: Could not uniquely identify a single column containing 'Name_'.")
  }
  return(colnames(df)[match_col])
}

# === Read BED files
bed1 <- read_bed(file1)
bed2 <- read_bed(file2)

# === Identify the column containing "Name_" in each
col1 <- find_name_column(bed1)
col2 <- find_name_column(bed2)

# === Rename both columns to a common name for merging
bed1 <- bed1 %>% rename(NameCol = all_of(col1))
bed2 <- bed2 %>% rename(NameCol = all_of(col2))

# === Merge based on the "Name_" identifier
merged <- merge(bed1, bed2, by = "NameCol", all = FALSE)

# === Write output
write.table(merged, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cat(paste("Matching rows (by 'Name_') written to:", output_file, "\n"))

# bed1_df <- read.delim("~/Documents/Vihervaara/hg19/chiapet/enhancer_ChIAPETLeft_overlap.bed", header=FALSE)
# bed2_df <- read.delim("~/Documents/Vihervaara/hg19/chiapet/enhancer_ChIAPETRight_overlap.bed", header=FALSE)
# merged <- merge(bed1_df, bed2_df, by = "V11", all = FALSE)
