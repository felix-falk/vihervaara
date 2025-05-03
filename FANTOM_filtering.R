# The script finds the datasets in the Human.sample_name2library_id file
# based on K562 cells and returns a dataset containing their annotations and
# filenames. 

# Creates a BED file of cell lines of interest from 

# Load required libraries
library(dplyr)
library(tidyr)
library(stringr) 

# Select experiment phrase of interest
count_matrix_path = "~/Documents/Vihervaara/hg19/cage-seq/human_permissive_enhancers_phase_1_and_2_expression_count_matrix.txt"
sample_library_path = "~/Documents/Vihervaara/hg19/cage-seq/Human.sample_name2library_id.txt"
bed_filename = "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562.bed"
RPK_matrix_out_path = "~/Documents/Vihervaara/hg19/cage-seq/CAGE_K562_RPK_matrix.txt"

# Fetch data set id data
sample_library <- read.delim(sample_library_path, header=FALSE)

# Return data sets of interest, those containing the phrase "K526"
filtered_sample_library <- sample_library %>% filter(grepl("K562", V1))

# Extract biological replicate information
filtered_sample_library$biol_rep <- str_extract(filtered_sample_library$V1, "biol_rep\\d+")
filtered_sample_library$biol_rep <- as.numeric(str_extract(filtered_sample_library$biol_rep, "\\d+"))

# Add row with hemin treatment time, and assign NA if hemin is missing
filtered_sample_library$hemin <- sub(".*hemin,\\s*([^,]+),.*", "\\1", filtered_sample_library$V1)
filtered_sample_library$hemin <- ifelse(grepl("hemin,\\s*[^,]+,", filtered_sample_library$V1), filtered_sample_library$hemin, NA)

# Filter out timepoints also found in PROseq data
filtered_sample_library <- filtered_sample_library %>%
  filter(hemin %in% c(NA, "00hr00min", "00hr15min", "00hr30min", "01hr00min", 
                      "24hr", "day02") | is.na(hemin))

### Create a mean count matrix and BED coordinate file for CAGE-seq data
# The output file is comparable to the PRO-seq HEMIN DATA file

# Import count matrix file
count_matrix <- read.delim(count_matrix_path, row.names = 1)

# Filter samples based on sample names list, remove rows with only 0:s
sample_names <- filtered_sample_library[[2]]
selected_count_matrix <- count_matrix %>% select(all_of(sample_names))
selected_count_matrix <- selected_count_matrix[rowSums(selected_count_matrix) > 0, ]

# Extract genomic coordinates
CAGE_coordinates <- rownames(selected_count_matrix) %>%
  str_extract_all("([^:-]+)") %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  setNames(c("chr", "eStart", "eEnd")) %>%
  mutate(across(c(eStart, eEnd), as.numeric))

# Calculate RPK (Reads Per Kilobase)
CAGE_RPK <- selected_count_matrix %>%
  mutate(length = CAGE_coordinates$eEnd - CAGE_coordinates$eStart) %>%
  mutate(across(everything(), ~ (. / length))) %>%
  select(-length)  # Remove length after normalization

# Compute mean RPK for time points
time_points <- list(
  min_0  = c("CNhs12458", "CNhs12684", "CNhs12786"),
  min_15 = c("CNhs12459", "CNhs12686", "CNhs12787"),
  min_30 = c("CNhs12460", "CNhs12687", "CNhs12788"),
  min_60 = c("CNhs12462", "CNhs12689", "CNhs12790"),
  hr_24  = c("CNhs12471", "CNhs12699", "CNhs12801"),
  hr_48 = c("CNhs12472", "CNhs12700", "CNhs12802")
)
for (tp in names(time_points)) {
  CAGE_RPK[[paste0(tp, "_mean_RPK")]] <- rowMeans(select(CAGE_RPK, all_of(time_points[[tp]])))
}

# Combine genomic coordinates and RPK data
CAGE_RPK <- cbind(CAGE_coordinates, CAGE_RPK)

# Remove unnecessary columns (adjust indices if needed)
CAGE_RPK <- CAGE_RPK[, -c(4:25)]

# Export final matrix
write.table(CAGE_RPK, RPK_matrix_out_path, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

write.table(CAGE_RPK %>% select(chr, eStart, eEnd, min_0_mean_RPK) %>% filter(min_0_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_0min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(CAGE_RPK %>% select(chr, eStart, eEnd, min_15_mean_RPK) %>% filter(min_15_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_15min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(CAGE_RPK %>% select(chr, eStart, eEnd, min_30_mean_RPK) %>% filter(min_30_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_30min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(CAGE_RPK %>% select(chr, eStart, eEnd, min_60_mean_RPK) %>% filter(min_60_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_60min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(CAGE_RPK %>% select(chr, eStart, eEnd, hr_24_mean_RPK) %>% filter(hr_24_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_1440min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(CAGE_RPK %>% select(chr, eStart, eEnd, hr_48_mean_RPK) %>% filter(hr_48_mean_RPK != 0), "~/Documents/Vihervaara/hg19/cage-seq/FANTOM5_K562_2880min.bed", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")



