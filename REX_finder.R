## REX-motif finding script
# This script returns a BED file of each REX motif found in a given FASTA file. 
# It takes about 5-15 minutes to run, depending on the size of the FASTA file. 

# Download hg19 chromosome FASTA files from this link: 
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/

# Download hs1 chromosome FASTA file from this link: 
# https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/

# Run in terminal as: 
# Rscript ~/Documents/GitHub/vihervaara/REX_finder.R ~/Documents/Vihervaara/hg19/hg19_chr/chr1.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr2.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr3.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr4.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr5.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr6.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr7.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr8.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr9.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr10.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr11.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr12.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr13.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr14.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr15.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr16.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr17.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr18.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr19.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr20.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr21.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr22.fa ~/Documents/Vihervaara/hg19/hg19_chr/chrX.fa ~/Documents/Vihervaara/hg19/hg19_chr/chrY.fa ~/Documents/Vihervaara/hg19/hg19_rex_bed/rex_new.bed

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)  

# Avoid scientific notation globally
options(scipen = 999)

if (length(args) < 2) {
  stop("Usage: Rscript script.R <fasta1.fa> <fasta2.fa> ... <output.bed>")
}

# Last argument is the output BED file
bed_filename <- args[length(args)]
fasta_files <- args[1:(length(args) - 1)]

# REX motifs
motifs <- c("taatta", "caatta", "taattg", "gaatta")

# Function to find motif positions in one FASTA
find_REX <- function(fasta, motifs) {
  fasta_lines <- readLines(fasta)
  header <- fasta_lines[1]
  chr <- sub("^>(\\S+).*", "\\1", header)
  
  sequence <- tolower(paste(fasta_lines[-1], collapse = ""))
  sequence_length <- nchar(sequence)
  
  positions <- 1:(sequence_length - 5)
  substrings <- substring(sequence, positions, positions + 5)
  
  match_indices <- positions[substrings %in% motifs]
  
  # Print summary to console
  cat(sprintf("Chromosome %s: %d motifs found\n", chr, length(match_indices)))
  
  return(data.frame(
    chr = rep(chr, length(match_indices)),
    start = match_indices - 1,
    end = match_indices + 5 
  ))
}

# Initialize empty data frame to collect all REX motifs
all_rex_df <- data.frame()
total_motifs <- 0  # Counter for total motif count

# Process each FASTA file
for (fasta in fasta_files) {
  rex_df <- find_REX(fasta, motifs)
  total_motifs <- total_motifs + nrow(rex_df)
  all_rex_df <- rbind(all_rex_df, rex_df)
}

# Report total number of motifs
cat(sprintf("\nTotal REX motifs found: %d\n", total_motifs))

# Write final combined BED file
write.table(all_rex_df, bed_filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

