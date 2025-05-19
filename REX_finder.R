## REX-motif finding script
# This script returns a .bed file of each REX motif found in a given .fasta file. 

# Download hg19 chromosome FASTA files from this link: 
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/

# Run in terminal as: 
# Rscript ~/Documents/GitHub/vihervaara/REX_finder.R ~/Documents/Vihervaara/hg19/hg19_chr/chr1.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr2.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr3.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr4.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr5.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr6.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr7.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr8.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr9.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr10.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr11.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr12.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr13.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr14.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr15.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr16.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr17.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr18.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr19.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr20.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr21.fa ~/Documents/Vihervaara/hg19/hg19_chr/chr22.fa ~/Documents/Vihervaara/hg19/hg19_chr/chrX.fa ~/Documents/Vihervaara/hg19/hg19_chr/chrY.fa ~/Documents/Vihervaara/hg19/hg19_chr/rex_new.bed

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Avoid scientific notation globally
options(scipen = 999)

# Check number of arguments
if (length(args) < 2) {
  stop("Usage: Rscript script.R <fasta1.fa> <fasta2.fa> ... <output.bed>")
}

# Last argument is the output BED file
bed_filename <- args[length(args)]
fasta_files <- args[-length(args)]

# REX motifs
motifs <- c("taatta", "caatta", "taattg")

# Function to find motif positions in one .fasta file
find_REX <- function(fasta, motifs) {
  fasta_lines <- readLines(fasta, warn = FALSE)
  chr <- sub("^>(\\S+).*", "\\1", fasta_lines[1])
  sequence <- tolower(paste(fasta_lines[-1], collapse = ""))
  len <- nchar(sequence)
  
  if (len < 6) return(NULL)
  
  starts <- 1:(len - 5)
  windows <- substring(sequence, starts, starts + 5)
  match_idx <- which(windows %in% motifs)
  
  cat(sprintf("Chromosome %s: %d motifs found\n", chr, length(match_idx)))
  
  if (length(match_idx) == 0) return(NULL)
  
  data.frame(
    chr = chr,
    start = match_idx - 1,
    end = match_idx + 5,
    stringsAsFactors = FALSE
  )
}

# Apply to all .fasta files and combine results
rex_list <- lapply(fasta_files, find_REX, motifs = motifs)
all_rex_df <- do.call(rbind, rex_list)

# Report total number of motifs
cat(sprintf("\nTotal REX motifs found: %d\n", nrow(all_rex_df)))

# Write .bed file
write.table(all_rex_df, bed_filename, quote = FALSE, row.names = FALSE, 
            col.names = FALSE, sep = "\t")

