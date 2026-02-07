#!/usr/bin/env Rscript

# DADA2 Pipeline for 16S rRNA Analysis
# Generates ASV table and representative sequences from paired-end FASTQ files
#
# Input: data/原始测序数据/*.fq.gz (paired-end reads)
# Output:
#   - figure3/dada2_output/rep_seqs.fasta (representative sequences)
#   - figure3/dada2_output/asv_table.csv (ASV abundance table)

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
})

# --- Configuration ---
fastq_dir <- "data/原始测序数据"
out_dir <- "figure3/dada2_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Primer trimming parameters (adjust if needed)
# V3-V4 region typically uses 341F/806R primers
truncLen_F <- 220  # Truncate forward reads
truncLen_R <- 200  # Truncate reverse reads

# --- Get file paths ---
fnFs <- sort(list.files(fastq_dir, pattern = "_R1.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(fastq_dir, pattern = "_R2.fq.gz", full.names = TRUE))

# Extract sample names
sample_names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)
# Clean up sample names to match existing data
sample_names <- gsub("ZTPSN24MW098-Control_T", "RDYS_1-", sample_names)
sample_names <- gsub("ZTPSN24MW099-Treat1_T", "RDYS_2-", sample_names)
sample_names <- gsub("ZTPSN24MW100-Treat2_T", "RDYS_3-", sample_names)
sample_names <- gsub("ZTPSN24MW101-Treat3_T", "RDYS_4-", sample_names)
sample_names <- gsub("ZTPSN24MW102-Treat4_T", "RDYS_5-", sample_names)
sample_names <- gsub("_W$", "", sample_names)

message("Found ", length(fnFs), " paired-end samples:")
print(data.frame(Sample = sample_names, Forward = basename(fnFs), Reverse = basename(fnRs)))

# --- Quality filtering ---
message("\n=== Step 1: Quality filtering ===")
filtFs <- file.path(out_dir, "filtered", paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(out_dir, "filtered", paste0(sample_names, "_R_filt.fastq.gz"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

filt_out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  truncLen = c(truncLen_F, truncLen_R),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)
message("Filtering complete. Reads passed:")
print(filt_out)

# --- Learn error rates ---
message("\n=== Step 2: Learning error rates ===")
# Use more bases and randomize across all samples for better error estimation
errF <- learnErrors(filtFs, multithread = TRUE, verbose = TRUE, nbases = 2e8, MAX_CONSIST = 20, randomize = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE, verbose = TRUE, nbases = 2e8, MAX_CONSIST = 20, randomize = TRUE)

# Save error plots
pdf(file.path(out_dir, "error_rates.pdf"), width = 10, height = 8)
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

# --- Dereplicate ---
message("\n=== Step 3: Dereplication ===")
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# --- Sample inference (DADA2 core algorithm) ---
message("\n=== Step 4: Sample inference (DADA2) ===")
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, verbose = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, verbose = TRUE)

# --- Merge paired reads ---
message("\n=== Step 5: Merging paired reads ===")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# --- Construct ASV table ---
message("\n=== Step 6: Constructing ASV table ===")
seqtab <- makeSequenceTable(mergers)
message("ASV table dimensions: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")

# --- Remove chimeras ---
message("\n=== Step 7: Removing chimeras ===")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
message("After chimera removal: ", ncol(seqtab_nochim), " ASVs")
message("Reads retained: ", round(sum(seqtab_nochim) / sum(seqtab) * 100, 1), "%")

# --- Track reads through pipeline ---
getN <- function(x) sum(getUniques(x))
track <- cbind(
  filt_out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample_names
message("\nRead tracking summary:")
print(track)
write.csv(track, file.path(out_dir, "read_tracking.csv"))

# --- Export representative sequences ---
message("\n=== Exporting outputs ===")
asv_seqs <- colnames(seqtab_nochim)
asv_headers <- paste0("ASV_", seq_along(asv_seqs))

# Write FASTA
fasta_out <- file.path(out_dir, "rep_seqs.fasta")
writeLines(
  paste0(">", asv_headers, "\n", asv_seqs),
  fasta_out
)
message("Representative sequences written to: ", fasta_out)

# Write ASV table with ASV names instead of sequences
asv_tab <- t(seqtab_nochim)
rownames(asv_tab) <- asv_headers
write.csv(asv_tab, file.path(out_dir, "asv_table.csv"))
message("ASV table written to: ", file.path(out_dir, "asv_table.csv"))

# Save sequence-to-ASV mapping
seq_map <- data.frame(ASV_ID = asv_headers, Sequence = asv_seqs)
write.csv(seq_map, file.path(out_dir, "asv_sequences.csv"), row.names = FALSE)

message("\n=== DADA2 pipeline complete! ===")
message("Total ASVs: ", length(asv_seqs))
message("Output directory: ", out_dir)
