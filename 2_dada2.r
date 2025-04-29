# DADA2 pipeline for reads resulting from cutadapt. The code is very minimally 
# modified from the tutorial at https://benjjneb.github.io/dada2/tutorial.html

# Load packages
library(dada2)
#packageVersion("dada2")
#‘1.26.0’
library(dplyr)

# Establish the file path to the directory containing the demultiplexed samples
path <- "demultiplex_cutadapt/trimmed_reads"
list.files(path) # Check that all the files are there 

# Forward and reverse fastq filenames have the format:  
# SAMPLENAME-R1.fastq and SAMPLENAME-R2.fastq
# Use this fact to sort them into two groups
fnFs <- sort(list.files(path, pattern = ".1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".2.fastq.gz", full.names = TRUE)) 

# Extract sample names (i.e., exclude the forward/reverse identifier), 
# assuming filenames have format SAMPLENAME-Rn.fastq and SAMPLENAME
# does not include any "-"
sample_names <- sapply(strsplit(basename(fnFs), "-R"), `[`, 1)

any(duplicated(sample_names)) # FALSE, so we can proceed

# Grab four samples, the reads of which will be examined in terms of 
# their quality profiles
ids <- round(runif(4,2,length(sample_names)))

# Output the quality profiles for the forward and reverse reads of 
# those samples
dir.create('figures')
pdf("figures/cutadapt_final_qualprofiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()

# Create a directory and names for the filtered files
filtFs <- file.path("reads/cutadapt_filtered_reads", 
                    paste0(sample_names, "-F-filt.fastq"))
filtRs <- file.path("reads/cutadapt_filtered_reads", 
                    paste0(sample_names, "-R-filt.fastq"))
names(filtFs) <- sample_names
names(filtRs) <- sample_names

# Perform quality filtering 
# I tested a bunch of truncLen values, 260,215 - 230,190
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,180),
                    maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     compress = TRUE, multithread = FALSE, verbose = TRUE)
# Kept ~60% with these parameters

# Learn the error rates
# slow
errF <- learnErrors(filtFs, nbases = 1e+09, multithread = TRUE)
errR <- learnErrors(filtRs, nbases = 1e+09, multithread = TRUE)

# Plot the error rates
pdf("figures/cutadapt_error_plot_final_F.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf("figures/cutadapt_error_plot_final_R.pdf")
plotErrors(errR, nominalQ = TRUE)
dev.off()

# Use the filtered files and error rates to perform
# sample inference (without pooling)
dadaFs <- dada(filtFs, err = errF, pool = FALSE, multithread = TRUE)
dadaFs[[1]]
# dada-class: object describing DADA2 denoising results
# 441 sequence variants were inferred from 32342 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

dadaRs <- dada(filtRs, err = errR, pool = FALSE, multithread = TRUE)
dadaRs[[1]]
# dada-class: object describing DADA2 denoising results
# 407 sequence variants were inferred from 19735 input unique sequences.
# Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

# Merge the forward and reverse paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# Build the counts table
counts_raw <- makeSequenceTable(mergers) 

# Check that all of the sequence lengths are within the expected range
# this will vary based on the parameters you used for quality filtering
table(nchar(getSequences(counts_raw)))
# Output:
#  388   402   403   404   405   406   407   408   409   410   411   412   413 
#     2   121  1570 11460   843   337   330    52   237   212   599    92   474 
#   414   415   417   418   419   420   421   422   423   424   425   426   428 
#     5     1     2     7    32    11    30    29   305 10173   161     2   571 
#   429   430   431 
#  4038   750     1 
# Most ASVs are acceptable lengths, but the others must be removed
# lots of variability in length of V3-V4 amplicons!

counts_trimmed <- counts_raw[, nchar(colnames(counts_raw)) %in% seq(402, 431)]

# Filter chimeras
counts_nochim <- removeBimeraDenovo(counts_trimmed, method = "consensus", 
                                    multithread = TRUE, verbose = TRUE)
# Output:
# Identified 29939 bimeras out of 32445 input sequences.


sum(counts_nochim)/sum(counts_trimmed)
#[1]	0.7510661


# Assign taxonomy
tax_nochim <- assignTaxonomy(counts_nochim,
                             "silva_nr99_v138.2_toSpecies_trainset.fa", 
                             multithread = TRUE)

# Filter by taxonomy
tax_filtered <- as.data.frame(tax_nochim) %>%
  filter(!is.na(Kingdom)) %>%
  filter(Kingdom != "Eukaryota") %>%
  filter(Family != "Mitochondria") %>%
  filter(Order != "Chloroplast")
  
counts_filtered <- counts_nochim[, rownames(tax_filtered)]

# Assign readable names
tax_filtered <- tax_filtered %>%
  mutate(Sequence = rownames(tax_filtered))
rownames(tax_filtered) <- paste0("SV_", 1:nrow(tax_filtered))

any(colnames(counts_filtered) != tax_filtered$Sequence) # FALSE, so the ASVs are in the same order

colnames(counts_filtered) <- paste0("SV_", 1:ncol(counts_filtered))

# Construct a table to summarize the removal of reads throughout the pipeline
getN <- function(x) {sum(getUniques(x))}
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(counts_trimmed), rowSums(counts_nochim), rowSums(counts_filtered))
colnames(track) <- c("Input", "Filtered", "DenoisedF", "DenoisedR", "Merged", "Trimmed", "Non-Chimeric", "Tax-Filtered")
rownames(track) <- sample_names

# Output the counts, tax, and tracking tables
# make data directory
dir.create('data')
write.table(counts_filtered, 
            file = "data/cutadapt_counts2.txt", 
            sep = "\t", col.names = NA, quote = F)

write.table(tax_filtered, 
            file = "data/cutadapt_tax2.txt",
            sep = "\t", col.names = NA, quote = F)

write.table(track, 
            file = "data/cutadapt_track2.txt",
            sep = "\t", col.names = NA, quote = F)

#save.image("scripts/dada2_workspace.RData")