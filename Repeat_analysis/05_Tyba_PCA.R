#!/usr/bin/env Rscript

# Install required packages if needed
# install.packages(c("BiocManager", "ggplot2", "seqinr"))
# BiocManager::install("Biostrings")

library(seqinr)      # For reading FASTA files
library(Biostrings)  # For sequence manipulation
library(ggplot2)     # For visualization
library(irlba)

sample_repeats = 1000

colors_vector <- c("#E69F00", 
                   "#D55E00", 
                   "#D95F02", 
                   "#FC8D62", 
                   "#E5C494", 
                   "#F0E442", 
                   "#FFD92F", 
                   "#A6D854", 
                   "#66C2A5", 
                   "#009E73", 
                   "#1B9E77", 
                   "#66A61E", 
                   "#CC79A7", 
                   "#E78AC3", 
                   "#E7298A", 
                   "#B3B3CC", 
                   "#8DA0CB", 
                   "#7570B3", 
                   "#56B4E9", 
                   "#0072B2")

custom_labels <- c(
  expression(italic("R. rugosa")), 
  expression(italic("R. alba")), 
  expression(italic("R. gaudichaudii")), 
  expression(italic("R. cephalotes")), 
  expression(italic("R. ridleyi")), 
  expression(italic("R. watsonii")), 
  expression(italic("R. radicans")), 
  expression(italic("R. pubera")), 
  expression(italic("R. breviscula")), 
  expression(italic("R. nervosa")), 
  expression(italic("R. ciliata")), 
  expression(italic("R. colorata")), 
  expression(italic("R. tenerrima")), 
  expression(italic("R. filiformis")), 
  expression(italic("R. austrobrasiliensis")), 
  expression(italic("R. tenuis")), 
  expression(italic("R. riparia")), 
  expression(italic("R. barbata")), 
  expression(italic("R. holoschoenoides")), 
  expression(italic("R. corymbosa"))
)

# Define point shapes (pch values) to improve distinction
pch_vector <- c(15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18)

custom_order <- c(
  "rugosa", 
  "alba", 
  "gaudichaudii", 
  "cephalotes", 
  "ridleyi", 
  "watsonii", 
  "radicans", 
  "pubera", 
  "breviuscula", 
  "nervosa", 
  "ciliata", 
  "colorata", 
  "tenerrima", 
  "filiformis", 
  "austrobrasiliensis", 
  "tenuis", 
  "riparia", 
  "barbata", 
  "holoschoenoides", 
  "corymbosa"
)


# Step 0: Load the repeats
# Sample if necessary

data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)

data_directories <- data_directories[grepl(pattern = paste("R_barbata_hap1_chrs.fa",
                                                           "R_holoschoenoides_hap1_chr.fasta",
                                                           "R_pubera_ref_2n10_v2.chr.fasta",
                                                           "Ralba.chr.fasta",
                                                           "Rbreviuscula.hap1.chr.fasta", 
                                                           "Rcephalotes.hap1.chr.fasta", 
                                                           "Rcolorata.hap1.chr.fasta", 
                                                           "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 
                                                           "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 
                                                           "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta", 
                                                           "Rhync_filiformis_5519G.hap1.chr.fasta", 
                                                           "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_riparia_5519C.hap1.chr.fasta", 
                                                           "Rtenuis.hap1.chr.fasta", 
                                                           "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta", 
                                                           "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta", 
                                                           "Rrugosa.chr.fasta", 
                                                           "Rtenerrima.chr.fasta",
                                                           sep = "|"), x = data_directories)]

tyba_all <- NULL
for(i in seq_along(data_directories)) {
  cat("reading in assembly from", data_directories[i], "\n")
  setwd(data_directories[i])
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  repeats = read.csv(file = repeat_file)
  sample_no <- sample_repeats
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  if(!is.na(sample_repeats) & sample_no > nrow(tyba)) sample_no <- nrow(tyba)
  if(!is.na(sample_repeats)) tyba <- tyba[sample(seq_len(nrow(tyba)), sample_no, replace = FALSE), ]
  tyba$assembly <- strsplit(data_directories[i], split = "/")[[1]][6]
  tyba$width <- tyba$end - tyba$start + 1
  tyba_all <- rbind(tyba_all, tyba[,c(1,8,14)])
}


# Step 1: Align repeats
setwd("/home/pwlodzimierz/Rhynchospora/upload_data/PCA")
write.fasta(sequences = as.list(tyba_all$sequence), names = c(paste(seq_len(nrow(tyba_all)))), file.out = paste0(sample_repeats, "_tyba_repeats.fasta"))

system(paste0("mafft --retree 2 --thread -1 --inputorder ", as.character(sample_repeats),  "_tyba_repeats.fasta > ", sample_repeats, "_aligned.fasta"))

fasta_file <- paste0(sample_repeats, "_aligned.fasta")
sequences <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)

# Convert to character vector
seq_vector <- sapply(sequences, function(x) toupper(as.character(x)))

# Check sequence length (should all be equal after MAFFT alignment)
seq_lengths <- nchar(seq_vector)
summary(seq_lengths)
if (length(unique(seq_lengths)) > 1) {
  stop("Sequences have varying lengths after alignment. Ensure MAFFT output is properly aligned.")
}
seq_length <- unique(seq_lengths)
cat("Aligned sequence length:", seq_length, "bp\n")

# Step 2: Feature Extraction
# Toggle between k-mer and positional features
feature_type <- "positional"  # Change to "kmer" to use k-mer frequencies

if (feature_type == "kmer") {
  # K-mer frequencies (4-mers)
  k <- 4
  kmers <- oligonucleotideFrequency(DNAStringSet(seq_vector), width = k, as.prob = FALSE)
  feature_matrix <- as.matrix(kmers)  # 1M rows x 256 columns (4^4)
  cat("Using 4-mer frequencies with", ncol(feature_matrix), "features\n")
} else if (feature_type == "positional") {
  # Positional nucleotide frequencies
  seq_matrix <- do.call(rbind, strsplit(seq_vector, ""))
  feature_matrix <- matrix(0, nrow = length(seq_vector), ncol = seq_length * 4)
  for (i in 1:seq_length) {
    # if(i%/%100 == 0) cat(i,"\n")
    idx <- (i - 1) * 4 + 1:4
    feature_matrix[, idx[1]] <- seq_matrix[, i] == "A"
    feature_matrix[, idx[2]] <- seq_matrix[, i] == "C"
    feature_matrix[, idx[3]] <- seq_matrix[, i] == "G"
    feature_matrix[, idx[4]] <- seq_matrix[, i] == "T"
    # Gaps (-) are left as 0,0,0,0
  }
  cat("Using positional frequencies with", ncol(feature_matrix), "features\n")
} else {
  stop("Invalid feature_type. Use 'kmer' or 'positional'.")
}

# Remove constant/zero-variance columns
variance <- apply(feature_matrix, 2, var)
constant_cols <- which(variance == 0)
if (length(constant_cols) > 0) {
  cat("Removing", length(constant_cols), "constant/zero-variance columns from feature_matrix\n")
  feature_matrix <- feature_matrix[, variance > 0]
}

# Check if any columns remain
if (ncol(feature_matrix) == 0) {
  stop("All columns have zero variance. Try a different feature extraction method.")
}

# Step 3: Run PCA
variance <- apply(feature_matrix, 2, var)
keep_cols <- which(variance > quantile(variance, 0.1))  # Keep top 90% variance
feature_matrix_reduced <- feature_matrix[, keep_cols]
system.time(pca_result <- prcomp(feature_matrix_reduced, center = TRUE, scale. = TRUE))
# system.time(pca_result <- prcomp(feature_matrix, center = TRUE, scale. = TRUE))
pca_df <- as.data.frame(pca_result$x[, 1:2])  # Take first 2 PCs
colnames(pca_df) <- c("PC1", "PC2")

# Proportion of variance explained
var_explained <- summary(pca_result)$importance[2, 1:2] * 100
cat("Variance explained by PC1 and PC2:", var_explained, "\n")

# Step 4: Visualize
pca_df$Factor <- tyba_all$assembly

pca_df$Factor <- factor(pca_df$Factor, levels = c("Rrugosa.chr.fasta",
                                                  "Ralba.chr.fasta",
                                                  "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta",
                                                  "Rcephalotes.hap1.chr.fasta",
                                                  "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta",
                                                  "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta",
                                                  "R_pubera_ref_2n10_v2.chr.fasta",
                                                  "Rbreviuscula.hap1.chr.fasta",
                                                  "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                  "Rcolorata.hap1.chr.fasta",
                                                  "Rtenerrima.chr.fasta",
                                                  "Rhync_filiformis_5519G.hap1.chr.fasta",
                                                  "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                  "Rtenuis.hap1.chr.fasta",
                                                  "Rhync_riparia_5519C.hap1.chr.fasta",
                                                  "R_barbata_hap1_chrs.fa",
                                                  "R_holoschoenoides_hap1_chr.fasta",
                                                  "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta"))

write.csv(x = pca_df, file = paste0(sample_repeats, "_pca_df.csv"))
# pca_df <- read.csv(paste0(sample_repeats, "_pca_df.csv"))
# var_explained <- c(2.7,2.5)

pca_df_sorted = NULL
for(i in seq_along(custom_order)) {
  temp = pca_df[grep(custom_order[i], pca_df$Factor), ]
  temp$group = custom_order[i]  # Use the Factor or custom_order as the grouping variable
  pca_df_sorted <- rbind(pca_df_sorted, temp)
}

pca_df_sorted <- pca_df_sorted[sample(nrow(pca_df_sorted)), ]

pca_df_sorted$group <- factor(pca_df_sorted$group, levels = custom_order)

p <- ggplot(pca_df_sorted, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = colors_vector, labels = custom_labels) +  # Use colors_vector for colors
  scale_shape_manual(values = pch_vector, labels = custom_labels) +    # Use pch_vector for shapes
  labs(
    title = "PCA of Satellite Repeats",
    x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_set(theme_bw()) +
  xlim(c(-2.5, 2.5)) + ylim(c(-2, 1))


# Save the plot
ggsave(paste0(sample_repeats, "_pca_satellite_repeats_colored_borders.png"), p, width = 12, height = 10, dpi = 300)
ggsave(paste0(sample_repeats, "_pca_satellite_repeats_colored_borders.pdf"), p, width = 12, height = 10, dpi = 300)










