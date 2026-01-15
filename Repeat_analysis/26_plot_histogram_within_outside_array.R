#!/usr/bin/env Rscript

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
# i = 4
print(i)

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))


setwd("/home/pwlodzimierz/Rhynchospora/git_Rhynchospora")
source("./aux_fun.R")
replace_existing_analysis = FALSE

data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))


data_directories <- data_directories[grep(paste(c("Rrugosa.chr.fasta",
                                                  "Ralba.chr.fasta",
                                                  "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta",
                                                  "Rcephalotes.hap2.chr.fasta", 
                                                  "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta",
                                                  "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta",
                                                  "R_pubera_ref_2n10_v2.chr.fasta",
                                                  "Rbreviuscula.hap1.chr.fasta", 
                                                  "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                  "Rcolorata.hap1.chr.fasta",
                                                  "Rtenerrima.chr.fasta",
                                                  # "Rhync_filiformis_5519G.hap1.chr.fasta",
                                                  "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                  "Rtenuis.hap1.chr.fasta",
                                                  "Rhync_riparia_5519C.hap1.chr.fasta",
                                                  "R_barbata_hap1_chrs.fa",
                                                  "R_holoschoenoides_hap1_chr.fasta",
                                                  "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta"), collapse = "|"), data_directories)]


for(i in seq_along(data_directories)) {
  print(paste0("A: ", i, " / ", length(data_directories)))
  
  setwd(data_directories[i])
  print(data_directories[i])
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); next}
  repeats = read.csv(file = repeat_file)
  repeats <- repeats[grep("Tyba", repeats$class),]
  
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  
  assembly_name_short = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  assembly_name_full = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]
  
  chr_no_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
  chr_no_sizes <- chr_no_sizes[chr_no_sizes$assembly.name == assembly_name_full, ]
  chromosomes <- chr_no_sizes$chromosome.name
  chromosomes <- sort(chromosomes)
  
  islands_genome_file <- list.files(pattern = "_islands_genome", full.names = TRUE)[1]
  if(length(islands_genome_file) != 1) {print(paste0(i, " no islands_genome_file!")); next}
  islands_genome <- read.csv(file = islands_genome_file)
  
  
  
  
  array_pairs_to_sample_for_similarity <- 1500
  repeats_to_sample_per_array <- 10
  
  
  repeats$arrayID = 0
  cat("====================================================================================================")
  for(j in 1 : nrow(islands_genome)) {
    if(j %in% round((1: 100) * nrow(islands_genome)/100)) cat("#")
    
    repeats$arrayID[repeats$seqID == islands_genome$chromosome[j] & repeats$end > islands_genome$start[j] & repeats$start < islands_genome$end[j]] = j
    
  }
  cat("\n")
  
  
  ### similarity within and between islands_genome ###
  {
    sample_arrays <- sample(repeats$arrayID, size = array_pairs_to_sample_for_similarity, replace = TRUE)
    
    
    similarity_data_frame <- data.frame(method = rep("same_array", array_pairs_to_sample_for_similarity),
                                        array_1 = sample_arrays,
                                        array_2 = sample_arrays,
                                        similarity_score = rep(0, array_pairs_to_sample_for_similarity))
    
    for(j in 1 : array_pairs_to_sample_for_similarity) {
      cat("array", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
      repeats_A <- repeats[repeats$arrayID == sample_arrays[j], ]
      repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
      dist_matrix <- adist(repeats_A_sample$sequence)
      similarity_data_frame$similarity_score[j] <- 100-100*mean(dist_matrix[upper.tri(dist_matrix)])/mean(repeats_A_sample$width)
    }
    
    
    # Get indices by chromosome
    chr_groups <- split(1:nrow(islands_genome), islands_genome$chromosome)
    
    # Initialize storage
    same_chr_pairs <- matrix(NA, nrow = array_pairs_to_sample_for_similarity, ncol = 2)
    pair_count <- 0
    
    # Keep sampling until we have array_pairs_to_sample_for_similarity valid pairs
    while (pair_count < array_pairs_to_sample_for_similarity) {
      # Randomly choose a chromosome group
      chr <- sample(names(chr_groups), 1)
      idxs <- chr_groups[[chr]]
      
      # Only sample if there are at least 2 elements
      if (length(idxs) >= 2) {
        pair <- sample(idxs, 2)
        pair_count <- pair_count + 1
        same_chr_pairs[pair_count, ] <- pair
      }
    }
    
    similarity_data_frame <- rbind(similarity_data_frame, data.frame(method = rep("different_array_same_chr", array_pairs_to_sample_for_similarity),
                                                                     array_1 = rep(0, array_pairs_to_sample_for_similarity),
                                                                     array_2 = rep(0, array_pairs_to_sample_for_similarity),
                                                                     similarity_score = rep(0, array_pairs_to_sample_for_similarity)))
    
    for(j in 1 : array_pairs_to_sample_for_similarity) {
      cat("chr", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
      repeats_A <- repeats[repeats$arrayID == same_chr_pairs[j,1], ]
      repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
      repeats_B <- repeats[repeats$arrayID == same_chr_pairs[j,2], ]
      repeats_B_sample <- repeats_B[sample(1 : nrow(repeats_B), replace = TRUE, size = repeats_to_sample_per_array),]
      similarity_data_frame$similarity_score[array_pairs_to_sample_for_similarity + j] = 
        100-100*mean(adist(repeats_A_sample$sequence, repeats_B_sample$sequence))/mean(c(repeats_A_sample$width, repeats_B_sample$width))
      similarity_data_frame$array_1[array_pairs_to_sample_for_similarity + j] = same_chr_pairs[j,1]
      similarity_data_frame$array_2[array_pairs_to_sample_for_similarity + j] = same_chr_pairs[j,2]
    }
    
    
    
    
    diff_chr_pairs <- matrix(NA, nrow = array_pairs_to_sample_for_similarity, ncol = 2)
    pair_count <- 0
    n <- nrow(islands_genome)
    
    while (pair_count < array_pairs_to_sample_for_similarity) {
      pair <- sample(1:n, 2)
      chr1 <- islands_genome$chromosome[pair[1]]
      chr2 <- islands_genome$chromosome[pair[2]]
      
      if (chr1 != chr2) {
        pair_count <- pair_count + 1
        diff_chr_pairs[pair_count, ] <- pair
      }
    }
    
    similarity_data_frame <- rbind(similarity_data_frame, data.frame(method = rep("different_array_different_chr", array_pairs_to_sample_for_similarity),
                                                                     array_1 = rep(0, array_pairs_to_sample_for_similarity),
                                                                     array_2 = rep(0, array_pairs_to_sample_for_similarity),
                                                                     similarity_score = rep(0, array_pairs_to_sample_for_similarity)))
    
    for(j in 1 : array_pairs_to_sample_for_similarity) {
      cat("genome", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
      repeats_A <- repeats[repeats$arrayID == diff_chr_pairs[j,1], ]
      repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
      repeats_B <- repeats[repeats$arrayID == diff_chr_pairs[j,2], ]
      repeats_B_sample <- repeats_B[sample(1 : nrow(repeats_B), replace = TRUE, size = repeats_to_sample_per_array),]
      similarity_data_frame$similarity_score[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = 
        100-100*mean(adist(repeats_A_sample$sequence, repeats_B_sample$sequence))/mean(c(repeats_A_sample$width, repeats_B_sample$width))
      similarity_data_frame$array_1[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = diff_chr_pairs[j,1]
      similarity_data_frame$array_2[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = diff_chr_pairs[j,2]
    }
    print("A1")
    xmin = 0
    xmax = 100
    colors <- c("#0072B2", "#E69F00", "#009E73")
    for(ymax in c(100,150,200,250,300,350,400,500)) {
      
      similarity_data_frame$similarity_score[similarity_data_frame$similarity_score < 0] = 0
      similarity_data_frame$similarity_score[similarity_data_frame$similarity_score > 100] = 100
      
      pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/26_histograms_similarity_within_outside/", assembly_name_short, "_", ymax, "_histograms_similarity_holocentric_within_outside.pdf"), width = 12, height = 12)
      par(mfrow = c(2,1))
      hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
           breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), border = "#0072B290", col = "#0072B290",
           xlab = "similarity distribution of repeats", main = paste0(assembly_name_short, " histogram of Tyba pairwise similarity within:"))
      hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
           breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#E69F0090", col = "#E69F0090")
      hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
           breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#009E7390", col = "#009E7390")
      legend(x = xmax * 0.9, y = ymax*0.95, legend = c("arrays", "chromosomes", "genomes"), fill = c("#0072B290", "#E69F0090", "#009E7390"))
      
      # boxplot(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
      #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
      #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
      #         names = c("arrays", "chromosomes", "genomes"),
      #         main = "per chromosome mean pairwise similarity of repeats within:",
      #         col = c("#0072B290", "#E69F0090","#009E7390", ylim = c(0,100)))
      
      arrays2      <- similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"]
      chromosomes <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"]
      genomes     <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"]
      
      # Plot the boxplot
      boxplot(arrays2, chromosomes, genomes,
              names = c("arrays", "chromosomes", "genomes"),
              main = "bocplot of similarities of Tyba repeats within:",
              col = c("#0072B290", "#E69F0090", "#009E7390"),
              ylim = c(0, 140))
      
      # Perform pairwise t-tests
      p1 <- t.test(arrays2, chromosomes)$p.value
      p2 <- t.test(arrays2, genomes)$p.value
      p3 <- t.test(chromosomes, genomes)$p.value
      
      # Add asterisks based on significance levels
      # Positioning settings
      y_max <- max(c(arrays2, chromosomes, genomes), na.rm = TRUE)
      offset <- 5
      
      # Function to determine asterisk level
      get_asterisks <- function(p) {
        if (p < 0.001) return("***")
        if (p < 0.01)  return("**")
        if (p < 0.05)  return("*")
        return("n.s.")  # not significant
      }
      
      # Draw lines and asterisks
      segments(1, y_max + offset, 2, y_max + offset)
      text(1.5, y_max + offset + 2, get_asterisks(p1))
      
      segments(1, y_max + 2*offset, 3, y_max + 2*offset)
      text(2, y_max + 2*offset + 2, get_asterisks(p2))
      
      segments(2, y_max + 3*offset, 3, y_max + 3*offset)
      text(2.5, y_max + 3*offset + 2, get_asterisks(p3))
      dev.off()
      
      
    }
    
  } ### similarity within and between arrays ###
  
}

















