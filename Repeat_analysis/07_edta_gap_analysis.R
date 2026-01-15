#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

suppressMessages(library(seqinr))
suppressMessages(library(msa))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

max_space_between_repeats_to_ignore = 50000
min_gap_size = 250


if(F) { # go to line 79
  repeats <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_repeats_filtered.csv")
  edta <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta_edta_filtered.csv")
  classes <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_classes_merged_filtered.csv")
  setwd("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data")
  source("C:/Users/Piotr Włodzimierz/Desktop/ToL/temp_data/aux_fun.R")
  
  max_space_between_repeats_to_ignore = 50000
  min_gap_size = 250
  
  suppressMessages(library(seqinr))
  suppressMessages(library(msa))
  
  assembly_name <- "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr"
  
  i=1
}


data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 30
print(i)

{
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  print(data_directories[i])
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); quit(save = "no")}
  
  no_edta = FALSE
  edta_file = list.files(pattern = "_edta_filtered.csv", full.names = TRUE)
  if(length(edta_file) != 1) {print(paste0(i, " no edta!")); no_edta = TRUE}
  
  assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  
  print(assembly_name)
  
  repeats = read.csv(file = repeat_file)
  
  if(!no_edta) edta = read.csv(file = edta_file)
  
  chromosomes <- unique(repeats$seqID)
  
  classes_all <- unique(repeats$class) 
  
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ,]
  
  chromosomes <- unique(repeats$seqID)
  
  islands_genome <- NULL
  
  for(j in seq_along(chromosomes)) {
    
    cat("/n", i, assembly_name, j, "/", length(chromosomes), chromosomes[j])
    
    tyba_chr <- tyba[tyba$seqID == chromosomes[j], ]
    tyba_chr <- tyba_chr[order(tyba_chr$start, decreasing = FALSE), ]
    if(!no_edta) edta_chr <- edta[edta$V1 == chromosomes[j], ]
    
    gaps <- tyba_chr$start[2 : nrow(tyba_chr)] - tyba_chr$end[1 : (nrow(tyba_chr) - 1)]
    
    # sort(tyba_chr$start) == tyba_chr$start
    # if(grepl("alba", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("austrobrasiliensis", assembly_name)) max_space_between_repeats_to_ignore <- 100000
    # if(grepl("barbata", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("breviuscula", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("cephalotes", assembly_name)) max_space_between_repeats_to_ignore <- 30000
    # if(grepl("ciliata", assembly_name)) max_space_between_repeats_to_ignore <- 80000
    # if(grepl("colorata", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("corymbosa", assembly_name)) max_space_between_repeats_to_ignore <- 30000
    # if(grepl("gaudichaudii", assembly_name)) max_space_between_repeats_to_ignore <- 100000
    # if(grepl("holoschoenoides", assembly_name)) max_space_between_repeats_to_ignore <- 70000
    # if(grepl("nervosa", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("pubera", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("radicans", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("ridleyi", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("riparia", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    # if(grepl("rugosa", assembly_name)) max_space_between_repeats_to_ignore <- 70000
    # if(grepl("tenerrima", assembly_name)) max_space_between_repeats_to_ignore <- 60000
    # if(grepl("tenuis", assembly_name)) max_space_between_repeats_to_ignore <- 100000
    # if(grepl("watsonii", assembly_name)) max_space_between_repeats_to_ignore <- 50000
    
    
    gaps[gaps < max_space_between_repeats_to_ignore] = 0
    
    
    island_ends <- which(gaps > 0)
    island_starts <- c(1, (island_ends + 1))
    island_ends <- c(island_ends, nrow(tyba_chr))
    
    if(length(island_starts) == 1) {
      islands_chr <- data.frame(genome = rep(assembly_name, length(island_starts)),
                                chromosome = rep(chromosomes[j], length(island_starts)),
                                start = tyba_chr$start[island_starts],
                                end = tyba_chr$end[island_ends],
                                distance_to_next = NA,
                                total_bp = tyba_chr$end[island_ends] - tyba_chr$start[island_starts] + 1,
                                tyba_no = island_ends - island_starts + 1,
                                tyba_total_bp = rep(NA, length(island_starts)),
                                tyba_mean_size = rep(NA, length(island_starts)),
                                tyba_mean_size_sd_perc = rep(NA, length(island_starts)),
                                tyba_consensus = rep("", length(island_starts)),
                                tyba_internal_ED_to_consensus_size_normalised = rep(NA, length(island_starts)),
                                no_of_gaps = rep(NA, length(island_starts)),
                                tatal_gap_bp = rep(NA, length(island_starts)),
                                no_of_gaps_with_any_TE_coverage = rep(NA, length(island_starts)),
                                longest_gap_bp = rep(NA, length(island_starts)),
                                gaps_TE_coverage.perc = rep(NA, length(island_starts)),
                                gaps_LTR_coverage.perc = rep(NA, length(island_starts)),
                                gaps_gypsy_coverage.perc = rep(NA, length(island_starts)),
                                gaps_copia_coverage.perc = rep(NA, length(island_starts)),
                                gaps_helitron_coverage.perc = rep(NA, length(island_starts)),
                                gaps_TIR_coverage.perc = rep(NA, length(island_starts)),
                                region_to_next_TE.bp = rep(NA, length(island_starts)))
    } else {
      islands_chr <- data.frame(genome = rep(assembly_name, length(island_starts)),
                                chromosome = rep(chromosomes[j], length(island_starts)),
                                start = tyba_chr$start[island_starts],
                                end = tyba_chr$end[island_ends],
                                distance_to_next = c(tyba_chr$start[island_ends[2:length(island_ends)]] - tyba_chr$end[island_starts[1:(length(island_starts) - 1)]] - 1, NA),
                                total_bp = tyba_chr$end[island_ends] - tyba_chr$start[island_starts] + 1,
                                tyba_mean_size = rep(NA, length(island_starts)),
                                tyba_mean_size_sd_perc = rep(NA, length(island_starts)),
                                tyba_no = island_ends - island_starts + 1,
                                tyba_total_bp = rep(NA, length(island_starts)),
                                tyba_consensus = rep("", length(island_starts)),
                                tyba_internal_ED_to_consensus_size_normalised = rep(NA, length(island_starts)),
                                no_of_gaps = rep(NA, length(island_starts)),
                                tatal_gap_bp = rep(NA, length(island_starts)),
                                no_of_gaps_with_any_TE_coverage = rep(NA, length(island_starts)),
                                longest_gap_bp = rep(NA, length(island_starts)),
                                gaps_TE_coverage.perc = rep(NA, length(island_starts)),
                                gaps_LTR_coverage.perc = rep(NA, length(island_starts)),
                                gaps_gypsy_coverage.perc = rep(NA, length(island_starts)),
                                gaps_copia_coverage.perc = rep(NA, length(island_starts)),
                                gaps_helitron_coverage.perc = rep(NA, length(island_starts)),
                                gaps_TIR_coverage.perc = rep(NA, length(island_starts)),
                                region_to_next_TE.bp = rep(NA, length(island_starts)))
      
    }
    islands_chr <- islands_chr[islands_chr$tyba_no >= 3,]
    
    cat(" Tyba:", nrow(tyba_chr), "Islands:", nrow(islands_chr), "")
    
    
    for(k in seq_len(nrow(islands_chr))) {
      
      if(k %in% round(nrow(islands_chr) / 10 * (1:9))) cat("#")
      
      tyba_island <- tyba_chr[island_starts[k] : island_ends[k], ]
      
      
      islands_chr$tyba_mean_size[k] <- mean(tyba_island$width)
      islands_chr$tyba_mean_size_sd_perc[k] <- 100 * sd(tyba_island$width) / islands_chr$tyba_mean_size[k]
      islands_chr$tyba_total_bp[k] = sum(tyba_island$width)
      
      if(T) { # consensus and adist
        sequences <- tyba_island$sequence
        if(length(sequences) > 100) sequences = sample(sequences, 100)
        if(length(sequences) > 1) {
          sequences_bs <- Biostrings::DNAStringSet(sequences)
          names(sequences_bs) <- seq_along(sequences)
          msa_message <- capture.output({alignment <- msa::msa(sequences_bs, method = "ClustalOmega", type = "dna",
                                                               verbose = FALSE)})
          
          if (msa_message != "using Gonnet") {
            print(msa_message)
          }
          
          consensus <- consensusMatrix(alignment)
          coverages <- unlist(lapply(seq_len(ncol(consensus)), function(X) sum(consensus[(1:4),X])))
          consensus <- consensus[(1:4), sort(order(coverages, decreasing = T)[1 : round(median(nchar(sequences)))])]
          consensus <- names(unlist(lapply(seq_len(ncol(consensus)), function(X)  which.max(consensus[,X]))))
          consensus <- paste(consensus, collapse = "")
          
          islands_chr$tyba_consensus[k] = consensus
          islands_chr$tyba_internal_ED_to_consensus_size_normalised[k] = 100 * mean(adist(consensus, tyba_island$sequence, ignore.case = TRUE, costs = c(1,1,1))[1, ]) / nchar(consensus)
          
          remove(consensus, sequences)
        }
        
      }
      if(!no_edta) {
        if(k < nrow(islands_chr)) {
          overlaps <- overlapsRanges(IRanges(edta_chr$V4, edta_chr$V5), IRanges(islands_chr$end[k], islands_chr$start[k + 1]))
          overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
          islands_chr$region_to_next_TE.bp[k] = overlap_coords
        }
      }
      
      
      small_gaps <- tyba_island$start[2 : nrow(tyba_island)] - tyba_island$end[1 : (nrow(tyba_island) - 1)] - 1
      
      islands_chr$no_of_gaps[k] = sum(small_gaps >= min_gap_size)
      islands_chr$longest_gap_bp[k] = max(small_gaps)
      
      if(!is.na(islands_chr$no_of_gaps[k])) {
        if(islands_chr$no_of_gaps[k] > 0) {
          gaps_starts <- tyba_island$end[1 : (nrow(tyba_island) - 1)] + 1
          gaps_ends <-  tyba_island$start[2 : nrow(tyba_island)] + 1
          gaps_size <- small_gaps
          
          gaps_starts <- gaps_starts[gaps_size >= min_gap_size]
          gaps_ends <- gaps_ends[gaps_size >= min_gap_size]
          gaps_size <- gaps_size[gaps_size >= min_gap_size]
          islands_chr$tatal_gap_bp[k] <- sum(gaps_size)
          
          if(!no_edta) {
            
            islands_chr$no_of_gaps_with_any_TE_coverage[k] <- sum(overlapsAny(IRanges(gaps_starts, gaps_ends), IRanges(edta_chr$V4, edta_chr$V5)))
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4, edta_chr$V5), IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_TE_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4[grepl("(?i)^(?!.*non[-]?LTR).*LTR|long_terminal_repeat", edta_chr$V3, ignore.case = T, perl = T)], edta_chr$V5[grepl("(?i)^(?!.*non[-]?LTR).*LTR|long_terminal_repeat", edta_chr$V3, ignore.case = T, perl = T)]), 
                                       IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_LTR_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4[grepl("copia", edta_chr$V3, ignore.case = T)], edta_chr$V5[grepl("copia", edta_chr$V3, ignore.case = T)]), 
                                       IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_copia_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4[grepl("gypsy", edta_chr$V3, ignore.case = T)], edta_chr$V5[grepl("gypsy", edta_chr$V3, ignore.case = T)]), 
                                       IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_gypsy_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4[grepl("helitron", edta_chr$V3, ignore.case = T)], edta_chr$V5[grepl("helitron", edta_chr$V3, ignore.case = T)]), 
                                       IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_helitron_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            overlaps <- overlapsRanges(IRanges(edta_chr$V4[grepl("TIR", edta_chr$V3, ignore.case = T)], edta_chr$V5[grepl("TIR", edta_chr$V3, ignore.case = T)]), 
                                       IRanges(gaps_starts, gaps_ends))
            overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
            islands_chr$gaps_TIR_coverage.perc[k] = 100 *  overlap_coords / sum(gaps_size)
            
            
          }
          
          
        }
      }
      
      remove(small_gaps)
      
    }
    
    islands_genome <- rbind(islands_genome, islands_chr)
    
  }
  
  if(!no_edta) write.csv(x = islands_genome, file = paste0(assembly_name, "_islands_genome.csv"), row.names = FALSE) else write.csv(x = islands_genome, file = paste0(assembly_name, "_islands_genome_no_edta.csv"), row.names = FALSE)
  if(!no_edta) write.csv(x = islands_genome, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/islands/", assembly_name, "_islands_genome.csv"), row.names = FALSE) else write.csv(x = islands_genome, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/islands/", assembly_name, "_islands_genome_no_edta.csv"), row.names = FALSE)
  
}






















