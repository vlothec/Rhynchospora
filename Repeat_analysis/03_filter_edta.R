#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))
suppressMessages(library(msa))
suppressMessages(library(GenomicRanges))

setwd("/home/pwlodzimierz/Rhynchospora/git_Rhynchospora")
source("./aux_fun.R")
replace_existing_analysis = FALSE

edta_overlap_max_perc = 0.8 # if an edta annotation overlaps more than this with repeat ARRAY coordinates, it is discarded


data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

# for(i in seq_along(data_directories)) 
{
  print(paste0(i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, "No repeats!")); setwd(".."); next}
  
  array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
  if(length(array_file) != 1) {print(paste0(i, " no arrays!")); setwd(".."); next}
  
  classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
  if(length(classes_file) != 1) {print(paste0(i, " no classes!")); setwd(".."); next}
  
  edta_file = list.files(pattern = "_edta_modified.csv", full.names = TRUE)
  if(length(edta_file) != 1) {print(paste0(i, " no edta!")); setwd(".."); next}
  
  assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  assembly_file = grep(strsplit(strsplit(repeat_file, split = "_repeats")[[1]][1], split = "/")[[1]][2], assembly_files)
  print(assembly_file)
  print(assembly_files[assembly_file])
  
  if(!replace_existing_analysis) {
    if(file.exists(paste0(assembly_name, "_edta_filtered.csv"))) quit(save = "no", status = 1)
  }
  
  if(length(assembly_file) == 1) {
    assembly = read.fasta(assembly_files[assembly_file], seqtype = "DNA", forceDNAtolower = TRUE)
  } else {
    print(" No assembly found, cannot proceed")
    next
  }
  
  repeats = read.csv(file = repeat_file)
  arrays = read.csv(file = array_file)
  classes = read.csv(file = classes_file)
  edta = read.csv(file = edta_file)
  
  if(assembly_name == "R_holoschoenoides_hap1_chr") { # 3
    for(j in 1 : length(unique(edta$V1))) {
      edta$V1[edta$V1 == unique(edta$V1)[j]] = paste0(unique(edta$V1)[j], "1")
    } 
  }
  if(assembly_name == "Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr") { # 3
    for(j in 1 : length(unique(edta$V1))) {
      edta$V1[edta$V1 == 1] = "HiC_scaffold_1_h2"
      edta$V1[edta$V1 == 2] = "HiC_scaffold_2_h2"
      edta$V1[edta$V1 == 3] = "HiC_scaffold_3_h2"
      edta$V1[edta$V1 == 4] = "HiC_scaffold_4_h2"
      edta$V1[edta$V1 == 5] = "HiC_scaffold_5_h2"
      edta$V1[edta$V1 == 6] = "HiC_scaffold_6_h2"
    } 
  }
  if(assembly_name == "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr") { # 47
    edta$V1[edta$V1 == 1] = "HiC_scaffold_1"
    edta$V1[edta$V1 == 2] = "HiC_scaffold_2"
    edta$V1[edta$V1 == 3] = "HiC_scaffold_3"
    edta$V1[edta$V1 == 4] = "HiC_scaffold_4"
    edta$V1[edta$V1 == 5] = "HiC_scaffold_5"
  }
  
  if(assembly_name == "R_pubera_ref_2n10_v2.chr") { # 5
    arrays$seqID <- unlist(lapply(arrays$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  }
  
  # repeats = repeats[repeats$new_class %in% classes$]
  
  arrays$overlapping_bp = 0
  edta$overlapping_bp = 0
  
  edta_filtered_total = NULL
  if(!dir.exists("sequence_repeat_edta_plots")) dir.create("sequence_repeat_edta_plots")
  # for each sequence
  cen_rep_classes_total = NULL
  cen_edta_v3_total = NULL
  cen_edta_class_total = NULL
  for (j in seq_along(assembly)) {  
    sequence_arrays = arrays[arrays$seqID == names(assembly)[j], ]
    sequence_edta = edta[edta$V1 == names(assembly)[j], ]
    sequence_repeats = repeats[repeats$seqID == names(assembly)[j], ]
    
    if(nrow(sequence_arrays) + nrow(sequence_edta) == 0) next
    
    ### Filter EDTA if needed ==================================================
    if(nrow(sequence_arrays) != 0 & nrow(sequence_edta) != 0) {
      gr1 <- with(sequence_arrays, GRanges(names(assembly)[j], IRanges(start, end)))
      gr2 <- with(sequence_edta, GRanges(names(assembly)[j], IRanges(V4, V5)))
      
      overlaps <- as.data.frame(findOverlaps(gr1, gr2))
      
      for(k in seq_len(nrow(overlaps))) {
        print(paste(names(assembly)[j], k, "/", nrow(overlaps)))
        overlap_bp = width(pintersect(gr1[overlaps$queryHits[k]], gr2[overlaps$subjectHits[k]]))
        sequence_arrays$overlapping_bp[overlaps$queryHits[k]] = sequence_arrays$overlapping_bp[overlaps$queryHits[k]] + overlap_bp
        sequence_edta$overlapping_bp[overlaps$subjectHits[k]] = sequence_edta$overlapping_bp[overlaps$subjectHits[k]] + overlap_bp
      }
      sequence_arrays$width = sequence_arrays$end - sequence_arrays$start + 1
      sequence_edta$width = sequence_edta$V5 - sequence_edta$V4 + 1
      sequence_edta$overlapping_percentage = sequence_edta$overlapping_bp / sequence_edta$width
      sequence_arrays$overlapping_percentage = sequence_arrays$overlapping_bp / sequence_arrays$width
      
      edta_filtered = sequence_edta[sequence_edta$overlapping_percentage <= edta_overlap_max_perc,]
      
      edta_filtered_total = rbind(edta_filtered_total, sequence_edta[sequence_edta$overlapping_percentage <= edta_overlap_max_perc,])
      
    } else {
      edta_filtered = sequence_edta
    }
  }
  edta <- edta_filtered_total
  str(edta)
  export_gff(annotations.data.frame = edta[,-1], output = ".", 
             file.name = paste0(assembly_name, "_edta_filtered"), seqid = 1, source = 2, type = 3, 
             start = 4, end = 5, strand = 7, score = ".")
  
  write.csv(edta[,-1], paste0(assembly_name, "_edta_filtered.csv"))
  
}