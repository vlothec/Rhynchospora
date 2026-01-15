#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))
suppressMessages(library(msa))

setwd("/home/pwlodzimierz/Rhynchospora/git_Rhynchospora")
source("./aux_fun.R")

replace_existing_analysis = FALSE

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
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  
  assembly_name = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]
  
  edta_file = list.files(path = ".", pattern = "TEanno", recursive = F, full.names = T)[grep(assembly_name, list.files(path = ".", pattern = "TEanno", recursive = F, full.names = T))]
  
  if(!replace_existing_analysis) {
    if(file.exists(paste0(assembly_name, "_edta_modified.csv"))) quit(save = "no", status = 1)
  }
  if(length(edta_file) == 1) {
    edta_og = read.table(edta_file, header = FALSE, sep = "\t", skip = 6, comment.char = "#", stringsAsFactors = FALSE)
  } else {
    print(" No edta found, cannot proceed")
    next
  }
  ### Reformat EDTA gff table ==================================================
  v10 = lapply(edta_og$V9, function(X) strsplit(X, split = ";")[[1]])
  new_cols = lapply(v10, function(X) unlist(lapply(X, function(x) strsplit(x, split = "=")[[1]][1])))
  new_cols = unique(unlist(new_cols))
  edta_full = edta_og[, 1:8]
  for (j in seq_along(new_cols)) {
    cat(j, "/", length(new_cols), "\n")
    
    new_data = lapply(edta_og$V9, function(X) strsplit(X, split = paste0(new_cols[j], "=", collapse = ""))[[1]][2])
    new_data = unlist(lapply(new_data, function(X) strsplit(X, split = ";")[[1]][1]))
    
    edta_full = cbind(edta_full, new_data)
    names(edta_full)[ncol(edta_full)] = new_cols[j]
  }
  
  edta_full$oldV3 = edta_full$V3
  
  edta_repeat_region = edta_full[edta_full$V3 == "repeat_region", ]
  edta_non_rep_reg = edta_full[edta_full$V3 != "repeat_region", ]
  
  names = unique(edta_repeat_region$Name)
  names = sort(names)
  
  names = data.frame(names)
  names$short = unlist(lapply(names$names, function(X) {
    pos = grep("[^A-Za-z]", strsplit(X, split = "")[[1]])
    if(length(pos) == 0 ) return(X)
    return(substr(X, 1, min(pos) - 1))
    } ))
  
  edta_repeat_region$oldV3 = edta_repeat_region$V3
  edta_repeat_region$V3 = unlist(lapply(seq_len(nrow(edta_repeat_region)), function(X) names$short[names$names == edta_repeat_region$Name[X]]))
  
  edta_modified = rbind(edta_non_rep_reg, edta_repeat_region)
  
  if(assembly_name == "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta") {
    edta_modified$V1[edta_modified$V1 == 1] = "HiC_scaffold_1_h1"
    edta_modified$V1[edta_modified$V1 == 2] = "HiC_scaffold_2_h1"
    edta_modified$V1[edta_modified$V1 == 3] = "HiC_scaffold_3_h1"
    edta_modified$V1[edta_modified$V1 == 4] = "HiC_scaffold_4_h1"
    edta_modified$V1[edta_modified$V1 == 5] = "HiC_scaffold_5_h1"
  }
  
  if(assembly_name == "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta") {
    edta_modified$V1[edta_modified$V1 == 1] = "HiC_scaffold_1"
    edta_modified$V1[edta_modified$V1 == 2] = "HiC_scaffold_2"
    edta_modified$V1[edta_modified$V1 == 3] = "HiC_scaffold_3"
    edta_modified$V1[edta_modified$V1 == 4] = "HiC_scaffold_4"
    edta_modified$V1[edta_modified$V1 == 5] = "HiC_scaffold_5"
  }
  
  if(assembly_name == "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta") {
    edta_modified$V1[edta_modified$V1 == 1] = "HiC_scaffold_1_h1"
    edta_modified$V1[edta_modified$V1 == 2] = "HiC_scaffold_2_h1"
    edta_modified$V1[edta_modified$V1 == 3] = "HiC_scaffold_3_h1"
    edta_modified$V1[edta_modified$V1 == 4] = "HiC_scaffold_4_h1"
    edta_modified$V1[edta_modified$V1 == 5] = "HiC_scaffold_5_h1"
    edta_modified$V1[edta_modified$V1 == 6] = "HiC_scaffold_6_h1"
  }
  
  if(assembly_name == "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta") {
    edta_modified$V1[edta_modified$V1 == 1] = "HiC_scaffold_1"
    edta_modified$V1[edta_modified$V1 == 2] = "HiC_scaffold_2"
    edta_modified$V1[edta_modified$V1 == 3] = "HiC_scaffold_3"
    edta_modified$V1[edta_modified$V1 == 4] = "HiC_scaffold_4"
    edta_modified$V1[edta_modified$V1 == 5] = "HiC_scaffold_5"
    edta_modified$V1[edta_modified$V1 == 6] = "HiC_scaffold_6"
    edta_modified$V1[edta_modified$V1 == 7] = "HiC_scaffold_7"
  }
  
  if(assembly_name == "R_holoschoenoides_hap1_chr") {
    for(j in 1 : length(unique(edta_modified$V1))) {
      edta_modified$V1[edta_modified$V1 == unique(edta_modified$V1)[j]] = paste0(unique(edta_modified$V1)[j], "1")
    } 
  }
  
  if(assembly_name == "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr") {
    edta_modified$V1[edta_modified$V1 == 1] = "HiC_scaffold_1"
    edta_modified$V1[edta_modified$V1 == 2] = "HiC_scaffold_2"
    edta_modified$V1[edta_modified$V1 == 3] = "HiC_scaffold_3"
    edta_modified$V1[edta_modified$V1 == 4] = "HiC_scaffold_4"
    edta_modified$V1[edta_modified$V1 == 5] = "HiC_scaffold_5"
  }
  
  
  
  export_gff(annotations.data.frame = edta_modified, output = ".", 
             file.name = paste0(assembly_name, "_edta_modified"), seqid = 1, source = 2, type = 3, 
             start = 4, end = 5, strand = 7, score = ".")
  
  write.csv(edta_modified, paste0(assembly_name, "_edta_modified.csv"))
  
  
}












