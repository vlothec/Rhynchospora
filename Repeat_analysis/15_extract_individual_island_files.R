#!/usr/bin/env Rscript

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

### load R libraries
library(seqinr)


### Load data ==================================================================
data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)

data_directories <- data_directories[grepl(paste(c("R_barbata_hap1_chrs.fa",
                                                   "R_holoschoenoides_hap1_chr.fasta",
                                                   "R_pubera_ref_2n10_v2.chr.fasta",
                                                   "Ralba.chr.fasta",
                                                   "Rbreviuscula.hap1.chr.fasta", 
                                                   "Rcephalotes.hap1.chr.fasta", 
                                                   "Rcolorata.hap1.chr.fasta",
                                                   "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                   "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                   "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta",
                                                   "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                   "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                   "Rhync_riparia_5519C.hap1.chr.fasta",
                                                   "Rhync_tenuis_ref.hap1.chr.fasta",
                                                   "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta",
                                                   "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta",
                                                   "Rradicans.016.FINAL.chr.fasta",
                                                   "Rrugosa.chr.fasta",
                                                   "Rtenerrima.chr.fasta"), collapse = "|"), data_directories)]

assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))


for(i in seq_along(data_directories)) {
  print(paste0("A: ", i, " / ", length(data_directories)))
  
  setwd(data_directories[i])
  print(data_directories[i])
  
  if(!dir.exists("./tyba_island_split")) dir.create("./tyba_island_split")
  if(!dir.exists("./tyba_island_split")) dir.create("./tyba_island_split_HORs")
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); quit(save = "no")}
  repeats = read.csv(file = repeat_file)
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  tyba$width <- tyba$end - tyba$start + 1
  tyba_name <- tyba$new_class[1]
  
  assembly_name_short = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  assembly_name_full = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]
  
  chr_no_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
  chr_no_sizes <- chr_no_sizes[chr_no_sizes$assembly.name == assembly_name_full, ]
  chromosomes <- chr_no_sizes$chromosome.name
  
  islands_genome_file <- list.files(pattern = "_islands_genome", full.names = TRUE)
  if(length(islands_genome_file) != 1) {print(paste0(i, " no islands_genome_file!")); quit(save = "no")}
  islands_genome <- read.csv(file = islands_genome_file)
  ### Load data ##################################################################
  
  
  
  chromosomes <- sort(chromosomes)
  
  dir.create(path = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/tyba_per_island/", assembly_name_short))
  setwd(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/tyba_per_island/", assembly_name_short))
  
  for(j in seq_along(chromosomes)) {
    
    tyba_chr <- tyba[tyba$seqID == chromosomes[j], ]
    islands_chr <- islands_genome[islands_genome$chromosome == chromosomes[j], ]
    
    for(k in seq_len(nrow(islands_chr))) {
      cat(assembly_name_short, j, "/", length(chromosomes), k, "/", nrow(islands_chr), "\n")
      
      tyba_island <- tyba_chr[(tyba_chr$end - 5) %in% ((islands_chr$start[k] - 200) : (islands_chr$end[k] + 200)), ]
      if(nrow(tyba_island) < 10) next
      
      tyba_island$seqID <- unlist(lapply(tyba_island$seqID, function(X) paste0(X, "_island", k)))
      
      write.fasta(sequences = as.list(tyba_island$sequence), 
                  names = paste0(tyba_island$seqID, rep("_", nrow(tyba_island)), 1 : nrow(tyba_island), sep = ""),  
                  file.out = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/tyba_per_island/", assembly_name_short, "/", assembly_name_short, "_chr_", j, "_isl_", k, ".fasta"))
      
      
    }
  }
}










































