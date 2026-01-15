#!/usr/bin/env Rscript
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 30
print(i)

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")


data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)

assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))



# for(i in 1 : 52) 
{
  chr_no_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ==================================================================
  setwd(data_directories[i])
  print(data_directories[i])
  
  islands_genome_file <- list.files(pattern = "_islands_genome", full.names = TRUE)
  if(length(islands_genome_file) == 0) {print(paste0(i, " no islands_genome_file!")); quit(save = "no")}
  if(length(islands_genome_file) > 1) {
    islands_genome_file <- islands_genome_file[grep("islands_genome.csv", islands_genome_file)]
  }
  islands_genome <- read.csv(file = islands_genome_file)
  
  if(!dir.exists("./tyba_island_split")) dir.create("./tyba_island_split")
  if(!dir.exists("./tyba_island_split_HORs")) dir.create("./tyba_island_split_HORs")
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); quit(save = "no")}
  repeats = read.csv(file = repeat_file)
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  tyba_name <- tyba$new_class[1]
  
  assembly_name_short = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  assembly_name_full = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]
  
  chr_no_sizes <- chr_no_sizes[chr_no_sizes$assembly.name == assembly_name_full, ]
  chromosomes <- chr_no_sizes$chromosome.name
  ### Load data ##################################################################
  hor_scirpt = "/home/pwlodzimierz/Rhynchospora/TRASH_v2_HORs/islands_hors_submit_slurm.sh" 
  

  if(file.exists("./hor_submission.sh")) {
    file.remove("./hor_submission.sh")
  }
  
  for(j in seq_along(chromosomes)) {
    
    tyba_chr <- tyba[tyba$seqID == chromosomes[j], ]
    islands_chr <- islands_genome[islands_genome$chromosome == chromosomes[j], ]
    
    for(k in seq_len(nrow((islands_chr)))) {
      cat(assembly_name_short, j, "/", length(chromosomes), k, "/", nrow(islands_chr), "\n")
      
      tyba_island <- tyba_chr[tyba_chr$start %in% (islands_chr$start[k] : islands_chr$end[k]), ]
      if(nrow(tyba_island) < 10) next
      
      tyba_island$seqID <- unlist(lapply(tyba_island$seqID, function(X) paste0(X, "_island", k)))
      
      output_dir <- paste0(getwd(), "/tyba_island_split_HORs")
      
      ### save islands repeat file
      repeats_file <- paste0(getwd(), "/tyba_island_split/", assembly_name_short, "_", chromosomes[j], "_array_", k, ".csv")
      write.csv(x = tyba_island, file = repeats_file, row.names = FALSE)
      
      ### make hort submission script

      cat(hor_scirpt, output_dir, repeats_file, tyba_island$seqID[1], tyba_name, repeats_file, tyba_island$seqID[1], tyba_name, assembly_name_full, assembly_name_full, "\n",
          file = "./hor_submission.sh", append = TRUE)
      cat("sleep 0.01", "\n",
          file = "./hor_submission.sh", append = TRUE)

      # cat(hor_scirpt, output_dir, repeats_file, tyba_island$seqID[1], tyba_name, repeats_file, tyba_island$seqID[1], tyba_name, assembly_name_full, assembly_name_full, "\n",
      #     file = "../hor_submission.sh", append = TRUE)
      # cat("sleep 0.01", "\n",
      #     file = "../hor_submission.sh", append = TRUE)

      
    }
    
  }
}






































































