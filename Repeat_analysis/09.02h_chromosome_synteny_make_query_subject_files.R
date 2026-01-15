


input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 1 : 206


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

suppressMessages(library(seqinr))

runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table.csv")

runs_todo2 <- data.frame(genome = vector(mode = "character"),
                         chromosome = vector(mode = "character"),
                         arrays = vector(mode = "character"))
                         
for(i in 1 : nrow(runs_todo)) {
  runs_todo2 <- rbind(runs_todo2, data.frame(genome = runs_todo$genome_A[i],
                                             chromosome = runs_todo$chromosome_A[i],
                                             arrays = runs_todo$arrays_file_A[i]))
  
  runs_todo2 <- rbind(runs_todo2, data.frame(genome = runs_todo$genome_B[i],
                                             chromosome = runs_todo$chromosome_B[i],
                                             arrays = runs_todo$arrays_file_B[i]))
}
runs_todo <- unique(runs_todo2)



genome_ex1 <- runs_todo$genome[input_id]

chr_ex1 <- runs_todo$chromosome[input_id]

flanking_range = 7500
flanging_skip = 2500


query_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_data/query_", flanking_range, "_", flanging_skip, "_", genome_ex1, "_", chr_ex1, ".fa")
subject_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_data/subject_", genome_ex1, "_", chr_ex1, ".fa")

arrays_1 <- read.csv(runs_todo$arrays[input_id])
arrays_1 <- arrays_1[arrays_1$chromosome == chr_ex1, ]

print("reading query assembly")
assemblies <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", recursive = T, full.names = T)
assemblies <- assemblies[grep(genome_ex1, assemblies)[1]]
assembly <- read.fasta(file = assemblies)
assembly <- assembly[[which(names(assembly) == chr_ex1)]]

print("saving query sequences")

for(i in 1 : nrow(arrays_1)) {
  write.fasta(sequences = assembly[(arrays_1$start[i] - flanking_range - flanging_skip) : (arrays_1$start[i] - flanging_skip)],
              names = paste0("up_", i), file.out = query_file, open = "a")
  write.fasta(sequences = assembly[(arrays_1$end[i] + flanging_skip) : (arrays_1$end[i] + flanging_skip + flanking_range)],
              names = paste0("down_", i), file.out = query_file, open = "a")
  
}

if(!file.exists(subject_file)) { # create the subject chromosome file
  
  print("saving subject sequences")
  
  write.fasta(sequences = assembly, names = paste0(genome_ex1, "_", chr_ex1), file.out = subject_file, open = "w")
  
}






