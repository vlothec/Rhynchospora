


input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 3

keep_chrs_austr <- c(1,4,7)
keep_chrs_austr2 <- c(1,4,7)


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

suppressMessages(library(seqinr))

runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv")

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

chr_ex1 <- ""

flanking_range = 7500
flanging_skip = 2500

chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
chr_1_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_ex1]

chr_1_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_ex1]


query_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_data/query_", flanking_range, "_", flanging_skip, "_", genome_ex1, "_", chr_ex1, ".fa")
subject_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_data/subject_", genome_ex1, "_", chr_ex1, ".fa")

arrays_1 <- read.csv(runs_todo$arrays[input_id])



print("reading assembly")
assemblies <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", recursive = T, full.names = T)
assemblies <- assemblies[grep(genome_ex1, assemblies)[1]]
assembly <- read.fasta(file = assemblies)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
chr_1_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_ex1]
chr_1_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_ex1]

if(genome_ex1 == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
  keep_chrs <- keep_chrs_austr
  chr_1_sizes <- chr_1_sizes[keep_chrs]
  chr_1_names <- chr_1_names[keep_chrs]
  assembly <- assembly[keep_chrs]
  arrays_1 <- arrays_1[arrays_1$chromosome %in% names(assembly),]
}

if(genome_ex1 == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
  keep_chrs <- keep_chrs_austr2
  chr_1_sizes <- chr_1_sizes[keep_chrs]
  chr_1_names <- chr_1_names[keep_chrs]
  assembly <- assembly[keep_chrs]
  arrays_1 <- arrays_1[arrays_1$chromosome %in% names(assembly),]
}


assembly_merged <- NULL
for(i in 1 : length(assembly)) {
  assembly_merged <- c(assembly_merged, assembly[[i]])
}


chr_1_starts_adjust <- 0
for(i in 1 : (length(chr_1_sizes) - 1)) {
  chr_1_starts_adjust <- c(chr_1_starts_adjust, sum(chr_1_sizes[1 : (i)]))
}


for(i in 1 : length(chr_1_sizes)) {
  arrays_1$start[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$start[arrays_1$chromosome == chr_1_names[i]] + chr_1_starts_adjust[i]
  arrays_1$end[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$end[arrays_1$chromosome == chr_1_names[i]] + chr_1_starts_adjust[i]
}





print("saving query sequences")

for(i in 1 : nrow(arrays_1)) {
  if((arrays_1$start[i] - flanking_range - flanging_skip) < 0) {
    write.fasta(sequences = assembly_merged[(arrays_1$end[i] + flanging_skip) : (arrays_1$end[i] + flanging_skip + flanking_range)],
                names = paste0("down_", i), file.out = query_file, open = "a")
  } else if(arrays_1$end[i] + flanging_skip + flanking_range > length(assembly_merged)) {
    write.fasta(sequences = assembly_merged[(arrays_1$start[i] - flanking_range - flanging_skip) : (arrays_1$start[i] - flanging_skip)],
                names = paste0("up_", i), file.out = query_file, open = "a")
  } else {
    write.fasta(sequences = assembly_merged[(arrays_1$start[i] - flanking_range - flanging_skip) : (arrays_1$start[i] - flanging_skip)],
                names = paste0("up_", i), file.out = query_file, open = "a")
    write.fasta(sequences = assembly_merged[(arrays_1$end[i] + flanging_skip) : (arrays_1$end[i] + flanging_skip + flanking_range)],
                names = paste0("down_", i), file.out = query_file, open = "a")
  }
  
  
}
print("saving subject sequences")

write.fasta(sequences = assembly_merged, names = paste0("merged_", genome_ex1), file.out = subject_file, open = "w")





















