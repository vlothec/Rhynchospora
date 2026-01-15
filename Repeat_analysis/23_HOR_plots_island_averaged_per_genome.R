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



chr_no_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
print(paste0("A: ", i, " / ", length(data_directories)))
### Load data ==================================================================
setwd(data_directories[i])
print(data_directories[i])

assembly_name_short = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
assembly_name_full = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]

hor_files <- list.files(path = "./tyba_island_split_HORs/", pattern = "HORs_Tyba_", full.names = T)
repeat_files <- list.files(path = "./tyba_island_split_HORs/", pattern = "repeats_with_hors_Tyba_", full.names = T)

windows <- 25
score_per_window <- rep(0, windows)
repeats_per_window <- rep(0, windows)

cat("/n")
cat("=========================/n")
for(j in seq_along(hor_files)) {
  if(j %in% round((1 : 25) * (length(hor_files) / 25))) cat("#")
  
  hors <- read.csv(file = hor_files[j])
  island_ID <- strsplit(strsplit(hor_files[j], split = "//HORs_")[[1]][2], split = "sv")[[1]][1]
  repeat_file <- repeat_files[grep(island_ID, repeat_files)]
  if(length(repeat_file) != 1) next
  repeats <- read.csv(file = repeat_file)
  
  repeats$novel_score <- 0
  
  for(k in seq_len(nrow(hors))) {
    repeats$novel_score[hors$start_A[k] : hors$end_A[k]] = repeats$novel_score[hors$start_A[k] : hors$end_A[k]] + 1
    repeats$novel_score[hors$start_B[k] : hors$end_B[k]] = repeats$novel_score[hors$start_B[k] : hors$end_B[k]] + 1
  }
  repeats$novel_score = repeats$novel_score / nrow(repeats)
  
  repeats$novel_score <- repeats$novel_score * 100
  
  
  indices <- 1 : length(repeats$novel_score)
  bin_size <- floor(length(repeats$novel_score) / windows)
  remainder <- length(repeats$novel_score) %% windows
  bin_sizes <- rep(bin_size, windows)
  if (remainder > 0) {
    extra_bins <- sample(1 : windows, remainder, replace = FALSE)
    bin_sizes[extra_bins] <- bin_sizes[extra_bins] + 1
  }
  bin_assignments <- rep(1 : windows, times = bin_sizes)
  
  bin_assignments <- bin_assignments[1 : length(repeats$novel_score)]
  
  bins <- split(indices, bin_assignments)
  
  # score_per_window <- rep(0, length(bins))
  # repeats_per_window <- rep(0, length(bins))
  
  for(k in seq_along(bins)) {
    score_per_window[k] <- score_per_window[k] + sum(repeats$novel_score[bins[k][[1]]])
    repeats_per_window[k] <- repeats_per_window[k] + length(bins[k][[1]])
  }
}
cat("/n")

averaged_scores <- score_per_window / repeats_per_window

if(sum(unique(averaged_scores) == 0) == windows) stop("Every score is 0")

write.csv(x = data.frame(repeat_number_in_bin = repeats_per_window,
                         hor_novel_score_per_bin = score_per_window,
                         genome = rep(assembly_name_short, windows)), 
          file = paste0("./", assembly_name_short, "_HOR_island_averaged_values.csv"), 
          row.names = FALSE)

write.csv(x = data.frame(repeat_number_in_bin = repeats_per_window,
                         hor_novel_score_per_bin = score_per_window,
                         genome = rep(assembly_name_short, windows)), 
          file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/23_hor_csvs_averaged_island_per_genome/", assembly_name_short, "_HOR_island_averaged_values.csv"), 
          row.names = FALSE)

pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/23_hor_plots_averaged_island_per_genome/", assembly_name_short, "_HOR_island_averaged_values.pdf"), 
    width = 7, height = 4)
plot(x = 1 : windows, y = averaged_scores, 
     main = paste0(assembly_name_short, " HOR scores per island"), 
     xlab = paste0(windows, " windows averaged scores of ", length(hor_files), " islands, totaling ", sum(repeats_per_window), " repeats"),
     ylab = "HOR score: % of other repeats a repeat forms a HOR with", 
     cex.lab = 0.7, cex.axis = 0.7, cex.main = 1, cex.sub = 0.5, 
     type = "b")
dev.off()









