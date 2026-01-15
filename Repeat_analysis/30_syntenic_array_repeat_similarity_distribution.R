
input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 40


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
library(seqinr)
# library(Biostrings)

split_indexes <- function(n, z) {
  base_size <- n %/% z   # minimum size of each bin
  remainder <- n %% z    # number of bins that get +1 element
  
  L <- vector("list", z)
  start <- 1
  
  for (i in seq_len(z)) {
    size <- base_size + if (i <= remainder) 1 else 0
    end <- start + size - 1
    L[[i]] <- start:end
    start <- end + 1
  }
  
  return(L)
}


###
print("set up data")
###
flanking_range = 7500
flanging_skip = 2500

chr_bins <- 40

padding <- flanging_skip * 2 + 500000 # padding is double the flanging_skip and extra bp of 1 Mb, the extra can be modified

max_unpaired_distance = flanging_skip * 2 + 50000


runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table.csv")


repeats_file_A <- runs_todo$repeats_file_A[input_id]
repeats_file_B <- runs_todo$repeats_file_B[input_id]

arrays_file_A <- runs_todo$arrays_file_A[input_id]
arrays_file_B <- runs_todo$arrays_file_B[input_id]

genome_A <- runs_todo$genome_A[input_id]
genome_B <- runs_todo$genome_B[input_id]

chromosome_A <- runs_todo$chromosome_A[input_id]
chromosome_B <- runs_todo$chromosome_B[input_id]

chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
chr_A_size <- chr_sizes$size[chr_sizes$assembly.name == genome_A & chr_sizes$chromosome.name == chromosome_A][1]
chr_B_size <- chr_sizes$size[chr_sizes$assembly.name == genome_B & chr_sizes$chromosome.name == chromosome_B][1]

arrays_A <- read.csv(arrays_file_A)
arrays_A <- arrays_A[arrays_A$chromosome == chromosome_A, ]

arrays_B <- read.csv(arrays_file_B)
arrays_B <- arrays_B[arrays_B$chromosome == chromosome_B, ]

repeats_A <- read.csv(repeats_file_A)
repeats_A <- repeats_A[repeats_A$seqID == chromosome_A, ]
repeats_A <- repeats_A[grep("Tyba", repeats_A$new_class), ]

repeats_B <- read.csv(repeats_file_B)
repeats_B <- repeats_B[repeats_B$seqID == chromosome_B, ]
repeats_B <- repeats_B[grep("Tyba", repeats_B$new_class), ]

repeats_A$seqID = unlist(lapply(repeats_A$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
repeats_B$seqID = unlist(lapply(repeats_B$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))

if(!file.exists(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                       flanging_skip, "_", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))) stop("File A does not exist")
if(!file.exists(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                       flanging_skip, "_", genome_B, "_", chromosome_B, "_", genome_A, "_", chromosome_A, ".csv"))) stop("File B does not exist")

match_data_A <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                                       flanging_skip, "_", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
match_data_B <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                                       flanging_skip, "_", genome_B, "_", chromosome_B, "_", genome_A, "_", chromosome_A, ".csv"))



unpaired_DF_A = read.csv(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.03_data/unpaired_DF_A", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
unpaired_DF_B = read.csv(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.03_data/unpaired_DF_B", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
paired_DF_A = read.csv(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.03_data/paired_DF_A", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
paired_DF_B = read.csv(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.03_data/paired_DF_B", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))



paired_DF_pairs <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.04_paired_synteny_dfs/", 
                                          genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))


min_repeats_to_count <- 10

counted_arrays <- 0

bins = 10
bins_results <- rep(0, bins)

for(j in 1 : nrow(paired_DF_pairs)) {
  
  cat(j, "/", nrow(paired_DF_pairs), paired_DF_pairs$region_a_size[j]/1000, "kb vs", paired_DF_pairs$region_b_size[j]/1000, "kb,")
  timea = Sys.time()
  
  repeats_A_array <- repeats_A[repeats_A$start >= paired_DF_pairs$array_start_A[j] & repeats_A$end <= paired_DF_pairs$array_end_A[j], ]
  if(nrow(repeats_A_array) < min_repeats_to_count) {
    print(difftime(Sys.time(), timea))
    next
  }
  
  repeats_B_array <- repeats_B[repeats_B$start >= paired_DF_pairs$array_start_B[j] & repeats_B$end <= paired_DF_pairs$array_end_B[j], ]
  if(nrow(repeats_A_array) < min_repeats_to_count) {
    print(difftime(Sys.time(), timea))
    next
  }
  
  counted_arrays <- counted_arrays + 1
  
  bins_A <- split_indexes(nrow(repeats_A_array), bins)
  bins_B <- split_indexes(nrow(repeats_B_array), bins)
  
  for(k in seq_along(bins_results)) {
    if(paired_DF_pairs$reversed[j]) {
      edit_matrix <- adist(repeats_A_array$sequence[bins_A[[k]]], repeats_B_array$sequence[bins_B[[length(bins_results) - k + 1]]])
    } else {
      edit_matrix <- adist(repeats_A_array$sequence[bins_A[[k]]], repeats_B_array$sequence[bins_B[[k]]])
    }
    edit_matrix <- edit_matrix / mean(c(repeats_A_array$width, repeats_B_array$width))
    
    if(is.na(mean(edit_matrix))) {
      edit_matrix = 0
      warning("NA found")
    }
    
    bins_results[k] <- bins_results[k] + mean(edit_matrix)
    if(is.na(bins_results[k])) warning("NA found")
  }
  
  print(difftime(Sys.time(), timea))
  
}

bins_results <- bins_results / counted_arrays

similarity_along_arrays <- 100 * (1 - bins_results)

similarity_along_arrays <- data.frame(similarity_along_arrays)
names(similarity_along_arrays) = input_id

write.csv(x = similarity_along_arrays, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/30_syntenic_array_similarity_along_array/", 
                                                     genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"), row.names = F)

### summarize when finished
if(F) {
  files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/upload_data/30_syntenic_array_similarity_along_array/", full.names = T, pattern = ".csv")
  data <- read.csv(files[1])
  for(file in files[-1]) {
    data <- cbind(data, read.csv(file))
  }
  
  pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/30_syntenic_array_similarity_along_array/30_Summary_chromosomes.pdf"))
  plot(x = 1 : nrow(data), rowSums(data)/ncol(data), type = "l", ylim = c(70,90), main = "Similarity along syntenic arrays", xlab = "Repeats similarity in bins", ylab = "Bins along both arrays")
  for(j in 1 : ncol(data)) {
    par(new = T)
    plot(x = 1 : nrow(data), data[,j], type = "l", lwd = 0.5, axes = FALSE, xlab = "", ylab = "", ylim = c(70,90),)
  }
  par(new = T)
  plot(x = 1 : nrow(data), rowSums(data)/ncol(data), type = "l", col = "orange", axes = FALSE, lwd = 2, ylim = c(70,90), xlab = "", ylab = "")
  dev.off()
  
  
  
  
}












