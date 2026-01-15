input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 15

###
print("set up data")
###

sample_array_pairs <- 250
sample_repeats_per_array <- 50

# keep_chrs_austr <- c(1,4,7)
# keep_chrs_austr2 <- c(1,4,7)


colors_vector <- c("#E69F00", 
                   "#D55E00", 
                   "#D95F02", 
                   "#FC8D62", 
                   "#E5C494", 
                   "#F0E442", 
                   "#FFD92F", 
                   "#A6D854", 
                   "#66C2A5", 
                   "#009E73", 
                   "#1B9E77", 
                   "#66A61E", 
                   "#CC79A7", 
                   "#E78AC3", 
                   "#E7298A", 
                   "#B3B3CC", 
                   "#8DA0CB", 
                   "#7570B3", 
                   "#56B4E9", 
                   "#0072B2")
colors_vector_transp <- c("#E69F0095", 
                          "#D55E0095", 
                          "#D95F0295", 
                          "#FC8D6295", 
                          "#E5C49495", 
                          "#F0E44295", 
                          "#FFD92F95", 
                          "#A6D85495", 
                          "#66C2A595", 
                          "#009E7395", 
                          "#1B9E7795", 
                          "#66A61E95", 
                          "#CC79A795", 
                          "#E78AC395", 
                          "#E7298A95", 
                          "#B3B3CC95", 
                          "#8DA0CB95", 
                          "#7570B395", 
                          "#56B4E995", 
                          "#0072B295")


runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv")

repeats_file_A <- runs_todo$repeats_file_A[input_id]
repeats_file_B <- runs_todo$repeats_file_B[input_id]

arrays_file_A <- runs_todo$arrays_file_A[input_id]
arrays_file_B <- runs_todo$arrays_file_B[input_id]

genome_A <- runs_todo$genome_A[input_id]
genome_B <- runs_todo$genome_B[input_id]


arrays_1 <- read.csv(arrays_file_A)

arrays_2 <- read.csv(arrays_file_B)

repeats_1 <- read.csv(repeats_file_A)

repeats_2 <- read.csv(repeats_file_B)

repeats_1 <- repeats_1[grepl("yba", repeats_1$new_class),]
repeats_2 <- repeats_2[grepl("yba", repeats_2$new_class),]

synteny_pairs_df <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.031_data/synteny_pairs_", genome_A, "_NA_", genome_B, "_NA.csv"))


species <- c('rugosa', 'alba', 'cephalotes', 'ridleyi', 'gaudichaudii', 
             'watsonii', 'radicans', 'pubera', 'breviuscula', 'nervosa', 
             'ciliata', 'colorata', 'tenerrima', 'filiformis', 'austrobrasiliensis', 
             'tenuis', 'riparia', 'barbata', 'holoschoenoides', 'corymbosa')






chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
chr_1_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_A]
chr_2_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_B]
chr_1_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_A]
chr_2_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_B]

chr_1_sum_sizes <- c(0,unlist(lapply(seq_along(chr_1_sizes), function(X) sum(chr_1_sizes[1:X]))))
chr_2_sum_sizes <- c(0,unlist(lapply(seq_along(chr_2_sizes), function(X) sum(chr_2_sizes[1:X]))))
chr_1_sum_starts <- chr_1_sum_sizes[1:(length(chr_1_sum_sizes) - 1)]
chr_2_sum_starts <- chr_2_sum_sizes[1:(length(chr_2_sum_sizes) - 1)]
chr_1_sum_ends <- chr_1_sum_sizes[2:(length(chr_1_sum_sizes) - 0)]
chr_2_sum_ends <- chr_2_sum_sizes[2:(length(chr_2_sum_sizes) - 0)]

for(k in 1 : (length(chr_1_sum_sizes) - 1)) {
  synteny_pairs_df$chromosome_A[synteny_pairs_df$flanked_reg_start_A > chr_1_sum_sizes[k] & synteny_pairs_df$flanked_reg_start_A < chr_1_sum_sizes[k + 1]] = chr_1_names[k]
}
for(k in 1 : (length(chr_2_sum_sizes) - 1)) {
  synteny_pairs_df$chromosome_B[synteny_pairs_df$flanked_reg_start_B > chr_2_sum_sizes[k] & synteny_pairs_df$flanked_reg_start_B < chr_2_sum_sizes[k + 1]] = chr_2_names[k]
}

# if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
#   keep_chrs <- keep_chrs_austr
#   arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
#   repeats_1 <- repeats_1[repeats_1$seqID %in% chr_1_names[keep_chrs],]
#   synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names[keep_chrs],]
#   chr_1_sizes <- chr_1_sizes[keep_chrs]
#   chr_1_names <- chr_1_names[keep_chrs]
#   chr_1_sum_sizes <- chr_1_sum_sizes[keep_chrs]
#   chr_1_sum_starts <- chr_1_sum_starts[keep_chrs]
#   chr_1_sum_ends <- chr_1_sum_ends[keep_chrs]
# } 
# if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
#   keep_chrs <- keep_chrs_austr
#   arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
#   repeats_2 <- repeats_2[repeats_2$seqID %in% chr_2_names[keep_chrs],]
#   synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names[keep_chrs],]
#   chr_2_sizes <- chr_2_sizes[keep_chrs]
#   chr_2_names <- chr_2_names[keep_chrs]
#   chr_2_sum_sizes <- chr_2_sum_sizes[keep_chrs]
#   chr_2_sum_starts <- chr_2_sum_starts[keep_chrs]
#   chr_2_sum_ends <- chr_2_sum_ends[keep_chrs]
# } 
# if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
#   keep_chrs <- keep_chrs_austr2
#   arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
#   repeats_1 <- repeats_1[repeats_1$seqID %in% chr_1_names[keep_chrs],]
#   synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names[keep_chrs],]
#   chr_1_sizes <- chr_1_sizes[keep_chrs]
#   chr_1_names <- chr_1_names[keep_chrs]
#   chr_1_sum_sizes <- chr_1_sum_sizes[keep_chrs]
#   chr_1_sum_starts <- chr_1_sum_starts[keep_chrs]
#   chr_1_sum_ends <- chr_1_sum_ends[keep_chrs]
# } 
# if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
#   keep_chrs <- keep_chrs_austr2
#   arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
#   repeats_2 <- repeats_2[repeats_2$seqID %in% chr_2_names[keep_chrs],]
#   synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names[keep_chrs],]
#   chr_2_sizes <- chr_2_sizes[keep_chrs]
#   chr_2_names <- chr_2_names[keep_chrs]
#   chr_2_sum_sizes <- chr_2_sum_sizes[keep_chrs]
#   chr_2_sum_starts <- chr_2_sum_starts[keep_chrs]
#   chr_2_sum_ends <- chr_2_sum_ends[keep_chrs]
# } 



### rearrange data
pair_A <- which(unlist(lapply(species, function(x) {grepl(x, genome_A)})))
pair_B <- which(unlist(lapply(species, function(x) {grepl(x, genome_B)})))


# adjust starts of arrays (using rearranged chr sizes vectors)

for(i in seq_along(chr_1_sizes)) {
  arrays_1$start[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$start[arrays_1$chromosome == chr_1_names[i]] + chr_1_sum_sizes[i]
  arrays_1$end[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$end[arrays_1$chromosome == chr_1_names[i]] + chr_1_sum_sizes[i]
}
for(i in seq_along(chr_2_sizes)) {
  arrays_2$start[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$start[arrays_2$chromosome == chr_2_names[i]] + chr_2_sum_sizes[i]
  arrays_2$end[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$end[arrays_2$chromosome == chr_2_names[i]] + chr_2_sum_sizes[i]
}
# adjust starts of repeats (using rearranged chr sizes vectors)

for(i in seq_along(chr_1_sizes)) {
  repeats_1$start[repeats_1$seqID == chr_1_names[i]] <- repeats_1$start[repeats_1$seqID == chr_1_names[i]] + chr_1_sum_sizes[i]
  repeats_1$end[repeats_1$seqID == chr_1_names[i]] <- repeats_1$end[repeats_1$seqID == chr_1_names[i]] + chr_1_sum_sizes[i]
}
for(i in seq_along(chr_2_sizes)) {
  repeats_2$start[repeats_2$seqID == chr_2_names[i]] <- repeats_2$start[repeats_2$seqID == chr_2_names[i]] + chr_2_sum_sizes[i]
  repeats_2$end[repeats_2$seqID == chr_2_names[i]] <- repeats_2$end[repeats_2$seqID == chr_2_names[i]] + chr_2_sum_sizes[i]
}

### order data so that better matches are at the front


synteny_pairs_df$relative_array_diff <- 100 * abs(synteny_pairs_df$Tyba_bp_if_island_A - synteny_pairs_df$Tyba_bp_if_island_B) / 
  (synteny_pairs_df$Tyba_bp_if_island_A + synteny_pairs_df$Tyba_bp_if_island_B)/2
synteny_pairs_df$relative_array_diff[is.na(synteny_pairs_df$relative_array_diff)] = max(synteny_pairs_df$relative_array_diff[!is.na(synteny_pairs_df$relative_array_diff)])

synteny_pairs_df <- synteny_pairs_df[order(synteny_pairs_df$relative_array_diff, decreasing = T), ]


# 1. array number
arr_num_A <- nrow(arrays_1)
arr_num_B <- nrow(arrays_2)

# 2. arrays found in the other one as arrays
arr_A_found_as_array_in_B <- sum(!is.na(unique(synteny_pairs_df$island_A_ID[synteny_pairs_df$is_island_B])))
arr_B_found_as_array_in_A <- sum(!is.na(unique(synteny_pairs_df$island_B_ID[synteny_pairs_df$is_island_A])))

# 3. arrays found in the other one but with missing repeats
arr_A_found_as_missing_in_B <- sum(!is.na(unique(synteny_pairs_df$island_A_ID[!synteny_pairs_df$is_island_B])))
arr_B_found_as_missing_in_A <- sum(!is.na(unique(synteny_pairs_df$island_B_ID[!synteny_pairs_df$is_island_A])))

# 4. similarity of syntenic arrays against the homologous array
sample_synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B,]
similarity_values <- NULL

if(nrow(sample_synteny_pairs_df) > 0) {
  sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
  for(i in 1 : nrow(sample_synteny_pairs_df)) {
    repeats_sample_a <- which(repeats_1$seqID == sample_synteny_pairs_df$chromosome_A[i] & repeats_1$start >= sample_synteny_pairs_df$flanked_reg_start_A[i] & repeats_1$end <= sample_synteny_pairs_df$flanked_reg_end_A[i])
    if(!length(repeats_sample_a)) next
    repeats_sample_b <- which(repeats_2$seqID == sample_synteny_pairs_df$chromosome_B[i] & repeats_2$start >= sample_synteny_pairs_df$flanked_reg_start_B[i] & repeats_2$end <= sample_synteny_pairs_df$flanked_reg_end_B[i])
    if(!length(repeats_sample_b)) next
    repeats_sample_a <- sample(repeats_sample_a, sample_repeats_per_array, replace = T)
    repeats_sample_b <- sample(repeats_sample_b, sample_repeats_per_array, replace = T)
    
    repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
    repeats_sample_b <- repeats_2$sequence[repeats_sample_b]
    
    similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
    similarity_values <- c(similarity_values, similarity)
  }
} else {
  similarity_values <- 0
}
similarity_values_A <- similarity_values

sample_synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B,]
similarity_values <- NULL

if(nrow(sample_synteny_pairs_df) > 0) {
  sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
  for(i in 1 : nrow(sample_synteny_pairs_df)) {
    repeats_sample_a <- which(repeats_1$seqID == sample_synteny_pairs_df$chromosome_A[i] & repeats_1$start >= sample_synteny_pairs_df$flanked_reg_start_A[i] & repeats_1$end <= sample_synteny_pairs_df$flanked_reg_end_A[i])
    if(!length(repeats_sample_a)) next
    repeats_sample_b <- which(repeats_2$seqID == sample_synteny_pairs_df$chromosome_B[i] & repeats_2$start >= sample_synteny_pairs_df$flanked_reg_start_B[i] & repeats_2$end <= sample_synteny_pairs_df$flanked_reg_end_B[i])
    if(!length(repeats_sample_b)) next
    repeats_sample_a <- sample(repeats_sample_a, sample_repeats_per_array, replace = T)
    repeats_sample_b <- sample(repeats_sample_b, sample_repeats_per_array, replace = T)
    
    repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
    repeats_sample_b <- repeats_2$sequence[repeats_sample_b]
    
    similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
    similarity_values <- c(similarity_values, similarity)
  }
} else {
  similarity_values <- 0
}
similarity_values_B <- similarity_values

# # 5. similarity of syntenic arrays against itself
# sample_synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B,]
# similarity_values <- NULL
# 
# if(nrow(sample_synteny_pairs_df) > 0) {
#   sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_synteny_pairs_df)) {
#     repeats_sample <- which(repeats_1$seqID == sample_synteny_pairs_df$chromosome_A[i] & repeats_1$start >= sample_synteny_pairs_df$flanked_reg_start_A[i] & repeats_1$end <= sample_synteny_pairs_df$flanked_reg_end_A[i])
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_1$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_values_A <- similarity_values
# 
# sample_synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B,]
# similarity_values <- NULL
# 
# if(nrow(sample_synteny_pairs_df) > 0) {
#   sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_synteny_pairs_df)) {
#     repeats_sample <- which(repeats_2$seqID == sample_synteny_pairs_df$chromosome_B[i] & repeats_2$start >= sample_synteny_pairs_df$flanked_reg_start_B[i] & repeats_2$end <= sample_synteny_pairs_df$flanked_reg_end_B[i])
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_2$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_2$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_values_B <- similarity_values
# 
# # 6. similarity of syntenic arrays that went missing against itself
# sample_synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & !synteny_pairs_df$is_island_B,]
# similarity_values <- NULL
# 
# if(nrow(sample_synteny_pairs_df) > 0) {
#   sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_synteny_pairs_df)) {
#     repeats_sample <- which(repeats_1$seqID == sample_synteny_pairs_df$chromosome_A[i] & repeats_1$start >= sample_synteny_pairs_df$flanked_reg_start_A[i] & repeats_1$end <= sample_synteny_pairs_df$flanked_reg_end_A[i])
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_1$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_missing_values_A <- similarity_values
# 
# sample_synteny_pairs_df <- synteny_pairs_df[!synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B,]
# similarity_values <- NULL
# 
# if(nrow(sample_synteny_pairs_df) > 0) {
#   sample_synteny_pairs_df <- sample_synteny_pairs_df[sample(1 : nrow(sample_synteny_pairs_df), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_synteny_pairs_df)) {
#     repeats_sample <- which(repeats_2$seqID == sample_synteny_pairs_df$chromosome_B[i] & repeats_2$start >= sample_synteny_pairs_df$flanked_reg_start_B[i] & repeats_2$end <= sample_synteny_pairs_df$flanked_reg_end_B[i])
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_2$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_2$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_missing_values_B <- similarity_values
# 
# # 7. similarity of non-syntenic arrays against itself
# sample_arrays <- arrays_1[!((1 : nrow(arrays_1)) %in% synteny_pairs_df$island_A_ID),]
# similarity_values <- NULL
# 
# if(nrow(sample_arrays) > 0) {
#   sample_arrays <- sample_arrays[sample(1 : nrow(sample_arrays), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_arrays)) {
#     repeats_sample <- which(repeats_1$seqID == sample_arrays$chromosome[i] & repeats_1$start >= (sample_arrays$start[i] - 2500) & repeats_1$end <= (sample_arrays$end[i] + 2500))
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_1$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_nonsyntenic_values_A <- similarity_values
# 
# 
# sample_arrays <- arrays_2[!((1 : nrow(arrays_2)) %in% synteny_pairs_df$island_B_ID),]
# similarity_values <- NULL
# 
# if(nrow(sample_arrays) > 0) {
#   sample_arrays <- sample_arrays[sample(1 : nrow(sample_arrays), sample_array_pairs, replace = T), ]
#   for(i in 1 : nrow(sample_arrays)) {
#     repeats_sample <- which(repeats_2$seqID == sample_arrays$chromosome[i] & repeats_2$start >= (sample_arrays$start[i] - 2500) & repeats_2$end <= (sample_arrays$end[i] + 2500))
#     if(!length(repeats_sample)) next
#     
#     repeats_sample_a <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     repeats_sample_b <- sample(repeats_sample, sample_repeats_per_array, replace = T)
#     
#     repeats_sample_a <- repeats_1$sequence[repeats_sample_a]
#     repeats_sample_b <- repeats_1$sequence[repeats_sample_b]
#     
#     similarity <- 100 * (1 - mean(adist(repeats_sample_a, repeats_sample_b)) / mean(nchar(c(repeats_sample_a, repeats_sample_b))))
#     similarity_values <- c(similarity_values, similarity)
#   }
# } else {
#   similarity_values <- 0
# }
# selfsimilarity_nonsyntenic_values_B <- similarity_values


arr_num_A
arr_num_B
arr_A_found_as_array_in_B
arr_B_found_as_array_in_A
arr_A_found_as_missing_in_B
arr_B_found_as_missing_in_A
similarity_values_A
similarity_values_B
# selfsimilarity_values_A
# selfsimilarity_values_B
# selfsimilarity_missing_values_A
# selfsimilarity_missing_values_B
# selfsimilarity_nonsyntenic_values_A
# selfsimilarity_nonsyntenic_values_B



# max_bar_value <- max(c(arr_num_A, arr_num_B, arr_A_found_as_array_in_B, 
#                        arr_B_found_as_array_in_A, arr_A_found_as_missing_in_B, 
#                        arr_B_found_as_missing_in_A))
# scale_factor <- 100 / max_bar_value

# # Prepare data for bar plots (scaled to -100:100)
# bar_data <- data.frame(
#   category = rep(c("Number", "Found as array", "Found as missing"), each = 2),
#   group = rep(c("A", "B"), 3),
#   scaled_value = c(arr_num_A * scale_factor, -arr_num_B * scale_factor,
#                    arr_A_found_as_array_in_B * scale_factor, -arr_B_found_as_array_in_A * scale_factor,
#                    arr_A_found_as_missing_in_B * scale_factor, -arr_B_found_as_missing_in_A * scale_factor),
#   original_value = c(arr_num_A, arr_num_B,
#                      arr_A_found_as_array_in_B, arr_B_found_as_array_in_A,
#                      arr_A_found_as_missing_in_B, arr_B_found_as_missing_in_A),
#   x_pos = rep(1:3, each = 2)
# )

# # Prepare data for boxplots
# box_data <- data.frame(
#   category = rep(c("Similarity", "Self-similarity", "Self-sim missing", "Self-sim nonsyntenic"), each = 2),
#   group = rep(c("A", "B"), 4),
#   x_pos = rep(4:7, each = 2),
#   values = I(list(
#     similarity_values_A, similarity_values_B
#     # selfsimilarity_values_A, selfsimilarity_values_B,
#     # selfsimilarity_missing_values_A, selfsimilarity_missing_values_B,
#     # selfsimilarity_nonsyntenic_values_A, selfsimilarity_nonsyntenic_values_B
#   ))
# )

save_data <- data.frame(species_A = genome_A,
                        species_B = genome_B,
                        haplotype_comparison = pair_A == pair_B,
                        arr_num_A = mean(arr_num_A),
                        arr_num_B = mean(arr_num_B),
                        arr_A_found_as_array_in_B = mean(arr_A_found_as_array_in_B),
                        arr_B_found_as_array_in_A = mean(arr_B_found_as_array_in_A),
                        arr_A_found_as_missing_in_B = mean(arr_A_found_as_missing_in_B),
                        arr_B_found_as_missing_in_A = mean(arr_B_found_as_missing_in_A),
                        similarity_values_A = mean(similarity_values_A),
                        similarity_values_B = mean(similarity_values_B)
                        # selfsimilarity_values_A = mean(selfsimilarity_values_A),
                        # selfsimilarity_values_B = mean(selfsimilarity_values_B),
                        # selfsimilarity_missing_values_A = mean(selfsimilarity_missing_values_A),
                        # selfsimilarity_missing_values_B = mean(selfsimilarity_missing_values_B),
                        # selfsimilarity_nonsyntenic_values_A = mean(selfsimilarity_nonsyntenic_values_A),
                        # selfsimilarity_nonsyntenic_values_B = mean(selfsimilarity_nonsyntenic_values_B)
                        )
write.csv(x = save_data, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.091_mirrored_boxplots/data_mirrored_plot_", genome_A, "_", genome_B, ".csv"))

# bar_cols <- c("steelblue", "steelblue", "coral","coral", "forestgreen", "forestgreen")
# box_cols <- c("purple","purple", "gold","gold", "cyan3","cyan3", "magenta", "magenta")
# 
# pdf(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.091_mirrored_boxplots/mirrored_plot_", genome_A, "_", genome_B, ".pdf"), width = 4, height = 6)
# 
# # Create the plot with dual y-axes
# par(mar = c(8, 5, 4, 5))
# plot(NULL, xlim = c(0.5, 7.5), ylim = c(-100, 100),
#      xlab = "", ylab = "", xaxt = "n", yaxt = "n",
#      main = "")
# 
# abline(h = c(10*(1:9)), lty = 6)
# abline(h = c(-10*(1:9)), lty = 6)
# 
# # Add horizontal line at y = 0
# abline(h = 0, col = "gray30", lwd = 2)
# 
# # Plot bars
# for (i in 1:nrow(bar_data)) {
#   rect(xleft = bar_data$x_pos[i] - 0.35, 
#        xright = bar_data$x_pos[i] + 0.35,
#        ybottom = 0, 
#        ytop = bar_data$scaled_value[i],
#        col = bar_cols[i],
#        border = "black")
# }
# 
# # Plot boxplots (scale to -100:100 for display)
# for (i in 1:nrow(box_data)) {
#   vals <- box_data$values[[i]]
#   if (box_data$group[i] == "B") {
#     vals <- -vals  # Flip B values to negative side
#   }
#   boxplot(vals, at = box_data$x_pos[i], add = TRUE, 
#           col = box_cols[i],
#           boxwex = 0.7, outline = FALSE, xaxt = "n", yaxt = "n")
# }
# 
# # Add x-axis labels
# axis(1, at = 1:7, 
#      labels = c("Number of arrays", "Arra-array matches", "Array-missing matches",
#                 "Array-array similarity", "Array array matched self similarity", "Array missing matched self similarity", "Array not matched self similarity"),
#      las = 2, cex.axis = 0.5)
# 
# # Add left y-axis (original scale for bars)
# left_ticks <- seq(-100, 100, by = 20)
# left_labels <- round(left_ticks / scale_factor)
# axis(2, at = left_ticks, labels = abs(left_labels), las = 1, cex.axis = 0.8)
# mtext("Count (bars)", side = 2, line = 3.5, cex = 0.9)
# 
# # Add right y-axis (percentage scale for boxplots)
# right_ticks <- seq(-100, 100, by = 20)
# right_labels <- abs(right_ticks)
# axis(4, at = right_ticks, labels = paste0(right_labels, "%"), las = 1, cex.axis = 0.8)
# mtext("Percentage (boxplots)", side = 4, line = 3.5, cex = 0.9)
# 
# # Add legend
# # legend("topright", legend = c("Group A (positive)", "Group B (negative)"),
# #        fill = c("steelblue", "coral"), border = "black", cex = 0.9)
# 
# # Add text labels for clarity
# text(x = 0.4, y = 80, labels = genome_A, col = "black", font = 2, cex = 0.5, srt = 90)
# text(x = 0.4, y = -80, labels = genome_B, col = "black", font = 2, cex = 0.5, srt = 90)
# 
# 
# 
# dev.off()



