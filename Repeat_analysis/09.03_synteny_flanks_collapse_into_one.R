
input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 671

###
print("set up data")
###
flanking_range = 7500
flanging_skip = 2500


max_unpaired_distance = flanging_skip * 2 + 500000


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
chr_1_size <- chr_sizes$size[chr_sizes$assembly.name == genome_A & chr_sizes$chromosome.name == chromosome_A][1]
chr_2_size <- chr_sizes$size[chr_sizes$assembly.name == genome_B & chr_sizes$chromosome.name == chromosome_B][1]

arrays_1 <- read.csv(arrays_file_A)
arrays_1 <- arrays_1[arrays_1$chromosome == chromosome_A, ]

arrays_2 <- read.csv(arrays_file_B)
arrays_2 <- arrays_2[arrays_2$chromosome == chromosome_B, ]

arrays_1$array_middle <- arrays_1$start + (arrays_1$end - arrays_1$start)/2
arrays_2$array_middle <- arrays_2$start + (arrays_2$end - arrays_2$start)/2

match_data_A <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                                       flanging_skip, "_", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
match_data_B <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_out_data/", flanking_range, "_", 
                                       flanging_skip, "_", genome_B, "_", chromosome_B, "_", genome_A, "_", chromosome_A, ".csv"))



###
print("Set up merged data frame")
###


match_data_A <- match_data_A[,-1]
match_data_B <- match_data_B[,-1]
match_data_A <- unique(match_data_A)
match_data_B <- unique(match_data_B)

# filter short alignment lengths (less than 25% of query length)
match_data_A <- match_data_A[match_data_A$alignment_length > flanking_range*0.25,]
match_data_B <- match_data_B[match_data_B$alignment_length > flanking_range*0.25,]


synteny_pairs_df <- data.frame(genome_A = vector(mode = "character"),
                               chromosome_A = vector(mode = "character"),
                               is_island_A = vector(mode = "logical"),
                               island_A_ID = vector(mode = "numeric"),
                               Tyba_bp_if_island_A = vector(mode = "numeric"),
                               flanked_reg_start_A = vector(mode = "numeric"),
                               flanked_reg_end_A = vector(mode = "numeric"),
                               genome_B = vector(mode = "character"),
                               chromosome_B = vector(mode = "character"),
                               is_island_B = vector(mode = "logical"),
                               island_B_ID = vector(mode = "numeric"),
                               Tyba_bp_if_island_B = vector(mode = "numeric"),
                               flanked_reg_start_B = vector(mode = "numeric"),
                               flanked_reg_end_B = vector(mode = "numeric"),
                               reversed = vector(mode = "logical"))



for(i in 1 : nrow(arrays_1)) {
  hits_AB <- match_data_A[match_data_A$index == i, ]
  
  hits_down <- hits_AB[hits_AB$type == "down", ]
  hits_up <- hits_AB[hits_AB$type == "up", ]
  
  potential_hits <- data.frame(hit_up_id = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_id = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_start = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_end = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_direction = vector(mode = "character", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_alilength = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_start = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_end = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_direction = vector(mode = "character", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_alilength = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)))
  
  if(nrow(potential_hits) == 0) next
  
  df_id <- 0
  for(j in 1 : nrow(hits_up)) {
    hitup <- hits_up[j,]
    for(k in 1 : nrow(hits_down)) {
      hitdown <- hits_down[k,]
      df_id <- df_id + 1
      potential_hits$hit_up_id[df_id] = j
      potential_hits$hit_down_id[df_id] = k
      potential_hits$hit_up_start[df_id] = hitup$subject_start
      potential_hits$hit_up_end[df_id] = hitup$subject_end
      potential_hits$hit_up_direction[df_id] = hitup$strand
      potential_hits$hit_up_alilength[df_id] = hitup$alignment_length
      potential_hits$hit_down_start[df_id] = hitdown$subject_start
      potential_hits$hit_down_end[df_id] = hitdown$subject_start
      potential_hits$hit_down_direction[df_id] = hitdown$strand
      potential_hits$hit_down_alilength[df_id] = hitdown$alignment_length
    }
  }
  
  # check the Tyba between and potential second array and the distance
  potential_hits$Tyba_between <- 0
  potential_hits$arrays_between <- ""
  potential_hits$distance_between <- ""
  wrong_directions <- NULL
  for(j in 1 : nrow(potential_hits)) {
    # ignore and later remove those that are the wrong direcition
    if(potential_hits$hit_up_direction[j] != potential_hits$hit_down_direction[j]) {
      wrong_directions <- c(wrong_directions, j)
      next
    }
    
    if(potential_hits$hit_up_direction[j] == "+" & potential_hits$hit_down_direction[j] == "+") {
      if(potential_hits$hit_up_start[j] > potential_hits$hit_down_end[j]) {
        wrong_directions <- c(wrong_directions, j)
        next
      }
    }
    
    if(potential_hits$hit_up_direction[j] == "-" & potential_hits$hit_down_direction[j] == "-") {
      if(potential_hits$hit_up_end[j] < potential_hits$hit_down_start[j]) {
        wrong_directions <- c(wrong_directions, j)
        next
      }
    }
    which_arr2 <- NA
    if(potential_hits$hit_up_direction[j] == "+") {
      which_arr2 <- which(arrays_2$array_middle > potential_hits$hit_up_end[j] & arrays_2$array_middle < potential_hits$hit_down_start[j])
      potential_hits$distance_between[j] <- potential_hits$hit_down_start[j] - potential_hits$hit_up_end[j]
    } else {
      which_arr2 <- which(arrays_2$array_middle > potential_hits$hit_down_end[j] & arrays_2$array_middle <  potential_hits$hit_up_start[j])
      potential_hits$distance_between[j] <- potential_hits$hit_up_start[j] - potential_hits$hit_down_end[j]
    }
    if(length(which_arr2) == 0) next
    
    potential_hits$Tyba_between[j] <- sum(arrays_2$tyba_total_bp[which_arr2])
    arr_size_diff <- abs(arrays_2$tyba_total_bp[which_arr2] - arrays_1$tyba_total_bp[j])
    
    potential_hits$arrays_between[j] <- which_arr2[which.min(arr_size_diff)]
    
  }
  if(length(wrong_directions) != 0) potential_hits <- potential_hits[-wrong_directions, ]
  
  if(nrow(potential_hits) == 0) next
  
  potential_hits$alignment_total_length <- potential_hits$hit_up_alilength + potential_hits$hit_down_alilength
  
  # filter out hits that would surround a region too long to reasonably be still the same syntenic one
  potential_hits <- potential_hits[potential_hits$distance_between < max_unpaired_distance, ]
  if(nrow(potential_hits) == 0) next
  
  # choose the best one 
  
  if(sum(potential_hits$Tyba_between > 0) > 0) {
    # more than zero with array within, choose one with the longest alignment length
    consider <- which(potential_hits$Tyba_between > 0)
    best_hit <- consider[which.max(potential_hits$alignment_total_length[consider])]
  } else {
    # choose one with the longest alignment length, doesn't have the array within
    best_hit <- which.max(potential_hits$alignment_total_length)
  }
  potential_hits$arrays_between <- as.numeric(potential_hits$arrays_between)
  
  if(!is.na(potential_hits$arrays_between[best_hit])) {
    flanked_reg_start_B <- arrays_2$start[potential_hits$arrays_between[best_hit]] - flanging_skip
    flanked_reg_end_B <- arrays_2$end[potential_hits$arrays_between[best_hit]] + flanging_skip
  } else {
    if(potential_hits$hit_up_direction[best_hit] == "+") {
      flanked_reg_start_B <- potential_hits$hit_up_end[best_hit]
      flanked_reg_end_B <- potential_hits$hit_down_start[best_hit]
    } else {
      flanked_reg_start_B <- potential_hits$hit_down_end[best_hit]
      flanked_reg_end_B <- potential_hits$hit_up_start[best_hit]
    }
    
  }
  
  
  
  synteny_pairs_df <- rbind(synteny_pairs_df, data.frame(genome_A = genome_A,
                                                         chromosome_A = chromosome_A,
                                                         is_island_A = T,
                                                         island_A_ID = i,
                                                         Tyba_bp_if_island_A = arrays_1$tyba_total_bp[i],
                                                         flanked_reg_start_A = arrays_1$start[i] - flanging_skip,
                                                         flanked_reg_end_A = arrays_1$end[i] + flanging_skip,
                                                         genome_B = genome_B,
                                                         chromosome_B = chromosome_B,
                                                         is_island_B = ifelse(test = !is.na(potential_hits$arrays_between[best_hit]), yes = T, no = F),
                                                         island_B_ID = potential_hits$arrays_between[best_hit],
                                                         Tyba_bp_if_island_B = ifelse(test = potential_hits$arrays_between[best_hit] != "", yes = arrays_2$tyba_total_bp[potential_hits$arrays_between[best_hit]], no = NA),
                                                         flanked_reg_start_B = flanked_reg_start_B,
                                                         flanked_reg_end_B = flanked_reg_end_B,
                                                         reversed = ifelse(test = (potential_hits$hit_up_direction[best_hit] == "+"), yes = F, no = T)))
  
}

# the other direction




for(i in 1 : nrow(arrays_2)) {
  hits_AB <- match_data_B[match_data_B$index == i, ]
  
  hits_down <- hits_AB[hits_AB$type == "down", ]
  hits_up <- hits_AB[hits_AB$type == "up", ]
  
  potential_hits <- data.frame(hit_up_id = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_id = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_start = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_end = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_direction = vector(mode = "character", length = nrow(hits_down) * nrow(hits_up)),
                               hit_up_alilength = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_start = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_end = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_direction = vector(mode = "character", length = nrow(hits_down) * nrow(hits_up)),
                               hit_down_alilength = vector(mode = "numeric", length = nrow(hits_down) * nrow(hits_up)))
  
  if(nrow(potential_hits) == 0) next
  
  df_id <- 0
  for(j in 1 : nrow(hits_up)) {
    hitup <- hits_up[j,]
    for(k in 1 : nrow(hits_down)) {
      hitdown <- hits_down[k,]
      df_id <- df_id + 1
      potential_hits$hit_up_id[df_id] = j
      potential_hits$hit_down_id[df_id] = k
      potential_hits$hit_up_start[df_id] = hitup$subject_start
      potential_hits$hit_up_end[df_id] = hitup$subject_end
      potential_hits$hit_up_direction[df_id] = hitup$strand
      potential_hits$hit_up_alilength[df_id] = hitup$alignment_length
      potential_hits$hit_down_start[df_id] = hitdown$subject_start
      potential_hits$hit_down_end[df_id] = hitdown$subject_start
      potential_hits$hit_down_direction[df_id] = hitdown$strand
      potential_hits$hit_down_alilength[df_id] = hitdown$alignment_length
    }
  }
  
  # check the Tyba between and potential second array and the distance
  potential_hits$Tyba_between <- 0
  potential_hits$arrays_between <- ""
  potential_hits$distance_between <- 0
  wrong_directions <- NULL
  for(j in 1 : nrow(potential_hits)) {
    # ignore and later remove those that are the wrong direcition
    if(potential_hits$hit_up_direction[j] != potential_hits$hit_down_direction[j]) {
      wrong_directions <- c(wrong_directions, j)
      next
    }
    
    if(potential_hits$hit_up_direction[j] == "+" & potential_hits$hit_down_direction[j] == "+") {
      if(potential_hits$hit_up_start[j] > potential_hits$hit_down_end[j]) {
        wrong_directions <- c(wrong_directions, j)
        next
      }
    }
    
    if(potential_hits$hit_up_direction[j] == "-" & potential_hits$hit_down_direction[j] == "-") {
      if(potential_hits$hit_up_end[j] < potential_hits$hit_down_start[j]) {
        wrong_directions <- c(wrong_directions, j)
        next
      }
    }
    which_arr2 <- NA
    if(potential_hits$hit_up_direction[j] == "+") {
      which_arr2 <- which(arrays_1$array_middle > potential_hits$hit_up_end[j] & arrays_1$array_middle < potential_hits$hit_down_start[j])
      potential_hits$distance_between[j] <- potential_hits$hit_down_start[j] - potential_hits$hit_up_end[j]
    } else {
      which_arr2 <- which(arrays_1$array_middle > potential_hits$hit_down_end[j] & arrays_1$array_middle <  potential_hits$hit_up_start[j])
      potential_hits$distance_between[j] <- potential_hits$hit_up_start[j] - potential_hits$hit_down_end[j]
    }
    if(length(which_arr2) == 0) next
    
    potential_hits$Tyba_between[j] <- sum(arrays_1$tyba_total_bp[which_arr2])
    arr_size_diff <- abs(arrays_1$tyba_total_bp[which_arr2] - arrays_2$tyba_total_bp[j])
    
    potential_hits$arrays_between[j] <- which_arr2[which.min(arr_size_diff)]
    
  }
  if(length(wrong_directions) != 0) potential_hits <- potential_hits[-wrong_directions, ]
  
  if(nrow(potential_hits) == 0) next
  
  potential_hits$alignment_total_length <- potential_hits$hit_up_alilength + potential_hits$hit_down_alilength
  
  # filter out hits that would surround a region too long to reasonably be still the same syntenic one
  potential_hits <- potential_hits[(potential_hits$distance_between - potential_hits$Tyba_between)  < max_unpaired_distance, ]
  if(nrow(potential_hits) == 0) next
  
  # choose the best one 
  
  if(sum(potential_hits$Tyba_between > 0) > 0) {
    # more than zero with array within, choose one with the longest alignment length
    consider <- which(potential_hits$Tyba_between > 0)
    best_hit <- consider[which.max(potential_hits$alignment_total_length[consider])]
  } else {
    # choose one with the longest alignment length, doesn't have the array within
    best_hit <- which.max(potential_hits$alignment_total_length)
  }
  potential_hits$arrays_between <- as.numeric(potential_hits$arrays_between)
  
  if(!is.na(potential_hits$arrays_between[best_hit])) {
    flanked_reg_start_B <- arrays_1$start[potential_hits$arrays_between[best_hit]] - flanging_skip
    flanked_reg_end_B <- arrays_1$end[potential_hits$arrays_between[best_hit]] + flanging_skip
  } else {
    if(potential_hits$hit_up_direction[best_hit] == "+") {
      flanked_reg_start_B <- potential_hits$hit_up_end[best_hit]
      flanked_reg_end_B <- potential_hits$hit_down_start[best_hit]
    } else {
      flanked_reg_start_B <- potential_hits$hit_down_end[best_hit]
      flanked_reg_end_B <- potential_hits$hit_up_start[best_hit]
    }
    
  }
  
  # check if already exists in data frame, if exists but doesnt match the array, issue a warning and double up the entry
  if(i %in% synteny_pairs_df$island_B_ID) {
    which_entry <- which(synteny_pairs_df$island_B_ID == i)
    # it exists, check if it found something and it's the same array
    if(!is.na(potential_hits$arrays_between[best_hit])) {
      if(potential_hits$arrays_between[best_hit] %in% synteny_pairs_df$island_A_ID[which_entry]) {
        # it exists and matches at least the correct array, no need to do anything and append the data
        next
      } else {
        # it exists but matches a different array, issue warning and add the new entry
        message(paste0(genome_B, "array", i, "found in", genome_A, "but matching a different array than the reciprocal search"))
      }
    } else {
      # it exists in the first df, but this did not found any match
      message(paste0(genome_B, "array", i, "found in", genome_A, "but matching a non-array region while reciprocal search found this one"))
    }
  } else {
    # it does not exist in the first df, even though this search found something, add it as normal with no warning
  }
  
  
  
  
  
  synteny_pairs_df <- rbind(synteny_pairs_df, data.frame(genome_A = genome_A,
                                                         chromosome_A = chromosome_A,
                                                         is_island_A = ifelse(test = !is.na(potential_hits$arrays_between[best_hit]), yes = T, no = F),
                                                         island_A_ID = potential_hits$arrays_between[best_hit],
                                                         Tyba_bp_if_island_A = ifelse(test = potential_hits$arrays_between[best_hit] != "", yes = arrays_1$tyba_total_bp[potential_hits$arrays_between[best_hit]], no = NA),
                                                         flanked_reg_start_A = flanked_reg_start_B,
                                                         flanked_reg_end_A = flanked_reg_end_B,
                                                         genome_B = genome_B,
                                                         chromosome_B = chromosome_B,
                                                         is_island_B = T,
                                                         island_B_ID = i,
                                                         Tyba_bp_if_island_B = arrays_2$tyba_total_bp[i],
                                                         flanked_reg_start_B = arrays_2$start[i] - flanging_skip,
                                                         flanked_reg_end_B = arrays_2$end[i] + flanging_skip,
                                                         reversed = ifelse(test = potential_hits$hit_up_direction[best_hit] == "+", yes = F, no = T)))
  
}

synteny_pairs_df <- synteny_pairs_df[order(synteny_pairs_df$flanked_reg_start_A, decreasing = F), ]



write.csv(x = synteny_pairs_df, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.031_data/synteny_pairs_", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"), row.names = F)
