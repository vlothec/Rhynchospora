input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 2

###
print("set up data")
###

keep_chrs_austr <- c(1,4,7)
keep_chrs_austr2 <- c(1,4,7)

add_lineages <- TRUE
min_lineage_length <- 3

remove_multiple_connections <- TRUE

{
  rearrangements_list <- list(c(7,1,6,9,11,14,8,10,2,3,4,12,16,5,15,17,18,13), # rugosa
                              c(11,8,10,9,12,2,3,13,5,1,7,4,6), # alba  OK
                              c(6,3,9,8,4,7,2,5,1), # cepha OK
                              c(2,3,4,1,5,6), # ridl OK
                              c(2,4,1,5,3), # gaudi OK
                              c(2,5,1,3,4), # wats 
                              c(2,5,1,4,3), # radic OK
                              c(2,5,1,4,3), # puber OK
                              c(2,5,1,4,3), # breviu OK
                              c(2,5,1,3,4), # nervos OK 
                              c(2,4,1,3,5), # ciliat ok   here
                              c(2,5,1,4,3), # colorata OK
                              c(3,4,6,5,1,2,9,10,8,7), # tenerrima OK
                              c(2,5,1,4,3), # filif OK
                              c(2,3,1), # austrob OK
                              c(2,1), # tenuis    OK
                              c(2,5,1,4,3), # riparia OK
                              c(1,2,3,4,5), # barbata OK
                              c(1,2,3,4,5), # holoschoenoides OK
                              c(2,5,6,7,4,8,1,3,9)) # corymbosa OK
  
  orientation_list <- list(c(T,F,F,T,T,F,T,F,T,F,F,T,T,F,T,F,T,T), # rugosa
                           c(T,T,T,T,F,F,F,T,F,F,F,F,T,F), # alba OK
                           c(T,T,F,T,T,T,F,F,T), # cepha OK
                           c(T,T,T,F,T,T), # ridl OK?
                           c(T,T,T,F,F), # gaudi OK
                           c(F,F,F,T,T), # wats OK?
                           c(T,F,T,T,F), # radic OK?
                           c(4,5,F,2,3), # puber OK?
                           c(T,F,T,T,F), # breviu OK
                           c(F,T,F,T,F), # nervos OK 
                           c(F,T,F,T,T), # ciliat ok  here
                           c(T,F,T,T,F), # colorata OK
                           c(T,T,F,T,F,T,F,T,F,F), # tenerrima OK
                           c(T,F,T,T,F), # filif OK
                           c(T,T,T), # austrob OK
                           c(T,T), # tenuis OK
                           c(T,T,T,F,T), # riparia OK
                           c(T,F,F,F,F), # barbata OK
                           c(T,T,T,T,F), # holoschoenoides OK
                           c(F,F,T,T,T,F,T,F,T)) # corymbosa OK
  
  colors_vector <- c("#E69F00", 
                     "#D55E00", 
                     "#FC8D62", 
                     "#E5C494", 
                     "#D95F02", 
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
                            "#FC8D6295", 
                            "#E5C49495", 
                            "#D95F0295", 
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
  
  
  gaps_between_chrs <- 10000000
  
  
  runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv")
  
  
  arrays_file_A <- runs_todo$arrays_file_A[input_id]
  arrays_file_B <- runs_todo$arrays_file_B[input_id]
  
  genome_A <- runs_todo$genome_A[input_id]
  genome_B <- runs_todo$genome_B[input_id]
  
  
  arrays_1 <- read.csv(arrays_file_A)
  
  arrays_2 <- read.csv(arrays_file_B)
  
  synteny_pairs_df <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.031_data/synteny_pairs_", genome_A, "_NA_", genome_B, "_NA.csv"))
  
  if(add_lineages) {
    synteny_pairs_df$colour_A = "#00000050"
    synteny_pairs_df$colour_B = "#00000050"
    
    lineages_df <- read.csv("/home/pwlodzimierz/Rhynchospora/upload_data/9.11_data/synteny_lineages_full_genomes.csv")
    lineages_summary <- read.csv("/home/pwlodzimierz/Rhynchospora/upload_data/9.11_data/synteny_lineages_summary_full_genomes.csv")
    lineages_summary <- lineages_summary[lineages_summary$Freq >= min_lineage_length, ]
    lineages_df <- lineages_df[lineages_df$lineage_id %in% lineages_summary$Var1, ]
    
    if(nrow(lineages_df) == 0) next
    
    set.seed(123)
    lineages_summary$colour <- hsv(h = runif(nrow(lineages_summary)), 
                                   s = runif(nrow(lineages_summary), 0.5, 1), 
                                   v = runif(nrow(lineages_summary), 0.5, 1))
    
    lineages_df$colour <- ""
    
    for(k in 1 : nrow(lineages_df)) {
      lineages_df$colour[k] <- lineages_summary$colour[lineages_summary$Var1 == lineages_df$lineage_id[k]]
    }
    
    
    synteny_pairs_df$lineage_id_A <- 0
    lineages_df_A <- lineages_df[lineages_df$genome == genome_A, ]
    for(k in seq_len(nrow(lineages_df_A))) {
      synteny_pairs_df$lineage_id_A[synteny_pairs_df$island_A_ID == lineages_df_A$array_id[k]] = lineages_df_A$lineage_id[k]
      synteny_pairs_df$colour_A[synteny_pairs_df$island_A_ID == lineages_df_A$array_id[k]] = lineages_df_A$colour[k]
    }
    
    
    synteny_pairs_df$lineage_id_B <- 0
    lineages_df_B <- lineages_df[lineages_df$genome == genome_B, ]
    for(k in seq_len(nrow(lineages_df_B))) {
      synteny_pairs_df$lineage_id_B[synteny_pairs_df$island_B_ID == lineages_df_B$array_id[k]] = lineages_df_B$lineage_id[k]
      synteny_pairs_df$colour_B[synteny_pairs_df$island_B_ID == lineages_df_B$array_id[k]] = lineages_df_B$colour[k]
    }
    
    
  }
  
  if(remove_multiple_connections) {
    remove_rows <- NULL
    for(k in unique(synteny_pairs_df$island_A_ID)) {
      if(is.na(k)) next
      connections <- which(synteny_pairs_df$island_A_ID == k & !is.na(synteny_pairs_df$island_A_ID))
      if(length(connections) > 1) remove_rows <- c(remove_rows, sample(connections,(length(connections) - 1)))
      
    }
    if(length(remove_rows) != 0) synteny_pairs_df <- synteny_pairs_df[-remove_rows,]
    
    remove_rows <- NULL
    for(k in unique(synteny_pairs_df$island_B_ID)) {
      if(is.na(k)) next
      connections <- which(synteny_pairs_df$island_B_ID == k & !is.na(synteny_pairs_df$island_B_ID))
      if(length(connections) > 1) remove_rows <- c(remove_rows, sample(connections,(length(connections) - 1)))
      
    }
    if(length(remove_rows) != 0) synteny_pairs_df <- synteny_pairs_df[-remove_rows,]
  }
  
  species <- c('rugosa', 'alba', 'cephalotes', 'ridleyi', 'gaudichaudii', 
               'watsonii', 'radicans', 'pubera', 'breviuscula', 'nervosa', 
               'ciliata', 'colorata', 'tenerrima', 'filiformis', 'austrobrasiliensis', 
               'tenuis', 'riparia', 'barbata', 'holoschoenoides', 'corymbosa') 
  
  
  
  
  
  
  chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
  chr_1_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_A]
  chr_2_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_B]
  chr_1_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_A]
  chr_2_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_B]
  
  
  
  if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
    
    chr_1_sizes <- chr_1_sizes[keep_chrs]
    chr_1_names <- chr_1_names[keep_chrs]
  }
  if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
    
    chr_2_sizes <- chr_2_sizes[keep_chrs]
    chr_2_names <- chr_2_names[keep_chrs]
  }
  if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
    
    chr_1_sizes <- chr_1_sizes[keep_chrs]
    chr_1_names <- chr_1_names[keep_chrs]
  }
  if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
    
    chr_2_sizes <- chr_2_sizes[keep_chrs]
    chr_2_names <- chr_2_names[keep_chrs]
  }
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
  
  if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names,]
    synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names,]
  }
  if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names,]
    synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names,]
  }
  if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names,]
    synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names,]
  }
  if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names,]
    synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names,]
  }
  
  
  ### rearrange data
  pair_A <- which(unlist(lapply(species, function(x) {grepl(x, genome_A)})))
  pair_B <- which(unlist(lapply(species, function(x) {grepl(x, genome_B)})))
  
  chr1_rearrangement = rearrangements_list[[pair_A]]
  chr1_keep_orientation = orientation_list[[pair_A]] 
  chr2_rearrangement = rearrangements_list[[pair_B]] 
  chr2_keep_orientation = orientation_list[[pair_B]] 
  
  if(genome_B == "Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr.fasta") {
    chr2_rearrangement = c(6,4,3,5,2,1)
    chr2_keep_orientation = c(T,F,F,T,F,T)
  }
  if(genome_B == "Rhync_tenuis_ref.hap2.chr.fasta") {
    chr2_rearrangement = c(1,2)
  }
  
  chr_1_sizes_rearranged <- chr_1_sizes[chr1_rearrangement]
  chr_2_sizes_rearranged <- chr_2_sizes[chr2_rearrangement]
  chr_1_sizes_rearranged <- chr_1_sizes_rearranged[!is.na(chr_1_sizes_rearranged)]
  chr_2_sizes_rearranged <- chr_2_sizes_rearranged[!is.na(chr_2_sizes_rearranged)]
  
  arrays_1_rearr <- NULL
  for(i in chr1_rearrangement) {
    if(!(chr_1_names[i] %in% arrays_1$chromosome)) next
    arr_temp <- arrays_1[arrays_1$chromosome == chr_1_names[i], ]
    if(!chr1_keep_orientation[which(chr1_rearrangement == i)]) {
      # reverse the array coordinates
      arr_temp$start <- chr_1_sizes[i] - arr_temp$start
      arr_temp$end <- chr_1_sizes[i] - arr_temp$end
    }
    arrays_1_rearr <- rbind(arrays_1_rearr, arr_temp)
  }
  arrays_1 <- arrays_1_rearr
  
  arrays_2_rearr <- NULL
  for(i in chr2_rearrangement) {
    if(!(chr_2_names[i] %in% arrays_2$chromosome)) next
    arr_temp <- arrays_2[arrays_2$chromosome == chr_2_names[i], ]
    if(!chr2_keep_orientation[which(chr2_rearrangement == i)]) {
      # reverse the array coordinates
      arr_temp$start <- chr_2_sizes[i] - arr_temp$start
      arr_temp$end <- chr_2_sizes[i] - arr_temp$end
    }
    arrays_2_rearr <- rbind(arrays_2_rearr, arr_temp)
  }
  arrays_2 <- arrays_2_rearr
  
  
  ### rearrange synteny_pairs_df
  
  
  chr_1_sum_sizes_rearranged <- c(0,unlist(lapply(seq_along(chr_1_sizes_rearranged), function(X) sum(chr_1_sizes_rearranged[1:X]))))
  chr_2_sum_sizes_rearranged <- c(0,unlist(lapply(seq_along(chr_2_sizes_rearranged), function(X) sum(chr_2_sizes_rearranged[1:X]))))
  
  synteny_pairs_df$middle_A_for_plotting <- synteny_pairs_df$flanked_reg_start_A + (synteny_pairs_df$flanked_reg_end_A - synteny_pairs_df$flanked_reg_start_A) / 2
  synteny_pairs_df$middle_A_for_plotting_rearr = NA
  for(i in chr1_rearrangement) {
    which_values <- which(synteny_pairs_df$chromosome_A == chr_1_names[i])
    values <- synteny_pairs_df$middle_A_for_plotting[which_values]
    
    values_original <- values - chr_1_sum_starts[i]
    if(!chr1_keep_orientation[which(chr1_rearrangement == i)]) {
      values_original <- chr_1_sizes[i] - values_original
      synteny_pairs_df$reversed[which_values] <- !synteny_pairs_df$reversed[which_values]
    }
    values_replaced <- values_original + chr_1_sum_sizes_rearranged[which(chr1_rearrangement == i)]
    synteny_pairs_df$middle_A_for_plotting_rearr[which_values] <- values_replaced
  }
  
  
  synteny_pairs_df$middle_B_for_plotting <- synteny_pairs_df$flanked_reg_start_B + (synteny_pairs_df$flanked_reg_end_B - synteny_pairs_df$flanked_reg_start_B) / 2
  synteny_pairs_df$middle_B_for_plotting_rearr = NA
  for(i in chr2_rearrangement) {
    which_values <- which(synteny_pairs_df$chromosome_B == chr_2_names[i])
    values <- synteny_pairs_df$middle_B_for_plotting[which_values]
    
    values_original <- values - chr_2_sum_starts[i]
    if(!chr2_keep_orientation[which(chr2_rearrangement == i)]) {
      values_original <- chr_2_sizes[i] - values_original
      synteny_pairs_df$reversed[which_values] <- !synteny_pairs_df$reversed[which_values]
    }
    values_replaced <- values_original + chr_2_sum_sizes_rearranged[which(chr2_rearrangement == i)]
    synteny_pairs_df$middle_B_for_plotting_rearr[which_values] <- values_replaced
    
  }
  
  synteny_pairs_df$middle_B_for_plotting <- synteny_pairs_df$middle_B_for_plotting_rearr
  synteny_pairs_df$middle_A_for_plotting <- synteny_pairs_df$middle_A_for_plotting_rearr
  # if("Rhync_austrobrasiliensis_6344D.hic.hap1.chr" %in% c(genome_A, genome_B)) {
  #   synteny_pairs_df <- synteny_pairs_df[!is.na(synteny_pairs_df$middle_A_for_plotting), ]
  #   synteny_pairs_df <- synteny_pairs_df[!is.na(synteny_pairs_df$middle_B_for_plotting), ]
  # }
  
  
  # adjust starts of arrays (using rearranged chr sizes vectors)
  
  for(i in chr1_rearrangement) {
    arrays_1$start[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$start[arrays_1$chromosome == chr_1_names[i]] + chr_1_sum_sizes_rearranged[which(chr1_rearrangement == i)]
    arrays_1$end[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$end[arrays_1$chromosome == chr_1_names[i]] + chr_1_sum_sizes_rearranged[which(chr1_rearrangement == i)]
  }
  for(i in chr2_rearrangement) {
    arrays_2$start[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$start[arrays_2$chromosome == chr_2_names[i]] + chr_2_sum_sizes_rearranged[which(chr2_rearrangement == i)]
    arrays_2$end[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$end[arrays_2$chromosome == chr_2_names[i]] + chr_2_sum_sizes_rearranged[which(chr2_rearrangement == i)]
  }
  
  ### add gaps between chromosomes
  chr_1_starts_plot <- 0
  chr_1_ends_plot <- NULL
  add_gap <- 0
  for(i in chr1_rearrangement) {
    arrays_1$start[arrays_1$chromosome == chr_1_names[i]] = arrays_1$start[arrays_1$chromosome == chr_1_names[i]] + add_gap
    arrays_1$end[arrays_1$chromosome == chr_1_names[i]] = arrays_1$end[arrays_1$chromosome == chr_1_names[i]] + add_gap
    synteny_pairs_df$middle_A_for_plotting[synteny_pairs_df$chromosome_A == chr_1_names[i]] <- synteny_pairs_df$middle_A_for_plotting[synteny_pairs_df$chromosome_A == chr_1_names[i]] + add_gap
    
    chr_1_ends_plot <- c(chr_1_ends_plot, chr_1_starts_plot[which(chr1_rearrangement == i)] + chr_1_sizes_rearranged[which(chr1_rearrangement == i)])
    add_gap <- add_gap + gaps_between_chrs
    chr_1_starts_plot <- c(chr_1_starts_plot, chr_1_ends_plot[which(chr1_rearrangement == i)] + gaps_between_chrs)
  }
  chr_1_starts_plot <- chr_1_starts_plot[-length(chr_1_starts_plot)]
  
  chr_2_starts_plot <- 0
  chr_2_ends_plot <- NULL
  add_gap <- 0
  for(i in chr2_rearrangement) {
    arrays_2$start[arrays_2$chromosome == chr_2_names[i]] = arrays_2$start[arrays_2$chromosome == chr_2_names[i]] + add_gap
    arrays_2$end[arrays_2$chromosome == chr_2_names[i]] = arrays_2$end[arrays_2$chromosome == chr_2_names[i]] + add_gap
    synteny_pairs_df$middle_B_for_plotting[synteny_pairs_df$chromosome_B == chr_2_names[i]] <- synteny_pairs_df$middle_B_for_plotting[synteny_pairs_df$chromosome_B == chr_2_names[i]] + add_gap
    
    chr_2_ends_plot <- c(chr_2_ends_plot, chr_2_starts_plot[which(chr2_rearrangement == i)] + chr_2_sizes_rearranged[which(chr2_rearrangement == i)])
    add_gap <- add_gap + gaps_between_chrs
    chr_2_starts_plot <- c(chr_2_starts_plot, chr_2_ends_plot[which(chr2_rearrangement == i)] + gaps_between_chrs)
  }
  chr_2_starts_plot <- chr_2_starts_plot[-length(chr_2_starts_plot)]
  
  
  ### shift the smaller chromosome values to the middle
  
  is_chr_1_bigger <- ifelse(max(chr_1_ends_plot) > max(chr_2_ends_plot), TRUE, FALSE)
  chr_1_middle <- max(chr_1_ends_plot)/2
  chr_2_middle <- max(chr_2_ends_plot)/2
  chr_mid_diff <- round(abs(chr_1_middle - chr_2_middle))
  if(is_chr_1_bigger) {
    # adjust chr 2
    arrays_2$start <- arrays_2$start + chr_mid_diff
    arrays_2$end <- arrays_2$end + chr_mid_diff
    synteny_pairs_df$middle_B_for_plotting <- synteny_pairs_df$middle_B_for_plotting + chr_mid_diff
    chr_2_starts_plot <- chr_2_starts_plot + chr_mid_diff
    chr_2_ends_plot <- chr_2_ends_plot + chr_mid_diff
  } else {
    # adjust chr 1
    arrays_1$start <- arrays_1$start + chr_mid_diff
    arrays_1$end <- arrays_1$end + chr_mid_diff
    synteny_pairs_df$middle_A_for_plotting <- synteny_pairs_df$middle_A_for_plotting + chr_mid_diff
    chr_1_starts_plot <- chr_1_starts_plot + chr_mid_diff
    chr_1_ends_plot <- chr_1_ends_plot + chr_mid_diff
  }
  arrays_1$array_middle[arrays_1$end > arrays_1$start] <- arrays_1$start[arrays_1$end > arrays_1$start] + (arrays_1$end[arrays_1$end > arrays_1$start] - arrays_1$start[arrays_1$end > arrays_1$start])/2
  arrays_1$array_middle[arrays_1$end < arrays_1$start] <- arrays_1$end[arrays_1$end < arrays_1$start] + (arrays_1$start[arrays_1$end < arrays_1$start] - arrays_1$end[arrays_1$end < arrays_1$start])/2
  arrays_2$array_middle[arrays_2$end > arrays_2$start] <- arrays_2$start[arrays_2$end > arrays_2$start] + (arrays_2$end[arrays_2$end > arrays_2$start] - arrays_2$start[arrays_2$end > arrays_2$start])/2
  arrays_2$array_middle[arrays_2$end < arrays_2$start] <- arrays_2$end[arrays_2$end < arrays_2$start] + (arrays_2$start[arrays_2$end < arrays_2$start] - arrays_2$end[arrays_2$end < arrays_2$start])/2
  
  ### order data so that better matches are at the front
  
  
  synteny_pairs_df$relative_array_diff <- 100 * abs(synteny_pairs_df$Tyba_bp_if_island_A - synteny_pairs_df$Tyba_bp_if_island_B) / 
    (synteny_pairs_df$Tyba_bp_if_island_A + synteny_pairs_df$Tyba_bp_if_island_B)/2
  synteny_pairs_df$relative_array_diff[is.na(synteny_pairs_df$relative_array_diff)] = max(synteny_pairs_df$relative_array_diff[!is.na(synteny_pairs_df$relative_array_diff)])
  
  synteny_pairs_df <- synteny_pairs_df[order(synteny_pairs_df$relative_array_diff, decreasing = T), ]
  
  
  
  
  
  
  
  ###
  print("plot")
  ###
  chromosome_A = ""
  chromosome_B = ""
  plot_pdf = T
  
  if(genome_A == "R_pubera_ref_2n10_v2.chr.fasta" | genome_B == "R_pubera_ref_2n10_v2.chr.fasta" ) {
    if(plot_pdf) {
      pdf(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.071_genome_plots_filtered_synteny/", genome_A, "_", genome_B, "_2.pdf"),
          width = round(max(c(chr_1_ends_plot, chr_2_ends_plot))/50000)/200, height = 3600/200)
      line_dist_to_arr = 0.010
      lwd_array <- 36
      lwd_chr <- 25
    } else {
      png(filename = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.071_genome_plots_filtered_synteny/", genome_A, "_", genome_B, "_2.png"), 
          width = round(max(c(chr_1_ends_plot, chr_2_ends_plot))/50000), height = 3600)#
      line_dist_to_arr = 0.010
      lwd_array <- 90
      lwd_chr <- 64
    }
    
  } else {
    if(plot_pdf) {
      pdf(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.071_genome_plots_filtered_synteny/", genome_A, "_", genome_B, "_2.pdf"),
          width = round(max(c(chr_1_ends_plot, chr_2_ends_plot))/50000)/200, height = 1200/200)
      line_dist_to_arr = 0.028
      lwd_array <- 33
      lwd_chr <- 25
    } else {
      png(filename = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.071_genome_plots_filtered_synteny/", genome_A, "_", genome_B, "_2.png"), 
          width = round(max(c(chr_1_ends_plot, chr_2_ends_plot))/50000), height = 1200)#
      line_dist_to_arr = 0.028
      lwd_array <- 90
      lwd_chr <- 64
    }
    
  }
  
  
  
  par(mar = c(0,0,0,0), oma = c(0,0,0,0))
  
  plot(x = NULL, y = NULL, xlim = c(0, max(c(chr_1_ends_plot, chr_2_ends_plot))), xlab = "", ylab = "", ylim = c(-2,-1),
       main = "", 
       xaxs = "i",
       yaxs = "i", axes = 0)
  
  # array thicker lines
  if(T) {
    for(i in 1 : nrow(arrays_1)) {
      if(arrays_1$start[i] < arrays_1$end[i]) {
        lines(x = c(arrays_1$start[i]-100000,arrays_1$end[i]+100000), y = c(-1,-1), lwd = lwd_array, col = "#00000090", lend = 1)
      } else {
        lines(x = c(arrays_1$end[i]-100000,arrays_1$start[i]+100000), y = c(-1,-1), lwd = lwd_array, col = "#00000090", lend = 1)
      }
      
    }
    for(i in 1 : nrow(arrays_2)) {
      if(arrays_2$start[i] < arrays_2$end[i]) {
        lines(x = c(arrays_2$start[i]-100000,arrays_2$end[i]+100000), y = c(-2,-2), lwd = lwd_array, col = "#00000090", lend = 1)
      } else {
        lines(x = c(arrays_2$end[i]-100000,arrays_2$start[i]+100000), y = c(-2,-2), lwd = lwd_array, col = "#00000090", lend = 1)
      }
    }
  }
  
  # connecting lines
  
  for(i in 1 : nrow(synteny_pairs_df)) {
    line_type <- 1
    
    if(T) { # colouring version 2, full connections coloured
      if(pair_A == pair_B) {
        if(synteny_pairs_df$is_island_A[i] & synteny_pairs_df$is_island_B[i]) {
          line_colourA <- colors_vector_transp[pair_A]
          line_colourB <- colors_vector_transp[pair_B]
          lwd_connecting =  ifelse(plot_pdf, 4, 12)
        } else {
          line_colourA <- "#00000050"
          line_colourB <- "#00000050"
          lwd_connecting =  ifelse(plot_pdf, 4, 12)
        }
      } else if(add_lineages) {
        if(synteny_pairs_df$colour_A[i] != "#00000050" & synteny_pairs_df$colour_B[i] != "#00000050" & synteny_pairs_df$colour_A[i] == synteny_pairs_df$colour_B[i]) {
          line_colourA <- "#ff222295"
          line_colourB <- "#ff222295"
          lwd_connecting =  ifelse(plot_pdf, 6, 18)
        } else if(synteny_pairs_df$is_island_A[i] & synteny_pairs_df$is_island_B[i]) {
          line_colourA <- "#FFAC1C95"
          line_colourB <- "#FFAC1C95"
          lwd_connecting =  ifelse(plot_pdf, 4, 12)
        } else {
          line_colourA <- "#00000050"
          line_colourB <- "#00000050" 
          lwd_connecting =  ifelse(plot_pdf, 4, 12)
        }
      } else if(synteny_pairs_df$is_island_A[i] & synteny_pairs_df$is_island_B[i]) {
        line_colourA <- "#FFAC1C95"
        line_colourB <- "#FFAC1C95"
        lwd_connecting =  ifelse(plot_pdf, 4, 12)
      } else {
        line_colourA <- "#00000050"
        line_colourB <- "#00000050"
        lwd_connecting =  ifelse(plot_pdf, 4, 12)
      }
      
      
    }
    
    if(synteny_pairs_df$reversed[i]) line_type <- 2
    
    middle_x_value <- (synteny_pairs_df$middle_B_for_plotting[i] + synteny_pairs_df$middle_A_for_plotting[i]) / 2
    
    lines(x = c(synteny_pairs_df$middle_A_for_plotting[i], middle_x_value), 
          y = c(-1 - line_dist_to_arr ,-1.5), 
          lwd = lwd_connecting, col = line_colourA, lend = 1, lty = line_type)
    
    lines(x = c(middle_x_value, synteny_pairs_df$middle_B_for_plotting[i]), 
          y = c(-1.5,-2 + line_dist_to_arr), 
          lwd = lwd_connecting, col = line_colourB, lend = 1, lty = line_type)
  }
  
  # main chromosome lines
  for(i in seq_along(chr_1_starts_plot)) lines(x = c(chr_1_starts_plot[i],chr_1_ends_plot[i]), y = c(-1,-1), lend = 1, lwd = lwd_chr, col = colors_vector[pair_A])
  for(i in seq_along(chr_2_starts_plot)) lines(x = c(chr_2_starts_plot[i],chr_2_ends_plot[i]), y = c(-2,-2), lend = 1, lwd = lwd_chr, col = colors_vector[pair_B])
  
  
  
  dev.off()
  
  
}

