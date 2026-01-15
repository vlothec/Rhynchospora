input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 2

###
print("set up data")
###

# keep_chrs_austr <- c(1,4,7)
# keep_chrs_austr2 <- c(1,4,7)
remove_multiple_connections <- TRUE

{
  
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
  
  
  gaps_between_chrs <- 10000000
  
  
  runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table.csv")
  
  repeats_file_A <- runs_todo$repeats_file_A[input_id]
  repeats_file_B <- runs_todo$repeats_file_B[input_id]
  
  arrays_file_A <- runs_todo$arrays_file_A[input_id]
  arrays_file_B <- runs_todo$arrays_file_B[input_id]
  
  genome_A <- runs_todo$genome_A[input_id]
  genome_B <- runs_todo$genome_B[input_id]

  chromosome_A <- runs_todo$chromosome_A[input_id]
  chromosome_B <- runs_todo$chromosome_B[input_id]
  
  reverseA = F
  if(genome_A == "Rhync_tenuis_ref.hap1.chr.fasta" & chromosome_A == "Chr2_h1") {
    reverseA = T
  }
  if(genome_A == "Rhync_tenuis_6344A.JGV.hap1.chr.fasta" & chromosome_A == "JGV_Chr1_h2_h1") {
    reverseA = T
  }
  if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta" & chromosome_A == "scaffold_16_b2") {
    reverseA = T
  }
  reverseB = F
  if(genome_B == "Rhync_tenuis_ref.hap1.chr.fasta" & chromosome_B == "Chr2_h1") {
    reverseB = T
  }
  if(genome_B == "Rhync_tenuis_6344A.JGV.hap1.chr.fasta" & chromosome_B == "JGV_Chr1_h2_h1") {
    reverseB = T
  }
  if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta" & chromosome_B == "scaffold_16_b2") {
    reverseB = T
  }
  
  
  arrays_1 <- read.csv(arrays_file_A)
  arrays_1 <- arrays_1[arrays_1$chromosome == chromosome_A, ]

  arrays_2 <- read.csv(arrays_file_B)
  arrays_2 <- arrays_2[arrays_2$chromosome == chromosome_B, ]
  
  
  arrays_1$array_middle <- arrays_1$start + (arrays_1$end - arrays_1$start)/2
  arrays_2$array_middle <- arrays_2$start + (arrays_2$end - arrays_2$start)/2 
  
  synteny_pairs_df <- read.csv(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.031_data/synteny_pairs_", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, ".csv"))
  
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
  chr_1_size <- chr_sizes$size[chr_sizes$assembly.name == genome_A & chr_sizes$chromosome.name == chromosome_A][1]
  chr_2_size <- chr_sizes$size[chr_sizes$assembly.name == genome_B & chr_sizes$chromosome.name == chromosome_B][1]
  
  if(reverseA) {
    arrays_1$start <- chr_1_size - arrays_1$start
    arrays_1$end <- chr_1_size - arrays_1$end
  }
  if(reverseB) {
    arrays_2$start <- chr_2_size - arrays_2$start
    arrays_2$end <- chr_2_size - arrays_2$end
  }
  
  
  pair_A <- which(unlist(lapply(species, function(x) {grepl(x, genome_A)})))
  pair_B <- which(unlist(lapply(species, function(x) {grepl(x, genome_B)})))
  
  
  
  synteny_pairs_df$middle_B_for_plotting <- synteny_pairs_df$flanked_reg_start_B + 
    ((synteny_pairs_df$flanked_reg_end_B - synteny_pairs_df$flanked_reg_start_B) / 2)
  synteny_pairs_df$middle_A_for_plotting <- synteny_pairs_df$flanked_reg_start_A + 
    ((synteny_pairs_df$flanked_reg_end_A - synteny_pairs_df$flanked_reg_start_A) / 2)

  synteny_pairs_df <- synteny_pairs_df[!is.na(synteny_pairs_df$middle_A_for_plotting), ]
  synteny_pairs_df <- synteny_pairs_df[!is.na(synteny_pairs_df$middle_B_for_plotting), ]
  
  if(reverseA) {
    synteny_pairs_df$middle_A_for_plotting <- chr_1_size - synteny_pairs_df$middle_A_for_plotting
    synteny_pairs_df$reversed <- !synteny_pairs_df$reversed
  }
  
  if(reverseB) {
    synteny_pairs_df$middle_B_for_plotting <- chr_2_size - synteny_pairs_df$middle_B_for_plotting
    synteny_pairs_df$reversed <- !synteny_pairs_df$reversed
  }
  
  
  ### shift the smaller chromosome values to the middle
  chr_2_starts_plot = 1
  chr_2_ends_plot = chr_2_size
  chr_1_starts_plot = 1
  chr_1_ends_plot = chr_1_size

  is_chr_1_bigger <- ifelse(chr_1_size > chr_2_size, TRUE, FALSE)
  chr_1_middle <- chr_1_size/2
  chr_2_middle <- chr_2_size/2
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
  plot_pdf = F
  
  if(genome_A == "R_pubera_ref_2n10_v2.chr.fasta" | genome_B == "R_pubera_ref_2n10_v2.chr.fasta" ) {
    if(plot_pdf) {
      pdf(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.07_chromosome_plots_filtered_synteny/", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, "_2.pdf"),
          width = round(max(c(chr_1_size, chr_2_size))/50000)/200, height = 3600/200)
      line_dist_to_arr = 0.010
      lwd_connecting <- 4
      lwd_array <- 36
      lwd_chr <- 25
    } else {
      png(filename = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.07_chromosome_plots_filtered_synteny/", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, "_2.png"), 
          width = round(max(c(chr_1_size, chr_2_size))/50000), height = 3600)#
      line_dist_to_arr = 0.010
      lwd_connecting <- 12
      lwd_array <- 90
      lwd_chr <- 64
    }
    
  } else {
    if(plot_pdf) {
      pdf(paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.07_chromosome_plots_filtered_synteny/", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, "_2.pdf"),
          width = round(max(c(chr_1_size, chr_2_size))/50000)/200, height = 1200/200)
      line_dist_to_arr = 0.028
      lwd_connecting <- 4
      lwd_array <- 33
      lwd_chr <- 25
    } else {
      png(filename = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.07_chromosome_plots_filtered_synteny/", genome_A, "_", chromosome_A, "_", genome_B, "_", chromosome_B, "_2.png"), 
          width = round(max(c(chr_1_size, chr_2_size))/50000), height = 1200)#
      line_dist_to_arr = 0.028
      lwd_connecting <- 12
      lwd_array <- 90
      lwd_chr <- 64
    }
    
  }
  
  
  
  par(mar = c(0,0,0,0), oma = c(0,0,0,0))
  
  plot(x = NULL, y = NULL, xlim = c(0, max(c(chr_1_size, chr_2_size))), xlab = "", ylab = "", ylim = c(-2,-1),
       main = "", 
       xaxs = "i",
       yaxs = "i", axes = 0)
  
  # connecting lines
  
  for(i in 1 : nrow(synteny_pairs_df)) {
    line_type <- 1
    
    if(T) { # colouring version 2, full connections coloured
      if(synteny_pairs_df$is_island_A[i] & synteny_pairs_df$is_island_B[i]) {
        # line_colour <- "#FFAC1C95"
        line_colourA <- colors_vector_transp[pair_A]
        line_colourB <- colors_vector_transp[pair_B]
      } else {
        line_colourA <- "#00000050"
        line_colourB <- "#00000050"
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
  
  # main chromosome lines
  for(i in seq_along(chr_1_starts_plot)) lines(x = c(chr_1_starts_plot[i],chr_1_ends_plot[i]), y = c(-1,-1), lend = 1, lwd = lwd_chr, col = colors_vector[pair_A])
  for(i in seq_along(chr_2_starts_plot)) lines(x = c(chr_2_starts_plot[i],chr_2_ends_plot[i]), y = c(-2,-2), lend = 1, lwd = lwd_chr, col = colors_vector[pair_B])
  
  text(x = chr_1_starts_plot, y = -1, labels = paste0(genome_A, "_", chromosome_A), cex = 2, adj = c(0,1))
  text(x = chr_2_starts_plot, y = -2, labels = paste0(genome_B, "_", chromosome_B), cex = 2, adj = c(0,0))
  
  dev.off()
  
  
}

