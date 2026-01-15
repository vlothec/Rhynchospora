

input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 27

# install.packages("devtools") 
# library(devtools)
# install_github("evolvedmicrobe/dotplot", build_vignettes = FALSE)

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
library(dotplot) # uses conda activate R-4.1.3
library(ggplot2)
library(seqinr)
# library(Biostrings)

revCompString = function(DNAstrVector) 
{
  DNAstrVector <- toupper(DNAstrVector)
  newDNAstrVector <- DNAstrVector
  newDNAstrVector[DNAstrVector == "A"] = "T"
  newDNAstrVector[DNAstrVector == "T"] = "A"
  newDNAstrVector[DNAstrVector == "C"] = "G"
  newDNAstrVector[DNAstrVector == "G"] = "C"
  return(rev(newDNAstrVector))
}

###
print("set up data")
###
flanking_range = 7500
flanging_skip = 2500

chr_bins <- 40

padding <- flanging_skip * 2 + 500000 # padding is double the flanging_skip and extra bp of 1 Mb, the extra can be modified

max_unpaired_distance = flanging_skip * 2 + 500000

sample_array_pairs <- 50
sample_repeats_per_array <- 20

keep_chrs_austr <- c(1,4,7)
keep_chrs_austr2 <- c(1,4,7)



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

assemblies <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", recursive = T, full.names = T)
assembly_A <- assemblies[grep(genome_A, assemblies)[1]]
assembly_A <- read.fasta(file = assembly_A)
assembly_B <- assemblies[grep(genome_B, assemblies)[1]]
assembly_B <- read.fasta(file = assembly_B)


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

if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
  keep_chrs <- keep_chrs_austr
  arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
  repeats_1 <- repeats_1[repeats_1$seqID %in% chr_1_names[keep_chrs],]
  synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names[keep_chrs],]
  chr_1_sizes <- chr_1_sizes[keep_chrs]
  chr_1_names <- chr_1_names[keep_chrs]
  chr_1_sum_sizes <- chr_1_sum_sizes[keep_chrs]
  chr_1_sum_starts <- chr_1_sum_starts[keep_chrs]
  chr_1_sum_ends <- chr_1_sum_ends[keep_chrs]
} 
if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
  keep_chrs <- keep_chrs_austr
  arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
  repeats_2 <- repeats_2[repeats_2$seqID %in% chr_2_names[keep_chrs],]
  synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names[keep_chrs],]
  chr_2_sizes <- chr_2_sizes[keep_chrs]
  chr_2_names <- chr_2_names[keep_chrs]
  chr_2_sum_sizes <- chr_2_sum_sizes[keep_chrs]
  chr_2_sum_starts <- chr_2_sum_starts[keep_chrs]
  chr_2_sum_ends <- chr_2_sum_ends[keep_chrs]
} 
if(genome_A == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
  keep_chrs <- keep_chrs_austr2
  arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names[keep_chrs],]
  repeats_1 <- repeats_1[repeats_1$seqID %in% chr_1_names[keep_chrs],]
  synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_A %in% chr_1_names[keep_chrs],]
  chr_1_sizes <- chr_1_sizes[keep_chrs]
  chr_1_names <- chr_1_names[keep_chrs]
  chr_1_sum_sizes <- chr_1_sum_sizes[keep_chrs]
  chr_1_sum_starts <- chr_1_sum_starts[keep_chrs]
  chr_1_sum_ends <- chr_1_sum_ends[keep_chrs]
} 
if(genome_B == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
  keep_chrs <- keep_chrs_austr2
  arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names[keep_chrs],]
  repeats_2 <- repeats_2[repeats_2$seqID %in% chr_2_names[keep_chrs],]
  synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$chromosome_B %in% chr_2_names[keep_chrs],]
  chr_2_sizes <- chr_2_sizes[keep_chrs]
  chr_2_names <- chr_2_names[keep_chrs]
  chr_2_sum_sizes <- chr_2_sum_sizes[keep_chrs]
  chr_2_sum_starts <- chr_2_sum_starts[keep_chrs]
  chr_2_sum_ends <- chr_2_sum_ends[keep_chrs]
} 



### rearrange data
pair_A <- which(unlist(lapply(species, function(x) {grepl(x, genome_A)})))
pair_B <- which(unlist(lapply(species, function(x) {grepl(x, genome_B)})))


# adjust starts of synteny data frame (using rearranged chr sizes vectors)

for(i in seq_along(chr_1_sizes)) {
  synteny_pairs_df$flanked_reg_start_A[synteny_pairs_df$chromosome_A == chr_1_names[i]] <- synteny_pairs_df$flanked_reg_start_A[synteny_pairs_df$chromosome_A == chr_1_names[i]] - chr_1_sum_sizes[i]
  synteny_pairs_df$flanked_reg_end_A[synteny_pairs_df$chromosome_A == chr_1_names[i]] <- synteny_pairs_df$flanked_reg_end_A[synteny_pairs_df$chromosome_A == chr_1_names[i]] - chr_1_sum_sizes[i]
}
for(i in seq_along(chr_2_sizes)) {
  synteny_pairs_df$flanked_reg_start_B[synteny_pairs_df$chromosome_B == chr_2_names[i]] <- synteny_pairs_df$flanked_reg_start_B[synteny_pairs_df$chromosome_B == chr_2_names[i]] - chr_2_sum_sizes[i]
  synteny_pairs_df$flanked_reg_end_B[synteny_pairs_df$chromosome_B == chr_2_names[i]] <- synteny_pairs_df$flanked_reg_end_B[synteny_pairs_df$chromosome_B == chr_2_names[i]] - chr_2_sum_sizes[i]
}


# make dotplots

wsize = 12
reversed = F

synteny_pairs_df$region_a_size = synteny_pairs_df$flanked_reg_end_A - synteny_pairs_df$flanked_reg_start_A + 1
synteny_pairs_df$region_b_size = synteny_pairs_df$flanked_reg_end_B - synteny_pairs_df$flanked_reg_start_B + 1

flank_plot = 7500
max_seq_to_plot <- 75000 # 7 is 226 kbp

# synteny_pairs_df = synteny_pairs_df[!is.na(synteny_pairs_df$island_B_ID) & synteny_pairs_df$island_B_ID == 13,]

for(j in 1 : nrow(synteny_pairs_df)) {
# for(j in 1 : 7) {
  
  
  start_A = synteny_pairs_df$flanked_reg_start_A[j] - flank_plot
  end_A = synteny_pairs_df$flanked_reg_end_A[j] + flank_plot
  
  start_B = synteny_pairs_df$flanked_reg_start_B[j] - flank_plot
  end_B = synteny_pairs_df$flanked_reg_end_B[j] + flank_plot
  
  seq_A_len <- end_A - start_A + 1 
  seq_B_len <- end_B - start_B + 1 
  
  len_A <- chr_1_sizes[which(chr_1_names == synteny_pairs_df$chromosome_A[j])]
  len_B <- chr_2_sizes[which(chr_2_names == synteny_pairs_df$chromosome_B[j])]
  
  cat(j, "/", nrow(synteny_pairs_df), seq_A_len/1000, "kb vs", seq_B_len/1000, "kb,")
  
  if(seq_A_len > max_seq_to_plot) {cat("seq_A_len is longer than allowed to calculate, not plotting \n") next}
  if(seq_B_len > max_seq_to_plot) {cat("seq_B_len is longer than allowed to calculate, not plotting \n") next}
  
  if(start_A < 1) start_A = 1
  if(end_A > len_A) end_A = len_A
  
  if(start_B < 1) start_B = 1
  if(end_B > len_B) end_B = len_B
  
  
  tyba_A_chr <- repeats_1[repeats_1$seqID == synteny_pairs_df$chromosome_A[j],]
  tyba_B_chr <- repeats_2[repeats_2$seqID == synteny_pairs_df$chromosome_B[j],]
  
  tyba_A_chr <- tyba_A_chr[tyba_A_chr$end > start_A & tyba_A_chr$start < end_A,]
  tyba_B_chr <- tyba_B_chr[tyba_B_chr$end > start_B & tyba_B_chr$start < end_B,]
  
  
  
  sequence_a = assembly_A[[synteny_pairs_df$chromosome_A[j]]][start_A : end_A]
  
  
  tyba_A_chr$start <- tyba_A_chr$start - start_A
  tyba_A_chr$end <- tyba_A_chr$end - start_A
  tyba_B_chr$start <- tyba_B_chr$start - start_B
  tyba_B_chr$end <- tyba_B_chr$end - start_B
  
  
  sequence_b = assembly_B[[synteny_pairs_df$chromosome_B[j]]][start_B : end_B]
  
  # write.fasta(sequences = sequence_a, 
  #             names = paste0(input_id, "_", j, "_", species[pair_A], "_", synteny_pairs_df$chromosome_A[j]), 
  #             file.out = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/33_dotplots/", input_id, "_", j, "_", species[pair_A], "_", synteny_pairs_df$chromosome_A[j], ".fasta"))
  # write.fasta(sequences = sequence_b, 
  #             names = paste0(input_id, "_", j, "_", species[pair_B], "_", synteny_pairs_df$chromosome_B[j]),  
  #             file.out = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/33_dotplots/", input_id, "_", j, "_", species[pair_B], "_", synteny_pairs_df$chromosome_B[j], ".fasta"))
  
  timea = Sys.time()
  
  # Build everything into p as before
  plot_forward <- dotPlotg(paste(toupper(sequence_a), collapse = ""), 
                           paste(toupper(sequence_b), collapse = ""), 
                           wsize = wsize)
  
  plot_reverse <- dotPlotg(paste(toupper(sequence_a), collapse = ""), 
                           paste(revCompString(sequence_b), collapse = ""), 
                           wsize = wsize)
  
  built_forward <- ggplot_build(plot_forward)
  built_reverse <- ggplot_build(plot_reverse)
  
  forward_data <- built_forward$data[[1]]
  reverse_data <- built_reverse$data[[1]]
  
  
  # Start with base plot
  p <- ggplot() +
    theme_bw() +
    labs(x = paste(synteny_pairs_df$chr_A[j], j), 
         y = paste(synteny_pairs_df$chr_B[j], j), 
         title = paste(synteny_pairs_df$genome_A[j], j)) + 
    geom_vline(xintercept = c(0, 0 + seq_A_len)) + 
    geom_hline(yintercept = c(0, 0 + seq_B_len)) +
    scale_x_continuous(limits = c(-length(sequence_a)/100, length(sequence_a))) +
    scale_y_continuous(limits = c(-length(sequence_b)/100, length(sequence_b))) +
    coord_equal()
  
  # Conditionally add forward data
  if (nrow(forward_data) > 0) {
    p <- p + geom_point(data = forward_data, aes(x = x, y = y), color = "blue", alpha = 0.5, size = 0.5)
  }
  
  # Conditionally add reverse data
  if (nrow(reverse_data) > 0) {
    reverse_data$y <- length(sequence_b) - reverse_data$y
    p <- p + geom_point(data = reverse_data, aes(x = x, y = y), color = "red", alpha = 0.5, size = 0.5)
  }
  
  # Conditionally add polygons
  if (nrow(tyba_A_chr) > 0) {
    p <- p + geom_polygon(data = data.frame(
      x = c(rbind(tyba_A_chr$start, tyba_A_chr$end, tyba_A_chr$end, tyba_A_chr$start)),
      y = c(rbind(rep(0, nrow(tyba_A_chr)), 
                  rep(0, nrow(tyba_A_chr)), 
                  seq_B_len, 
                  seq_B_len)),
      group = rep(1:nrow(tyba_A_chr), each = 4)),
      aes(x = x, y = y, group = group), 
      fill = "green", alpha = 0.1)
  }
  
  if (nrow(tyba_B_chr) > 0) {
    p <- p + geom_polygon(data = data.frame(
      y = c(rbind(tyba_B_chr$start, tyba_B_chr$end, tyba_B_chr$end, tyba_B_chr$start)),
      x = c(rbind(rep(0, nrow(tyba_B_chr)), 
                  rep(0, nrow(tyba_B_chr)), 
                  seq_A_len, 
                  seq_A_len)),
      group = rep(1:nrow(tyba_B_chr), each = 4)),
      aes(x = x, y = y, group = group), 
      fill = "green", alpha = 0.1)
  }
  
  
  
  heightp = 15
  widthp = 15
  # widthp <- heightp * seq_A_len/seq_B_len
  # if(widthp > 49) {
  #   widthp = 49
  #   heightp = widthp * seq_B_len/seq_A_len
  # }
  
  if(!synteny_pairs_df$is_island_A[j] | !synteny_pairs_df$is_island_B[j]) {
    ggsave(
      filename = paste0("uneven_", species[pair_A], "_", synteny_pairs_df$chromosome_A[j], "_", species[pair_B], "_", synteny_pairs_df$chromosome_B[j],
                        "_pair_", j, ".png"),
      device = NULL,
      path = "/home/pwlodzimierz/Rhynchospora/upload_data/33_dotplots/",
      scale = 1,
      width = widthp,
      height = heightp
    )
  } else {
    ggsave(
      filename = paste0("even_", species[pair_A], "_", synteny_pairs_df$chromosome_A[j], "_", species[pair_B], "_", synteny_pairs_df$chromosome_B[j],
                        "_pair_", j, ".png"),
      device = NULL,
      path = "/home/pwlodzimierz/Rhynchospora/upload_data/33_dotplots/",
      scale = 1,
      width = widthp,
      height = heightp
    )
  }
  
  
  
  
  print(difftime(Sys.time(), timea))
  
}













