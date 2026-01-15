# not used in the manuscript, but used to check island and spacing sizes

#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

suppressMessages(library(seqinr))
suppressMessages(library(msa))
suppressMessages(library(IRanges))


def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")


bins_per_feature <- 100


if(F) { # go to line 79
  repeats <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_repeats_filtered.csv")
  edta <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta_edta_filtered.csv")
  classes <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_classes_merged_filtered.csv")
  islands_genome <- read.csv(("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/upload_data/islands/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_islands_genome.csv"))
  helixer <- read.table(file = "C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/helixer_files/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.helixer.gff3", header = F, skip = 4, sep = "\t")
  chr_no_sizes <- read.csv(file = "C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/chr.no.and.sizes.full.csv")
  
  no_helixer = FALSE
  no_edta = FALSE
  
  
  setwd("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data")
  source("C:/Users/Piotr Włodzimierz/Desktop/ToL/temp_data/aux_fun.R")
  
  
  fasta_name <- "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta"
  assembly_name <- "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr"
  
  i=24
}



data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

chr_no_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 30
print(i)

{
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  print(data_directories[i])
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); quit(save = "no")}
  
  no_edta = FALSE
  edta_file = list.files(pattern = "_edta_filtered.csv", full.names = TRUE)
  if(length(edta_file) != 1) {print(paste0(i, " no edta!")); no_edta = TRUE}
  
  islands_genome_file <- list.files(pattern = "_islands_genome", full.names = TRUE)
  if(length(islands_genome_file) == 0) {print(paste0(i, " no islands!")); quit(save = "no")}
  if(length(islands_genome_file) > 1) {islands_genome_file <- list.files(pattern = "_islands_genome.csv", full.names = TRUE)}
  
  no_helixer = FALSE
  helixer_file <- list.files(pattern = "_helixer_filtered.gff", full.names = TRUE)
  if(length(helixer_file) == 0) {print(paste0(i, " no helixer!")); no_helixer = TRUE} 
  
  assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  fasta_name = strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2]
  
  print(assembly_name)
  
  repeats = read.csv(file = repeat_file)
  if(!no_edta) edta = read.csv(file = edta_file)
  islands_genome <- read.csv(file = islands_genome_file)
  if(!no_helixer) helixer <- read.table(file = helixer_file, header = F, skip = 4, sep = "\t")
    
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  
  chr_no_sizes <- chr_no_sizes[chr_no_sizes$assembly.name == fasta_name,]
  if(nrow(chr_no_sizes) == 0) {print(paste0(i, " no sequence lengths!")); quit(save = "no")}
  
  chromosomes <- chr_no_sizes$chromosome.name
  
  spacings_bins_bp <- rep(0, bins_per_feature)
  spacings_bins_genes_bp <- rep(0, bins_per_feature)
  spacings_bins_TEs_bp <- rep(0, bins_per_feature)
  spacings_bins_LTRs_bp <- rep(0, bins_per_feature)
  spacings_bins_repeats_bp <- rep(0, bins_per_feature)
  
  islands_bins_bp <- rep(0, bins_per_feature)
  islands_bins_genes_bp <- rep(0, bins_per_feature)
  islands_bins_TEs_bp <- rep(0, bins_per_feature)
  islands_bins_LTRs_bp <- rep(0, bins_per_feature)
  islands_bins_Tyba_bp <- rep(0, bins_per_feature)
  
  islands_bins_Tyba_distance_to_island_consensus_bp <- rep(0, bins_per_feature)
  
  spacing_sizes <- NULL
  
  for(j in seq_along(chromosomes)) {
    
    chr_islands <- islands_genome[islands_genome$chromosome == chromosomes[j],]
    chr_repeats <- repeats[repeats$seqID == chromosomes[j],]
    chr_tyba <- tyba[tyba$seqID == chromosomes[j],]
    if(!no_helixer) chr_helixer <- helixer[helixer$V1 == chromosomes[j],]
    if(!no_edta) chr_edta <- edta[edta$V1 == chromosomes[j],]
    if(!no_edta) chr_LTR <- chr_edta[grep("LTR", chr_edta$V3),]
    
    if(!no_helixer) chr_helixer <- chr_helixer[chr_helixer$V3 == "gene",]
    
    chr_tyba$midpoint <- round(chr_tyba$start + (chr_tyba$end - chr_tyba$start) / 2)
    
    spacings <- data.frame(start = c(1, (chr_islands$end + 1)))
    spacings$ends <- c((chr_islands$start - 1), chr_no_sizes$size[j])
    spacings$size <- spacings$ends - spacings$start + 1
    spacing_sizes <- c(spacing_sizes, spacings$size)
    
    for(k in 1 : bins_per_feature) {
      cat("chr", j, "/", length(chromosomes), "bin", k, "/", bins_per_feature, "\n")

      ### spacings
      spacings_bins_starts <- unlist(lapply(1 : nrow(spacings), function(X)
        round(spacings$start[X] + (spacings$size[X] / bins_per_feature) * (k-1))
      ))
      spacings_bins_ends <- unlist(lapply(1 : nrow(spacings), function(X)
        round(spacings$start[X] + (spacings$size[X] / bins_per_feature) * k)
      ))

      spacings_bins_bp[k] <- spacings_bins_bp[k] + sum(spacings_bins_ends - spacings_bins_starts + 1)
      
      if(!no_helixer) {
        overlaps <- overlapsRanges(IRanges(chr_helixer$V4, chr_helixer$V5), IRanges(spacings_bins_starts, spacings_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        spacings_bins_genes_bp[k] <- spacings_bins_genes_bp[k] + overlap_coords
      }
      
      
      if(!no_edta) {
        overlaps <- overlapsRanges(IRanges(chr_edta$V4, chr_edta$V5), IRanges(spacings_bins_starts, spacings_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        spacings_bins_TEs_bp[k] <- spacings_bins_TEs_bp[k] + overlap_coords
        
        overlaps <- overlapsRanges(IRanges(chr_LTR$V4, chr_LTR$V5), IRanges(spacings_bins_starts, spacings_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        spacings_bins_LTRs_bp[k] <- spacings_bins_LTRs_bp[k] + overlap_coords
      }
      
      overlaps <- overlapsRanges(IRanges(chr_repeats$start, chr_repeats$end), IRanges(spacings_bins_starts, spacings_bins_ends))
      overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
      spacings_bins_repeats_bp[k] <- spacings_bins_repeats_bp[k] + overlap_coords

      ### islands
      islands_bins_starts <- unlist(lapply(1 : nrow(chr_islands), function(X)
        round(chr_islands$start[X] + (chr_islands$total_bp[X] / bins_per_feature) * (k-1))
      ))
      islands_bins_ends <- unlist(lapply(1 : nrow(chr_islands), function(X)
        round(chr_islands$start[X] + (chr_islands$total_bp[X] / bins_per_feature) * k)
      ))

      islands_bins_bp[k] <- islands_bins_bp[k] + sum(islands_bins_ends - islands_bins_starts + 1)
      
      if(!no_helixer) {
        overlaps <- overlapsRanges(IRanges(chr_helixer$V4, chr_helixer$V5), IRanges(islands_bins_starts, islands_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        islands_bins_genes_bp[k] <- islands_bins_genes_bp[k] + overlap_coords
      }
      
      if(!no_edta) {
        overlaps <- overlapsRanges(IRanges(chr_edta$V4, chr_edta$V5), IRanges(islands_bins_starts, islands_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        islands_bins_TEs_bp[k] <- islands_bins_TEs_bp[k] + overlap_coords
        
        overlaps <- overlapsRanges(IRanges(chr_LTR$V4, chr_LTR$V5), IRanges(islands_bins_starts, islands_bins_ends))
        overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
        islands_bins_LTRs_bp[k] <- islands_bins_LTRs_bp[k] + overlap_coords
      }
      
      overlaps <- overlapsRanges(IRanges(tyba$start, tyba$end), IRanges(islands_bins_starts, islands_bins_ends))
      overlap_coords <- length(unique(unlist(lapply(seq_along(overlaps), function(X) return(start(overlaps)[X] : end(overlaps)[X])))))
      islands_bins_Tyba_bp[k] <- islands_bins_Tyba_bp[k] + overlap_coords

      ### island Tyba similarity

      dist_sum <- unlist(lapply(seq_along(islands_bins_starts), function(X) {
        tyba_island_bin <- chr_tyba[chr_tyba$midpoint %in% (islands_bins_starts[X] : islands_bins_ends[X]),]
        if(nrow(tyba_island_bin) == 0) return(0)
        consensus = tolower(chr_islands$tyba_consensus[X])
        distances <- stringdist::stringdistmatrix(consensus, tyba_island_bin$sequence, method = "lv")
        return(sum(distances))
      } ))

      islands_bins_Tyba_distance_to_island_consensus_bp[k] <- islands_bins_Tyba_distance_to_island_consensus_bp[k] + sum(dist_sum)

    }
    
  }
  
  
  island_sizes <- islands_genome$end - islands_genome$start
  
  mean_island_size <- mean(island_sizes)
  sd_island_size <- sd(island_sizes)
  
  max_island_size <- 150000
  
  mean_spacing_size <- mean(spacing_sizes)
  sd_spacing_size <- sd(spacing_sizes)
  
  max_spacing_size <- 1200000
  
  if(max_island_size < mean_island_size*1.1) max_island_size = mean_island_size*1.1
  if(max_spacing_size < mean_spacing_size*1.1) max_spacing_size = mean_spacing_size*1.1
  
  pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/average_island_spacing_properties/ave_isl_spacing_", assembly_name, "_plot.pdf"), width = 6, height = 2)
  # pdf(file = paste0("C:\\Users\\Piotr Włodzimierz\\Desktop\\Rhynchospora\\temp_data/ave_isl_spacing_", assembly_name, "_2plot.pdf"), width = 6, height = 2)
  par(mfrow = c(1,2), mar = c(3,1,2,0), oma = c(0,0,0,4), xpd=NA)
  
  # nf <- layout( # equal layout, islands and spacings take 50% of the plot width
  #   matrix(c(1,2,3,4), ncol=4, byrow=TRUE),
  #   widths=c(mean_spacing_size/max_spacing_size, (max_spacing_size-mean_spacing_size)/max_spacing_size,
  #            mean_island_size/max_island_size , (max_island_size-mean_island_size)/max_island_size), # this modifies how much space each plot will take
  #   heights=c(1)
  # )
  # nf <- layout( # proportional layout, islands and spacings take the same plot width fraction per bp
  #   matrix(c(1,2,3,4), ncol=4, byrow=TRUE),
  #   widths=c(mean_spacing_size, (max_spacing_size-mean_spacing_size),
  #            mean_island_size , (max_island_size-mean_island_size)), # this modifies how much space each plot will take
  #   heights=c(1)
  # )
  
  
  plot(100 * spacings_bins_genes_bp / spacings_bins_bp, type = "l", ylim = c(0,110), col = "#00ee0080", 
       ylab = "", xlab = "", lwd = 3, font = 2, xaxt = "n", yaxt = "n")
  # axis(1, labels = FALSE) 
  points(100 * spacings_bins_TEs_bp / spacings_bins_bp, type = "l", ylim = c(0,110), col = "#0000ee80", lwd = 3)
  # points(100 * spacings_bins_LTRs_bp / spacings_bins_bp, type = "l", ylim = c(0,110), col = "#eeee0080", lwd = 3)
  points(100 * spacings_bins_repeats_bp / spacings_bins_bp, type = "l", ylim = c(0,110), col = "#ee000080", lwd = 3) 
  # legend("topleft", legend = c("genes", "TEs", "LTRs", "repeats"), col = c("#00ee88", "#0088ee", "#880088", "#88ee00"), lty = 1, bty = "n")
  # mtext("Feature % coverage", side = 2, line = 2.2) 
  # mtext(paste0(bins_per_feature, " spacing bins, scaled according to mean spacing size"), side = 1, line = 2.2) 
  title <- paste0(assembly_name, ".   spacings mean Kbp: ", round(mean_spacing_size/1000), "+-", round(sd_spacing_size/1000), ".   Island mean Kbp: ", round(mean_island_size/1000), "+-", round(sd_island_size/1000))
 
  # mtext(title, 
        # side = 3, line = 0.1, adj = 0, font = 2) 
  
  lines(x = c(0,100/(round(mean_spacing_size/1000)/100)), y = c(103,103), lwd = 3, lend = 1)
  # text(x = 100/(round(mean_spacing_size/1000)/100), y = 103, substitute(paste(bold('100 Kbp'))), adj = 0, cex = 1.2)
  
  # plot(NA,NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
  
  plot(100 * islands_bins_genes_bp / islands_bins_bp, type = "l", ylim = c(0,110), col = "#00ee0080", 
       ylab = "", xlab = "", lwd = 3, font = 2, yaxt = "n", xaxt = "n")
  # axis(1, labels = FALSE) 
  axis(2, labels = FALSE) 
  points(100 * islands_bins_TEs_bp / islands_bins_bp, type = "l", ylim = c(0,110), col = "#0000ee80", lwd = 3)
  # points(100 * islands_bins_LTRs_bp / islands_bins_bp, type = "l", ylim = c(0,110), col = "#eeee0080", lwd = 3)
  points(100 * islands_bins_Tyba_bp / islands_bins_bp, type = "l", ylim = c(0,110), col = "#ee006680", lwd = 3)
  # mtext("Feature % coverage", side = 2, line = 2.2) 
  # mtext(paste0(bins_per_feature, " island bins, scaled according to mean spacing size"), side = 1, line = 2.2)
  
  lines(x = c(0,100/(round(mean_island_size/1000)/10)), y = c(105,105), lwd = 3, lend = 1)
  # text(x = 100/(round(mean_island_size/1000)/10), y = 105, substitute(paste(bold('10 Kbp'))), adj = 0, cex = 1.2)
  
  par(new = TRUE)
  
  plot(100 * islands_bins_Tyba_distance_to_island_consensus_bp / islands_bins_Tyba_bp, type = "l", ylim = c(5,40), axes = FALSE, 
       xlab = "", ylab = "", col = "#ee00ee", lty = 6, lwd = 2, font = 2, xaxt = "n")
  # axis(side = 4, at = pretty(range(100 * islands_bins_Tyba_distance_to_island_consensus_bp / islands_bins_Tyba_bp)))
  # axis(side = 4, col = "#ee00ee", font = 2)
  # mtext("Tyba dist to cons, %", side = 4, col = "#ee00ee", line = 2.2) 
  
  
  # plot(NA,NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
  # legend("topright", legend = c("repeats/Tyba", "TEs", "LTRs", "genes", "Tyba dist, right Yax"), col = c("#ee000080", "#ee00ee80", "#eeee0080", "#00ee88", "#ee0000"), 
  #        lty = c(1,1,1,1,2), bty = "n", cex = 0.5)
  
  dev.off()
  
  
  max_space_between_repeats_to_ignore = 50000
  # if(grepl("alba", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("austrobrasiliensis", assembly_name)) max_space_between_repeats_to_ignore <- 100000
  # if(grepl("barbata", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("breviuscula", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("cephalotes", assembly_name)) max_space_between_repeats_to_ignore <- 30000
  # if(grepl("ciliata", assembly_name)) max_space_between_repeats_to_ignore <- 80000
  # if(grepl("colorata", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("corymbosa", assembly_name)) max_space_between_repeats_to_ignore <- 30000
  # if(grepl("gaudichaudii", assembly_name)) max_space_between_repeats_to_ignore <- 100000
  # if(grepl("holoschoenoides", assembly_name)) max_space_between_repeats_to_ignore <- 70000
  # if(grepl("nervosa", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("pubera", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("radicans", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("ridleyi", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("riparia", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  # if(grepl("rugosa", assembly_name)) max_space_between_repeats_to_ignore <- 70000
  # if(grepl("tenerrima", assembly_name)) max_space_between_repeats_to_ignore <- 60000
  # if(grepl("tenuis", assembly_name)) max_space_between_repeats_to_ignore <- 100000
  # if(grepl("watsonii", assembly_name)) max_space_between_repeats_to_ignore <- 50000
  
  
  
  
  
  pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/average_island_spacing_properties/island_spacing_size_", assembly_name, "_histograms.pdf"), width = 10, height = 14)
  # pdf(file = paste0("C:\\Users\\Piotr Włodzimierz\\Desktop\\Rhynchospora\\temp_data/ave_isl_spacing_", assembly_name, "_histograms.pdf"), width = 10, height = 14)
  par(mfrow = c(5,1))
  hist(island_sizes, breaks = 500)
  hist(island_sizes[island_sizes < 100000 & island_sizes > 0], breaks = seq(0,100000, length.out = 100), xlim = c(0, 100000))
  hist(spacing_sizes, breaks = 500)
  abline(v = max_space_between_repeats_to_ignore, lty = 2, col = "red")
  hist(spacing_sizes[spacing_sizes < 1000000 & spacing_sizes > 0], breaks = seq(0,1000000, length.out = 100), xlim = c(0, 1000000))
  abline(v = max_space_between_repeats_to_ignore, lty = 2, col = "red")
  hist(spacing_sizes[spacing_sizes < 200000 & spacing_sizes > 0], breaks = seq(0,200000, length.out = 100), xlim = c(0, 200000))
  abline(v = max_space_between_repeats_to_ignore, lty = 2, col = "red")
  dev.off()
  
  
}



