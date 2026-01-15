#!/usr/bin/env Rscript

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")




data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)

assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

data_directories <- data_directories[grep(paste(c("Rrugosa.chr.fasta",
                                                  "Ralba.chr.fasta",
                                                  "Rcephalotes.hap2.chr.fasta", 
                                                  "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta",
                                                  "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta",
                                                  "R_pubera_ref_2n10_v2.chr.fasta",
                                                  "Rbreviuscula.hap1.chr.fasta", 
                                                  "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                  "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                  "Rcolorata.hap1.chr.fasta",
                                                  "Rtenerrima.chr.fasta",
                                                  "Rhync_filiformis_5519G.hap1.chr.fasta",
                                                  "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                  "Rhync_tenuis_ref.hap1.chr.fasta",
                                                  "Rhync_riparia_5519C.hap1.chr.fasta",
                                                  "R_barbata_hap1_chrs.fa",
                                                  "R_holoschoenoides_hap1_chr.fasta",
                                                  "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta"), collapse = "|"), data_directories)]


main_data <- data.frame()

for(i in 1 : length(data_directories))  {
  print(paste0(i, " / ", length(data_directories), ": ", basename(data_directories[i])))
  setwd(data_directories[i])
  all_repeats = read.csv(file = list.files(pattern = "repeats_filtered.csv"))

  arrays <- read.csv(file = list.files(pattern = "_arrays_filtered.csv"))

  hor_files <- list.files(path = "./tyba_island_split_HORs/", pattern = "HORs_Tyba_", full.names = T)
  repeat_files <- list.files(path = "./tyba_island_split_HORs/", pattern = "repeats_with_hors_Tyba_", full.names = T)

  windows <- 25
  score_per_window <- rep(0, windows)
  repeats_per_window <- rep(0, windows)
  new_scores <- NULL
  cat("\n")
  cat("=========================\n")
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

    new_scores <- c(new_scores, repeats$novel_score)

  }
  cat("\n")

  assembly <- seqinr::read.fasta(file = assembly_files[grepl(basename(data_directories[i]), assembly_files)])


  # assembly_name
  assembly_name = ""
  assembly_name = basename(data_directories[i])

  # haplotype_info
  haplotype_info = ""

  # genome_size
  genome_size = 0
  genome_size <- sum(unlist(lapply(assembly, length)))

  # chromosome_number
  chromosome_number = 0
  chromosome_number = length(assembly)

  # total_repeats_number
  total_repeats_number = 0
  total_repeats_number = nrow(all_repeats)

  # total_repeats_bp
  total_repeats_bp = 0
  total_repeats_bp = sum(all_repeats$width)

  # main_tyba_variant
  main_tyba_variant = ""

  # total_tyba_number
  total_tyba_number = 0
  
  tyba = all_repeats[grepl("yba", all_repeats$new_class),]

  total_tyba_number <- nrow(tyba)

  # tyba_total_bp
  tyba_total_bp = 0
  tyba_total_bp = sum(tyba$width)
  
  # tyba_monomer_mean_length 
  tyba_monomer_mean_length = 0
  tyba_monomer_mean_length = mean(tyba$width[tyba$width > 100 & tyba$width < 250])
  # * there is a number of repeats in non-monomer state, hence filtering
  

  # tyba_arrays_number
  tyba_arrays_number = 0
  tyba_arrays_number = nrow(arrays)

  # mean_array_bp
  mean_array_bp = 0
  arrays$width <- arrays$end - arrays$start + 1
  mean_array_bp = mean(arrays$width)


  # mean_spacing_bp
  mean_spacing_bp = 0
  arrays$dist_to_next <- 0
  arrays$dist_to_next[1: (nrow(arrays) - 1)] <- arrays$start[2 : nrow(arrays)] - arrays$end[1: (nrow(arrays) - 1)]

  mean_spacing_bp = mean(arrays$dist_to_next[arrays$dist_to_next > 0])

  # mean_Tyba_HOR_score
  mean_Tyba_HOR_score = 0
  mean_Tyba_HOR_score = mean(new_scores)

  main_data <- rbind(main_data, data.frame(assembly_name = assembly_name,
                                           haplotype_info = haplotype_info,
                                           genome_size = genome_size,
                                           chromosome_number = chromosome_number,
                                           total_repeats_number = total_repeats_number,
                                           total_repeats_bp = total_repeats_bp,
                                           main_tyba_variant = main_tyba_variant,
                                           total_tyba_number = total_tyba_number,
                                           tyba_total_bp = tyba_total_bp,
                                           tyba_monomer_mean_length = tyba_monomer_mean_length,
                                           tyba_arrays_number = tyba_arrays_number,
                                           mean_array_bp = mean_array_bp,
                                           mean_spacing_bp = mean_spacing_bp,
                                           mean_Tyba_HOR_score = mean_Tyba_HOR_score))
  
}

write.csv(x = main_data, file = "/home/pwlodzimierz/Rhynchospora/upload_data/S1_table.csv", row.names = F)











