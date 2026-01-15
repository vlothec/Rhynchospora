
colors_vector <- c("#E69F00", "#D55E00", "#D95F02", "#FC8D62", "#E5C494", 
                   "#F0E442", "#FFD92F", "#A6D854", "#66C2A5", "#009E73", 
                   "#1B9E77", "#66A61E", "#CC79A7", "#E78AC3", "#E7298A", 
                   "#B3B3CC", "#8DA0CB", "#7570B3", "#56B4E9", "#0072B2")

custom_labels <- c(
  expression(italic("R. rugosa")), 
  expression(italic("R. alba")), 
  expression(italic("R. gaudichaudii")), 
  expression(italic("R. cephalotes")), 
  expression(italic("R. ridleyi")), 
  expression(italic("R. watsonii")), 
  expression(italic("R. radicans")), 
  expression(italic("R. pubera")), 
  expression(italic("R. breviuscula")), 
  expression(italic("R. nervosa")), 
  expression(italic("R. ciliata")), 
  expression(italic("R. colorata")), 
  expression(italic("R. tenerrima")), 
  expression(italic("R. filiformis")), 
  expression(italic("R. austrobrasiliensis")), 
  expression(italic("R. tenuis")), 
  expression(italic("R. riparia")), 
  expression(italic("R. barbata")), 
  expression(italic("R. holoschoenoides")), 
  expression(italic("R. corymbosa"))
)
custom_order <- c(
  "rugosa", 
  "alba", 
  "gaudichaudii", 
  "cephalotes", 
  "ridleyi", 
  "watsonii", 
  "radicans", 
  "pubera", 
  "breviuscula", 
  "nervosa", 
  "ciliata", 
  "colorata", 
  "tenerrima", 
  "filiformis", 
  "austrobrasiliensis", 
  "tenuis", 
  "riparia", 
  "barbata", 
  "holoschoenoides", 
  "corymbosa"
)

# Define point shapes (pch values) to improve distinction
pch_vector <- c(15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18)


# Step 0: Load the repeats
# Sample if necessary

data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)

data_directories <- data_directories[grepl(pattern = paste("Rrugosa.chr.fasta", 
                                                           "R_barbata_hap1_chrs.fa",
                                                           "R_holoschoenoides_hap1_chr.fasta",
                                                           "R_pubera_ref_2n10_v2.chr.fasta",
                                                           "Ralba.chr.fasta",
                                                           "Rbreviuscula.hap1.chr.fasta", 
                                                           "Rcephalotes.hap1.chr.fasta", 
                                                           "Rcolorata.hap1.chr.fasta", 
                                                           "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 
                                                           "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 
                                                           "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta", 
                                                           "Rhync_filiformis_5519G.hap1.chr.fasta", 
                                                           "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_radicans.asm.bp.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 
                                                           "Rhync_riparia_5519C.hap1.chr.fasta", 
                                                           "Rtenuis.hap1.chr.fasta", 
                                                           "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta", 
                                                           "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta", 
                                                           "Rtenerrima.chr.fasta",
                                                           sep = "|"), x = data_directories)]
data_directories_sorted <- NULL

for(i in seq_along(custom_order)) {
  data_directories_sorted <- c(data_directories_sorted, data_directories[grep(custom_order[i], data_directories)])
}

tyba_all <- NULL
for(i in seq_along(data_directories_sorted)) {
  cat("reading in assembly from", data_directories_sorted[i], "\n")
  setwd(data_directories_sorted[i])
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  repeats = read.csv(file = repeat_file)
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  tyba$assembly <- strsplit(data_directories_sorted[i], split = "/")[[1]][6]
  tyba$width <- tyba$end - tyba$start + 1
  tyba_all <- rbind(tyba_all, tyba[,c(1,8,14)])
}
tyba_all_b <- tyba_all

setwd("/home/pwlodzimierz/Rhynchospora/upload_data/histogram_all_tyba")

xlims = c(150,200)

tyba_all <- tyba_all[tyba_all$width < xlims[2],]
tyba_all <- tyba_all[tyba_all$width > xlims[1],]

hist(tyba_all$width, breaks = xlims[1]:xlims[2])

for(i in seq_along(unique(tyba_all$assembly))) {
  hist_data <- hist(tyba_all$width[tyba_all$assembly == unique(tyba_all$assembly)[i]], breaks = xlims[1]:xlims[2], plot = FALSE)
  if(i == 1) {
    cum_hist_data <- list(hist_data$counts)
  } else {
    cum_hist_data[[i]] = (cum_hist_data[[i-1]] + hist_data$counts)
  }
}

pdf("tyba_histogram_one_per_species.pdf")
for(i in length(cum_hist_data) : 1) {
  if(i == 20)  {
    plot(hist_data$mids, cum_hist_data[[i]], type = "h", col = colors_vector[i], lend = 1, lwd = 10)
  } else {
    points(hist_data$mids, cum_hist_data[[i]], type = "h", col = colors_vector[i], lend = 1, lwd = 10)
  }
}
dev.off()


###############
#
###############
#
###############



data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories_sorted <- NULL

for(i in seq_along(custom_order)) {
  data_directories_sorted <- c(data_directories_sorted, data_directories[grep(custom_order[i], data_directories)])
}

tyba_all <- NULL
for(i in seq_along(data_directories_sorted)) {
  cat("reading in assembly from", data_directories_sorted[i], "\n")
  setwd(data_directories_sorted[i])
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  repeats = read.csv(file = repeat_file)
  repeats$seqID = unlist(lapply(repeats$seqID, function(X) strsplit(X, split = "\n")[[1]][1]))
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  
  for(j in seq_along(custom_order)) if(grepl(custom_order[j], data_directories_sorted[i])) assembly <- custom_order[j]
  
  tyba$assembly <- assembly
  tyba$width <- tyba$end - tyba$start + 1
  tyba_all <- rbind(tyba_all, tyba[,c(1,8,14)])
}
tyba_all_b <- tyba_all

setwd("/home/pwlodzimierz/Rhynchospora/upload_data/histogram_all_tyba")

xlims = c(150,200)

tyba_all <- tyba_all[tyba_all$width < xlims[2],]
tyba_all <- tyba_all[tyba_all$width > xlims[1],]

hist(tyba_all$width, breaks = xlims[1]:xlims[2])

for(i in seq_along(unique(tyba_all$assembly))) {
  hist_data <- hist(tyba_all$width[tyba_all$assembly == unique(tyba_all$assembly)[i]], breaks = xlims[1]:xlims[2], plot = FALSE)
  if(i == 1) {
    cum_hist_data <- list(hist_data$counts)
  } else {
    cum_hist_data[[i]] = (cum_hist_data[[i-1]] + hist_data$counts)
  }
}

pdf("tyba_histogram_all.pdf")
for(i in length(cum_hist_data) : 1) {
  if(i == 20)  {
    plot(hist_data$mids, cum_hist_data[[i]], type = "h", col = colors_vector[i], lend = 1, lwd = 10)
  } else {
    points(hist_data$mids, cum_hist_data[[i]], type = "h", col = colors_vector[i], lend = 1, lwd = 10)
  }
}
dev.off()




















