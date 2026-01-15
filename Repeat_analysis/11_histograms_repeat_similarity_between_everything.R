#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

suppressMessages(library(seqinr))
suppressMessages(library(msa))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")



if(F) { # go to line 
  repeats <- read.csv("C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/temp_data/Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr_repeats_filtered.csv")

  genome <- "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr"
  
  i=1
}


data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
data_directories <- data_directories[!grepl(pattern = "Rtenuis.hap1.chr.fasta", data_directories)]
data_directories <- data_directories[!grepl(pattern = "Rtenuis.hap2.chr.fasta", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))



### make three tables: one that includes chromosomes from non-haplotype phased assemblies, so no between haplotypes peak will be there,
#                      one that includes all chromosomes from each species (one assembly per species)
#                      one for tenuis specifically, to add comparisons between same chromomes from different species

### Each point of data is from a single chromosome
### TODO: consider calculating these for every array: much more calculation, but possibility of discovering more detailed dynamics

max_space_between_repeats_to_ignore = 10000
sample_repeats_per_chromosome <- 50
sample_repeat_pairs_within_arrays <- 200
sample_repeat_from_whole_genomes <- 5000

## Table 1
# consider each chromosome from the genomes below, compare only within the genomes below for following stats:
# a. within arrays
# b. between arrays, same chromosome
# d. between arrays, different chromosome
# f. between arrays, different genomes


data_directories_unique_genomes <- data_directories[grepl(paste(c("R_barbata_hap1_chrs.fa",
                                                                  "R_holoschoenoides_hap1_chr.fasta",
                                                                  "R_pubera_ref_2n10_v2.chr.fasta",
                                                                  "Ralba.chr.fasta",
                                                                  "Rbreviuscula.hap1.chr.fasta", 
                                                                  "Rcephalotes.hap1.chr.fasta", 
                                                                  "Rcolorata.hap1.chr.fasta",
                                                                  "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                                  "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                                  "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta",
                                                                  "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                                  "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                                  "Rhync_riparia_5519C.hap1.chr.fasta",
                                                                  "Rhync_tenuis_ref.hap1.chr.fasta",
                                                                  "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta",
                                                                  "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta",
                                                                  "Rradicans.016.FINAL.chr.fasta",
                                                                  "Rrugosa.chr.fasta",
                                                                  "Rtenerrima.chr.fasta"), collapse = "|"), data_directories)]

all_repeats_sequences_to_sample <- NULL


for(i in 1 : length(data_directories_unique_genomes)) {
  cat(i, length(data_directories_unique_genomes), "\n")
  
  repeats <- read.csv(file = list.files(path = data_directories_unique_genomes[i], pattern = "repeats_filtered.csv", full.names = T)[1])
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  
  tyba <- tyba[tyba$width < 220, ] # to exclude the dimer Tybas
  
  genome <- strsplit(data_directories_unique_genomes[i], split = "TRASH_2025/")[[1]][2]
  
  all_repeats_sequences_to_sample <- rbind(all_repeats_sequences_to_sample, tyba[(sample(1 : nrow(tyba), sample_repeat_from_whole_genomes, replace = T)),c(8,11)])
  
}


table_1 <- data.frame(genome = vector(mode = "character"),
                      chromosome = vector(mode = "character"),
                      a.within_arrays = vector(mode = "numeric"),
                      b.between_arrays_same_chr = vector(mode = "numeric"),
                      c.between_arrays_diff_chr = vector(mode = "numeric"),
                      d.between_arrays_diff_genome = vector(mode = "numeric"))

# repeats_samples <- NULL


for(i in 1 : length(data_directories_unique_genomes)) {
  repeats <- read.csv(file = list.files(path = data_directories_unique_genomes[i], pattern = "repeats_filtered.csv", full.names = T)[1])
  
  genome <- strsplit(data_directories_unique_genomes[i], split = "TRASH_2025/")[[1]][2]
  
  
  chromosomes <- unique(repeats$seqID)
  
  tyba <- repeats[repeats$new_class %in% c("Tyba_172", "Tyba_182"), ]
  
  tyba <- tyba[tyba$width < 220, ] # to exclude the dimer Tybas
  
  cat(genome, "tyba:", nrow(tyba), i, "/", length(data_directories_unique_genomes), "\n")
  
  # for(j in seq_along(chromosomes)) {
  #   tyba_chr <- tyba[tyba$seqID == chromosomes[j], ]
  #   tyba_chr$genome <- genome
  #   if(nrow(tyba_chr) > sample_repeats_per_chromosome) {
  #     repeats_samples <- rbind(repeats_samples, tyba_chr[sample(1 : nrow(tyba_chr), sample_repeats_per_chromosome), c(1,11,14)])
  #   } else {
  #     repeats_samples <- rbind(repeats_samples, tyba_chr[, c(1,11,14)])
  #   }
  # }
  
  for(j in seq_along(chromosomes)) {
    table_1[nrow(table_1) + 1,] <- list(genome, chromosomes[j], 0, 0, 0, 0)
    
    
    tyba_chr <- tyba[tyba$seqID == chromosomes[j], ]
    tyba_chr$genome <- genome
    
    gaps <- tyba_chr$start[2 : nrow(tyba_chr)] - tyba_chr$end[1 : (nrow(tyba_chr) - 1)]
    
    gaps[gaps < max_space_between_repeats_to_ignore] = 0
    
    
    island_ends <- which(gaps > 0)
    island_starts <- c(1, (island_ends + 1))
    island_ends <- c(island_ends, nrow(tyba_chr))
    
    if(length(island_starts) == 1) {
      islands_chr <- data.frame(genome = rep(genome, length(island_starts)),
                                chromosome = rep(chromosomes[j], length(island_starts)),
                                start = tyba_chr$start[island_starts],
                                end = tyba_chr$end[island_ends],
                                start.no = island_starts,
                                end.no = island_ends,
                                distance_to_next = NA,
                                total_bp = tyba_chr$end[island_ends] - tyba_chr$start[island_starts] + 1,
                                tyba_no = island_ends - island_starts + 1)
    } else {
      islands_chr <- data.frame(genome = rep(genome, length(island_starts)),
                                chromosome = rep(chromosomes[j], length(island_starts)),
                                start = tyba_chr$start[island_starts],
                                end = tyba_chr$end[island_ends],
                                start.no = island_starts,
                                end.no = island_ends,
                                distance_to_next = c(tyba_chr$start[island_ends[2:length(island_ends)]] - tyba_chr$end[island_starts[1:(length(island_starts) - 1)]] - 1, NA),
                                total_bp = tyba_chr$end[island_ends] - tyba_chr$start[island_starts] + 1,
                                tyba_no = island_ends - island_starts + 1)
      
    }
    islands_chr <- islands_chr[islands_chr$tyba_no >= 3,]
    
    # a. within arrays
    islands_chr$sample_number <- round(sample_repeat_pairs_within_arrays * islands_chr$tyba_no/sum(islands_chr$tyba_no))
    internal_similarities <- NULL
    for(k in seq_len(nrow(islands_chr))) {
      cat(genome, "chromosome", j, "/", length(chromosomes), "array", k, "out of", nrow(islands_chr), "\n")
      sampleA <- sample(tyba_chr$sequence[islands_chr$start.no[k] : islands_chr$end.no[k]], size = islands_chr$sample_number[k], replace = TRUE)
      sampleB <- sample(tyba_chr$sequence[islands_chr$start.no[k] : islands_chr$end.no[k]], size = islands_chr$sample_number[k], replace = TRUE)
      distances <- adist(sampleA, sampleB)
      internal_similarities <- c(internal_similarities, 100 * (1 - distances / mean(nchar(c(sampleA, sampleB)))))
    }
    table_1$a.within_arrays[nrow(table_1)] <- mean(internal_similarities)
    
    # b. between arrays
    sampleA <- sample(tyba_chr$sequence, size = sample_repeats_per_chromosome, replace = TRUE)
    sampleB <- sample(tyba_chr$sequence, size = sample_repeats_per_chromosome, replace = TRUE)
    distances <- adist(sampleA, sampleB)
    
    table_1$b.between_arrays_same_chr[nrow(table_1)] <- 100 * (1 - mean(distances) / mean(nchar(c(sampleA, sampleB))))
    
    
    # c. between chromosomes, same species
    sampleA <- sample(tyba_chr$sequence, size = sample_repeats_per_chromosome, replace = TRUE)
    # sampleB <- sample(repeats_samples$sequence[repeats_samples$seqID != chromosomes[j]], size = sample_repeats_per_chromosome, replace = TRUE)
    sampleB <- sample(tyba$sequence[tyba$seqID != chromosomes[j]], size = sample_repeats_per_chromosome, replace = TRUE)
    distances <- adist(sampleA, sampleB)
    
    table_1$c.between_arrays_diff_chr[nrow(table_1)] <- 100 * (1 - mean(distances) / mean(nchar(c(sampleA, sampleB))))
    
    
    
    # d. between chromosomes, different species
    
    
    sampleA <- sample(tyba_chr$sequence, size = sample_repeats_per_chromosome, replace = TRUE)
    sampleB <- sample(all_repeats_sequences_to_sample$sequence, size = sample_repeats_per_chromosome, replace = TRUE)
    distances <- adist(sampleA, sampleB)
    
    table_1$d.between_arrays_diff_genome[nrow(table_1)] <- 100 * (1 - mean(distances) / mean(nchar(c(sampleA, sampleB))))
  }
  
}

write.csv(x = table_1, file = "/home/pwlodzimierz/Rhynchospora/upload_data/data_for_histograms_repeat_similarity_within_between_table_1.csv", row.names = F)
# table_1 <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/upload_data/data_for_histograms_repeat_similarity_within_between_table_1.csv")
# table_1 <- read.csv(file = "C:/Users/Piotr Włodzimierz/Desktop/Rhynchospora/upload_data/data_for_histograms_repeat_similarity_within_between_table_1.csv")

xmax = ceiling(max(c(table_1$a.within_arrays, table_1$b.between_arrays_same_chr, table_1$c.between_arrays_diff_chr, table_1$d.between_arrays_diff_genome))) + 5
# xmax = 40
ymax = 30
xmin = 40


pdf(file = "Histogrmas_table_1_similarity_tyba_within_outside.pdf", width = 12, height = 12)
par(mfrow = c(2,1))
hist(table_1$a.within_arrays, breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), border = "#ee000080", col = "#ee000080",
     xlab = "per chromosome mean pairwise distance of Tyba", main = "Histogram of Tyba pairwise distances averaged per chromosome")
hist(table_1$b.between_arrays_same_chr, breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#eeee0080", col = "#eeee0080")
hist(table_1$c.between_arrays_diff_chr, breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#ee00ee80", col = "#ee00ee80")
hist(table_1$d.between_arrays_diff_genome, breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#00ee0080", col = "#00ee0080")
legend(x = xmax * 0.9, y = ymax*0.95, legend = c("arrays", "chromosomes", "genomes", "species"), fill = c("#ee000080", "#eeee0080", "#ee00ee80","#00ee0080"))

y_max <- max(c(table_1$a.within_arrays, 
               table_1$b.between_arrays_same_chr, 
               table_1$c.between_arrays_diff_chr, 
               table_1$d.between_arrays_diff_genome), na.rm = TRUE)
offset <- 5

boxplot(table_1$a.within_arrays, table_1$b.between_arrays_same_chr, table_1$c.between_arrays_diff_chr, table_1$d.between_arrays_diff_genome,
        names = c("arrays", "chromosomes", "genomes", "species"),
        main = "per chromosome mean pairwise distance of Tyba repeats within:",
        col = c("#ee000080", "#eeee0080", "#ee00ee80","#00ee0080"),
        ylim = c(0, y_max + 3 * offset))


p1 <- t.test(table_1$a.within_arrays, table_1$b.between_arrays_same_chr, var.equal = FALSE)$p.value
p2 <- t.test(table_1$b.between_arrays_same_chr, table_1$c.between_arrays_diff_chr)$p.value
p3 <- t.test(table_1$c.between_arrays_diff_chr, table_1$d.between_arrays_diff_genome)$p.value


# Function to determine asterisk level
get_asterisks <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")  # not significant
}

# Draw lines and asterisks
segments(1, y_max + offset, 2, y_max + offset)
text(1.5, y_max + offset + 2, get_asterisks(p1))

segments(2, y_max + 2*offset, 3, y_max + 2*offset)
text(2.5, y_max + 2*offset + 2, get_asterisks(p2))

segments(3, y_max + 3*offset, 4, y_max + 3*offset)
text(3.5, y_max + 3*offset + 2, get_asterisks(p3))

dev.off()



## Table 2
# consider each chromosome from the genomes below, compare within the same genomes for following stats:
# a. within arrays
# b. between arrays, same chromosome
# d. between arrays, different chromosome
# f. between arrays, different genomes

# and compare with genomes of the other haplotype for following stats:
# c. between arrays, different haplotype, same chromosome
# e. between arrays, different haplotype, different chromosome



data_directories_unique_genomes_first_hap <- data_directories[grepl(paste(c("R_barbata_hap1_chrs.fa",
                                                                            "R_holoschoenoides_hap1_chr.fasta",
                                                                            "Rbreviuscula.hap1.chr.fasta", 
                                                                            "Rcephalotes.hap1.chr.fasta", 
                                                                            "Rcolorata.hap1.chr.fasta",
                                                                            "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
                                                                            "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta",
                                                                            "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta",
                                                                            "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                                            "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta",
                                                                            "Rhync_riparia_5519C.hap1.chr.fasta",
                                                                            "Rhync_tenuis_ref.hap1.chr.fasta",
                                                                            "Rrugosa.chr.fasta"), collapse = "|"), data_directories)]
data_directories_unique_genomes_second_hap <- data_directories[grepl(paste(c("R_barbata_hap2_chrs.fa",
                                                                             "R_holoschoenoides_hap2_chr.fasta",
                                                                             "Rbreviuscula.hap2.chr.fasta", 
                                                                             "Rcephalotes.hap2.chr.fasta", 
                                                                             "Rcolorata.hap2.chr.fasta",
                                                                             "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta",
                                                                             "Rhync_ciliata_6041A.hic.hap2.out_JBAT.FINAL.chr.fasta",
                                                                             "Rhync_corymbosa_6179D.asm.hic.FINAL.hap2.chr.fasta",
                                                                             "Rhync_nervosa_6321B.asm.hic.hap2.p_ctg.FINAL.chr.fasta",
                                                                             "Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr.fasta",
                                                                             "Rhync_riparia_5519C.hap2.chr.fasta",
                                                                             "Rhync_tenuis_ref.hap2.chr.fasta",
                                                                             "Rrugosa.chr.fasta"), collapse = "|"), data_directories)]

## Table 3
# compare tenuis genomes, using all chromosomes from haplotype 1 for:
# a. within arrays
# b. between arrays, same chromosome
# c. between arrays, different haplotype, same chromosome
# d. between arrays, different chromosome
# e. between arrays, different haplotype, different chromosome
# f. between arrays, different genome, same species, same chromosome, same haplotype
# g. between arrays, different genome, same species, same chromosome, different haplotype
# h. between arrays, different genome, same species, different chromosome, same haplotype
# i. between arrays, different genome, same species, different chromosome, different haplotype



data_directories_unique_genomes_first_hap <- data_directories[grepl(paste(c("Rhync_tenuis_6228D.036.hap1.chr.fasta",
                                                                            "Rhync_tenuis_6344A.JGV.hap1.chr.fasta",
                                                                            "Rhync_tenuis_6344B.PECP2.hap1.chr.fasta", 
                                                                            "Rhync_tenuis_6344C.PECP3.hap1.chr.fasta", 
                                                                            "Rhync_tenuis_6365A.035-5.hap1.chr.fasta",
                                                                            "Rhync_tenuis_6365B.036-7.hap1.chr.fasta",
                                                                            "Rhync_tenuis_6523A.hap1.chr.fasta",
                                                                            "Rhync_tenuis_6524A.hap1.chr.fasta",
                                                                            "Rhync_tenuis_ref.hap1.chr.fasta"), collapse = "|"), data_directories)]
data_directories_unique_genomes_second_hap <- data_directories[grepl(paste(c("Rhync_tenuis_6228D.036.hap2.chr.fasta",
                                                                             "Rhync_tenuis_6344A.JGV.hap2.chr.fasta",
                                                                             "Rhync_tenuis_6344B.PECP2.hap2.chr.fasta", 
                                                                             "Rhync_tenuis_6344C.PECP3.hap2.chr.fasta", 
                                                                             "Rhync_tenuis_6365A.035-5.hap2.chr.fasta",
                                                                             "Rhync_tenuis_6365B.036-7.hap2.chr.fasta",
                                                                             "Rhync_tenuis_6523A.hap2.chr.fasta",
                                                                             "Rhync_tenuis_6524A.hap2.chr.fasta",
                                                                             "Rhync_tenuis_ref.hap2.chr.fasta"), collapse = "|"), data_directories)]




























































