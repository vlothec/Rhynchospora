setwd("/home/pwlodzimierz/Rhynchospora/synteny_analysis")

repeats.files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", 
                            recursive = TRUE, full.names = TRUE, pattern = "repeats_with_seq.csv")

repeats.files <- sort(repeats.files, decreasing = TRUE)

array.files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", 
                          recursive = TRUE, full.names = TRUE, pattern = "islands_genome")

array.files <- sort(array.files, decreasing = TRUE)
# 
# assemblies_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", 
#                                pattern = ".fa", full.names = TRUE, recursive = TRUE)
# 
# assemblies_files <- sort(assemblies_files, decreasing = TRUE)
#
if(F) {
  genomeA = c(rep("Rrugosa.chr.fasta", 18),
              rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 5),
              rep("Rcephalotes.hap1.chr.fasta", 9),
              rep("Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta", 5),
              rep("Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 5),
              rep("Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 9),
              rep("R_barbata_hap1_chrs.fa", 5),
              rep("Rcephalotes.hap1.chr.fasta", 9),
              rep("Rhync_riparia_5519C.hap1.chr.fasta", 5),
              rep("R_barbata_hap1_chrs.fa", 5),
              rep("R_holoschoenoides_hap1_chr.fasta", 5),
              rep("Rbreviuscula.hap1.chr.fasta", 5),
              rep("Rhync_tenuis_6523A.hap1.chr.fasta", 2),
              rep("Rhync_tenuis_6523A.hap2.chr.fasta", 2),
              rep("Rhync_tenuis_6524A.hap1.chr.fasta", 2),
              rep("Rhync_tenuis_6524A.hap2.chr.fasta", 2))
  
  
  genomeB = c(rep("Ralba.chr.fasta", 18),
              rep("Rhync_ciliata_6041A.hic.hap2.out_JBAT.FINAL.chr.fasta", 5),
              rep("Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr.fasta", 9),
              rep("Rradicans.016.FINAL.chr.fasta", 5),
              rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 5),
              rep("Rhync_tenuis_ref.hap1.chr.fasta", 9),
              rep("R_holoschoenoides_hap1_chr.fasta", 5),
              rep("Rcephalotes.hap2.chr.fasta", 9),
              rep("Rhync_riparia_5519C.hap2.chr.fasta", 5),
              rep("R_barbata_hap2_chrs.fa", 5),
              rep("R_holoschoenoides_hap2_chr.fasta", 5),
              rep("Rbreviuscula.hap2.chr.fasta", 5),
              rep("Rhync_tenuis_6524A.hap1.chr.fasta", 2),
              rep("Rhync_tenuis_6524A.hap2.chr.fasta", 2),
              rep("Rhync_tenuis_6524A.hap2.chr.fasta", 2),
              rep("Rhync_tenuis_6524A.hap1.chr.fasta", 2))
  
  
  chromosomeA = c(paste0("chr", 1:18),
                  paste0("Chr", c(1,2,3,4,5), "_h1"),
                  paste0("chr", 1:9, "_h1"),
                  paste0("HiC_scaffold_", 1:5),
                  paste0("HiC_scaffold_", 1:5, "_h1"),
                  "scaffold_1_a1", "scaffold_7_a1", "scaffold_13_a1",   "scaffold_3_b1", "scaffold_11_b1", "scaffold_15_b1",   "scaffold_5_a1", "scaffold_8_a1", "scaffold_17_a1",
                  paste0("chr", 1:5, "_RagTag"),
                  paste0("chr", 1:9, "_h1"),
                  paste0("chr", 1:5, "_h1"),
                  paste0("chr", 1:5, "_RagTag"),
                  paste0("Chr", 1:5, "_RagTag_h1"),
                  paste0("Chr", 1:5, "_h1"),
                  paste0("Chr", 1:2, "_h1"),
                  paste0("Chr", 1:2, "_h2"),
                  paste0("Chr", 1:2, "_h1"),
                  paste0("Chr", 1:2, "_h2"))
  
  
  chromosomeB = c(paste0("Chr", c(8,13,5,5,10,11,3,9,9,3,12,1,6,2,4,1,4,6)),
                  paste0("Chr", c(1,2,3,4,5), "_h2"),
                  paste0("HiC_scaffold_", c(6,5,3,2,5,2,1,1,1), "_h1"),
                  paste0("Chr", c(1,2,3,4,5)),
                  paste0("Chr", c(1,2,3,4,5), "_h1"),
                  paste0("Chr", c(1,1,1,2,2,2,2,2,1), "_h1"),
                  paste0("Chr", c(5,6,4,1,9), "_RagTag_h1"),
                  paste0("chr", 1:9, "_h2"),
                  paste0("chr", 1:5, "_h2"),
                  paste0("chr", 1:5, "_RagTag"),
                  paste0("Chr", 1:5, "_RagTag_h2"),
                  paste0("Chr", 1:5, "_h1"),
                  paste0("Chr", 1:2, "_h1"),
                  paste0("Chr", 1:2, "_h2"),
                  paste0("Chr", 1:2, "_h2"),
                  paste0("Chr", 1:2, "_h1"))
}

genomeAt = c(rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 1),
             rep("Rcephalotes.hap1.chr.fasta", 1),
             rep("Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 1),
             rep("R_barbata_hap1_chrs.fa", 1),
             rep("Rhync_riparia_5519C.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             rep("R_holoschoenoides_hap1_chr.fasta", 1),
             rep("Rhync_riparia_5519C.hap1.chr.fasta", 1),
             rep("Rbreviuscula.hap1.chr.fasta", 1),
             rep("Rcolorata.hap1.chr.fasta", 1),
             rep("Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             # down here pairs of assemblies
             rep("Rrugosa.chr.fasta", 1),
             rep("Ralba.chr.fasta", 1),
             rep("Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta", 1),
             rep("Rcephalotes.hap1.chr.fasta", 1),
             rep("Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta", 1),
             rep("Rradicans.016.FINAL.chr.fasta", 1),
             rep("R_pubera_ref_2n10_v2.chr.fasta", 1),
             rep("Rbreviuscula.hap1.chr.fasta", 1),
             rep("Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 1),
             rep("Rcolorata.hap1.chr.fasta", 1),
             rep("Rtenerrima.chr.fasta", 1),
             rep("Rhync_filiformis_5519G.hap1.chr.fasta", 1),
             rep("Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             rep("Rhync_riparia_5519C.hap1.chr.fasta", 1),
             rep("R_barbata_hap1_chrs.fa", 1),
             rep("R_holoschoenoides_hap1_chr.fasta", 1))


genomeBt = c(rep("Rhync_ciliata_6041A.hic.hap2.out_JBAT.FINAL.chr.fasta", 1),
             rep("Rcephalotes.hap2.chr.fasta", 1),
             rep("Rhync_ridleyi.asm.hic.hap2.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_nervosa_6321B.asm.hic.hap2.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta", 1),
             rep("R_barbata_hap2_chrs.fa", 1),
             rep("Rhync_riparia_5519C.hap2.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap2.chr.fasta", 1),
             rep("R_holoschoenoides_hap2_chr.fasta", 1),
             rep("Rhync_riparia_5519C.hap2.chr.fasta", 1),
             rep("Rbreviuscula.hap2.chr.fasta", 1),
             rep("Rcolorata.hap2.chr.fasta", 1),
             rep("Rhync_corymbosa_6179D.asm.hic.FINAL.hap2.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap2.chr.fasta", 1),
             # down here pairs of assemblies
             rep("Ralba.chr.fasta", 1),
             rep("Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta", 1),
             rep("Rcephalotes.hap1.chr.fasta", 1),
             rep("Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta", 1),
             rep("Rradicans.016.FINAL.chr.fasta", 1),
             rep("R_pubera_ref_2n10_v2.chr.fasta", 1),
             rep("Rbreviuscula.hap1.chr.fasta", 1),
             rep("Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta", 1),
             rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 1),
             rep("Rcolorata.hap1.chr.fasta", 1),
             rep("Rtenerrima.chr.fasta", 1),
             rep("Rhync_filiformis_5519G.hap1.chr.fasta", 1),
             rep("Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             rep("Rhync_riparia_5519C.hap1.chr.fasta", 1),
             rep("R_barbata_hap1_chrs.fa", 1),
             rep("R_holoschoenoides_hap1_chr.fasta", 1),
             rep("Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta", 1))



chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")

genomeA = NULL
genomeB = NULL
chromosomeA = NULL
chromosomeB = NULL


for(i in 1 : 13) {
  chrs_A <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genomeAt[i]]
  chrs_B <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genomeBt[i]]
  if(length(chrs_A) != length(chrs_B)) cat("check", i, "\n")
  genomeA <- c(genomeA, rep(genomeAt[i], length(chrs_A)))
  genomeB <- c(genomeB, rep(genomeBt[i], length(chrs_A)))
  chromosomeA <- c(chromosomeA, chrs_A)
  chromosomeB <- c(chromosomeB, chrs_B)
}
for(i in 14 : 31) {
  chrs_A <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genomeAt[i]]
  chrs_B <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genomeBt[i]]
  genomeA <- c(genomeA, rep(genomeAt[i], length(chrs_A) * length(chrs_B)))
  genomeB <- c(genomeB, rep(genomeBt[i], length(chrs_A) * length(chrs_B)))
  for(j in seq_along(chrs_A)) {
    for(k in seq_along(chrs_B)) {
      chromosomeA <- c(chromosomeA, chrs_A[j])
      chromosomeB <- c(chromosomeB, chrs_B[k])
    }
  }
}
# run only 1 to 79 to run haplotypes same chromosome


# tenuis pairs chr 1
genomeA <- c(genomeA, 
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344A.JGV.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6523A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6524A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344B.PECP2.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344C.PECP3.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6228D.036.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365A.035-5.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_ref.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344A.JGV.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6523A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6524A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344B.PECP2.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344C.PECP3.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6228D.036.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365A.035-5.hap1.chr.fasta", 1)
             )
genomeB <- c(genomeB,
             rep("Rhync_tenuis_6344A.JGV.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6523A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6524A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344B.PECP2.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344C.PECP3.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6228D.036.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365A.035-5.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365B.036-7.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344A.JGV.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6523A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6524A.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344B.PECP2.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6344C.PECP3.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6228D.036.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365A.035-5.hap1.chr.fasta", 1),
             rep("Rhync_tenuis_6365B.036-7.hap1.chr.fasta", 1)
             )


chromosomeA <- c(chromosomeA,
                 'Chr1_h2',
                 'JGV_Chr1_h1_h1',
                 'Chr1_h1',
                 'Chr1_h1',
                 'PECP2_Chr1_h1',
                 'PECP3_Chr1_h1',
                 '036_Chr1_h1',
                 '035-6_Chr1_h1',
                 'Chr2_h1',
                 'JGV_Chr1_h2_h1',
                 'Chr2_h1',
                 'Chr2_h1',
                 'PECP2_Chr2_h1',
                 'PECP3_Chr2_h1',
                 '036_Chr2_h1',
                 '035-6_Chr2_h1')

chromosomeB <- c(chromosomeB,
                 'JGV_Chr1_h1_h1',
                 'Chr1_h1',
                 'Chr1_h1',
                 'PECP2_Chr1_h1',
                 'PECP3_Chr1_h1',
                 '036_Chr1_h1',
                 '035-6_Chr1_h1',
                 '036-7_Chr1_h1',
                 'JGV_Chr1_h2_h1',
                 'Chr2_h1',
                 'Chr2_h1',
                 'PECP2_Chr2_h1',
                 'PECP3_Chr2_h1',
                 '036_Chr2_h1',
                 '035-6_Chr2_h1',
                 '036-7_Chr2_h1')


# autrobrasiliensis pairs
temp_gen <- c( 
"Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
"Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
"Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta",
"Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta",
"Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta",
"Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta"
)
temp_chr <- c(
  "scaffold_1_a1",
  "scaffold_3_b1",
  "scaffold_5_c1",
  "scaffold_2_a2",
  "scaffold_4_b2",
  "scaffold_6_c2"
)
for(k in 1 : (length(temp_gen)-1)) {
  for(l in (k + 1) : length(temp_gen)) {
    genomeA <- c(genomeA, temp_gen[k])
    genomeB <- c(genomeB, temp_gen[l])
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

temp_chr <- c(
  "scaffold_7_a1",
  "scaffold_8_c1",
  "scaffold_11_b1",
  "scaffold_9_c2",
  "scaffold_10_a2",
  "scaffold_12_b2"
)
for(k in 1 : (length(temp_gen)-1)) {
  for(l in (k + 1) : length(temp_gen)) {
    genomeA <- c(genomeA, temp_gen[k])
    genomeB <- c(genomeB, temp_gen[l])
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

temp_chr <- c(
  "scaffold_13_a1",
  "scaffold_15_b1",
  "scaffold_17_c1",
  "scaffold_16_b2",
  "scaffold_14_a2",
  "scaffold_18_c2"
)
for(k in 1 : (length(temp_gen)-1)) {
  for(l in (k + 1) : length(temp_gen)) {
    genomeA <- c(genomeA, temp_gen[k])
    genomeB <- c(genomeB, temp_gen[l])
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

# gaudichaudii pairs
genomeA <- c(genomeA, 
  rep("Rhync_gaudichaudii.chr.h1.fasta", 5),
  rep("Rhync_gaudichaudii.chr.h2.fasta", 5),
  rep("Rhync_gaudichaudii.chr.h3.fasta", 5)
)
genomeB <- c(genomeB, 
  rep("Rhync_gaudichaudii.chr.h2.fasta", 5),
  rep("Rhync_gaudichaudii.chr.h3.fasta", 5),
  rep("Rhync_gaudichaudii.chr.h1.fasta", 5)
)
chromosomeA <- c(chromosomeA,
  paste0("Chr0", 1:5, "_h1"),
  paste0("Chr0", 1:5, "_h2"),
  paste0("Chr0", 1:5, "_h3")
)
chromosomeB <- c(chromosomeB,
  paste0("Chr0", 1:5, "_h2"),
  paste0("Chr0", 1:5, "_h3"),
  paste0("Chr0", 1:5, "_h1")
)

# tenerrima pairs
temp_gen <- 'Rtenerrima.chr.fasta'
temp_chr <- paste0("chr", 1:10)
for(k in 1 : (length(temp_chr) - 1)) {
  for(l in (k + 1) : length(temp_chr)) {
    genomeA <- c(genomeA, temp_gen)
    genomeB <- c(genomeB, temp_gen)
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

# pubera pairs
temp_gen <- 'R_pubera_ref_2n10_v2.chr.fasta'
temp_chr <- paste0("Chr0", 1:5)
for(k in 1 : (length(temp_chr) - 1)) {
  for(l in (k + 1) : length(temp_chr)) {
    genomeA <- c(genomeA, temp_gen)
    genomeB <- c(genomeB, temp_gen)
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

# ridleyi pairs
temp_gen <- 'Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta'
temp_chr <- paste0("HiC_scaffold_", 1:6, "_h1")
for(k in 1 : (length(temp_chr) - 1)) {
  for(l in (k + 1) : length(temp_chr)) {
    genomeA <- c(genomeA, temp_gen)
    genomeB <- c(genomeB, temp_gen)
    chromosomeA <- c(chromosomeA, temp_chr[k])
    chromosomeB <- c(chromosomeB, temp_chr[l])
  }
}

table <- data.frame(no = 1 : length(genomeA),
                    genome_A = genomeA,
                    chromosome_A = chromosomeA, 
                    genome_B = genomeB,
                    chromosome_B = chromosomeB,
                    repeats_file_A = genomeA,
                    arrays_file_A = genomeA,
                    repeats_file_B = genomeB,
                    arrays_file_B = genomeB)


remove_rows <- NULL
for(i in 1 : nrow(table)) {
  
  start <- "/home/pwlodzimierz/Rhynchospora/TRASH_2025/"
  gen_a_sh <- strsplit(table$genome_A[i], split = ".fa")[[1]][1]
  gen_b_sh <- strsplit(table$genome_B[i], split = ".fa")[[1]][1]
  
  array_file_A <- paste0(start, table$genome_A[i], "/", gen_a_sh, "_islands_genome.csv")
  array_file_A_alt <- paste0(start, table$genome_A[i], "/", gen_a_sh, "_islands_genome_no_edta.csv")
  array_file_B <- paste0(start, table$genome_B[i], "/", gen_b_sh, "_islands_genome.csv")
  array_file_B_alt <- paste0(start, table$genome_B[i], "/", gen_b_sh, "_islands_genome_no_edta.csv")
  
  if(file.exists(array_file_A)) {
    table$arrays_file_A[i] = array_file_A
  } else if(file.exists(array_file_A_alt)) {
    table$arrays_file_A[i] = array_file_A_alt
  } else {
    remove_rows <- c(remove_rows, i)
    next
  }
  if(file.exists(array_file_B)) {
    table$arrays_file_B[i] = array_file_B
  } else if(file.exists(array_file_B_alt)) {
    table$arrays_file_B[i] = array_file_B_alt
  } else {
    remove_rows <- c(remove_rows, i)
    next
  }
  
  
  repeats_file_A <- paste0(start, table$genome_A[i], "/", gen_a_sh, "_repeats_filtered.csv")
  repeats_file_B <- paste0(start, table$genome_B[i], "/", gen_b_sh, "_repeats_filtered.csv")
  
  if(file.exists(repeats_file_A)) {
    table$repeats_file_A[i] = repeats_file_A
  } else {
    remove_rows <- c(remove_rows, i)
    next
  }
  if(file.exists(repeats_file_B)) {
    table$repeats_file_B[i] = repeats_file_B
  } else {
    remove_rows <- c(remove_rows, i)
    next
  }
  
  
}

# table <- table[-remove_rows, ]

write.csv(x = table, file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table.csv", row.names = F)



































