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
genomeA = c(rep("Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta", 1),
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


genomeB = c(rep("Rhync_ciliata_6041A.hic.hap2.out_JBAT.FINAL.chr.fasta", 1),
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

chromosomeA = rep("", length(genomeA))
chromosomeB = rep("", length(genomeA))


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

write.csv(x = table, file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv", row.names = F)



































