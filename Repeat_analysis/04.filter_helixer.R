#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))
suppressMessages(library(msa))
suppressMessages(library(GenomicRanges))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

replace_existing <- FALSE

helixer_overlap_max_perc = 0.8 # if an helixer annotation overlaps more than this with repeat ARRAY coordinates, it is discarded

data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
setwd(data_directories[i])
print(getwd())

helixer_file <- list.files(pattern = "helixer.gff3", full.names = TRUE)
print(helixer_file)


if(length(helixer_file) != 1) {warning(paste0(i, "No genes!")); setwd(".."); quit(save = "no", status = 1)}

repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {print(paste0(i, "No repeats!")); setwd(".."); next}

array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
if(length(array_file) != 1) {print(paste0(i, " no arrays!")); setwd(".."); next}

classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
if(length(classes_file) != 1) {print(paste0(i, " no classes!")); setwd(".."); next}

edta_file = list.files(pattern = "_edta_modified.csv", full.names = TRUE)
if(length(edta_file) != 1) {print(paste0(i, " no edta!")); setwd(".."); next}

assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
assembly_file = grep(strsplit(strsplit(repeat_file, split = "_repeats")[[1]][1], split = "/")[[1]][2], assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])

if(!replace_existing) {
  if(file.exists(paste0(assembly_name, "_helixer_filtered.csv"))) {
    print("Already done, not repeating")
    quit(save = "no", status = 0)
  }
}

repeats = read.csv(file = repeat_file)
arrays = read.csv(file = array_file)
classes = read.csv(file = classes_file)
classes$num_ID <- 1 : nrow(classes)
edta = read.csv(file = edta_file)
helixer = read.table(file = helixer_file, header = FALSE, sep = "\t", skip = 4)


helixer$overlapping_bp = 0
helixer_filtered_total = NULL

print("filter helixer")

all_chromosomes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
all_chromosomes <- all_chromosomes[all_chromosomes$assembly.name == strsplit(data_directories[i], split = "TRASH_2025/")[[1]][2], ]
chromosomes <- all_chromosomes$chromosome.name
chromosomes_lengths <- all_chromosomes$size

for (j in seq_along(chromosomes)) {  
  print(j)
  sequence_arrays = arrays[arrays$seqID == chromosomes[j], ]
  sequence_helixer <- helixer[helixer$V1 == chromosomes[j],]
  
  if(nrow(sequence_arrays) == 0 | nrow(sequence_helixer) == 0) next
  
  
  helixer_IDs <- unique(sequence_helixer$V9[sequence_helixer$V3 == "gene"])
  helixer_cluster_start <- which(sequence_helixer$V3 == "gene")
  helixer_cluster_end = c(helixer_cluster_start[-1], nrow(sequence_helixer))
  
  for(k in seq_along(helixer_IDs)) {
    helixer_IDs[k] = strsplit(helixer_IDs[k], split = "ID=")[[1]][2]
  }
  
  
  which_sequence_CDS <- which(sequence_helixer$V3 == "CDS")
  sequence_helixer_CDS = sequence_helixer[which_sequence_CDS, ]
  sequence_helixer_CDS$overlapping_bp <- 0
  
  gr1 <- with(sequence_arrays, GRanges(chromosomes[j], IRanges(start, end)))
  gr2 <- with(sequence_helixer_CDS, GRanges(chromosomes[j], IRanges(V4, V5)))
  
  overlaps <- as.data.frame(findOverlaps(gr1, gr2))
  
  for(k in seq_len(nrow(overlaps))) {
    print(paste(chromosomes[j], k, "/", nrow(overlaps)))
    overlap_bp = width(pintersect(gr1[overlaps$queryHits[k]], gr2[overlaps$subjectHits[k]]))
    sequence_helixer_CDS$overlapping_bp[overlaps$subjectHits[k]] = sequence_helixer_CDS$overlapping_bp[overlaps$subjectHits[k]] + overlap_bp
  }
  sequence_helixer$overlapping_bp[which_sequence_CDS] <- sequence_helixer_CDS$overlapping_bp
  
  sequence_helixer$width = sequence_helixer$V5 - sequence_helixer$V4 + 1
  
  genes_to_remove <- NULL
  
  for(k in seq_along(helixer_IDs)) { 
    gene_cluster <- sequence_helixer[helixer_cluster_start[k] : helixer_cluster_end[k], ]
    gene_CDS_total_size <- sum(gene_cluster$width[gene_cluster$V3 == "CDS"])
    gene_CDS_total_overlap <- sum(gene_cluster$overlapping_bp)
    if(gene_CDS_total_overlap == 0) next
    if((gene_CDS_total_overlap / gene_CDS_total_size) > helixer_overlap_max_perc) {
      genes_to_remove <- c(genes_to_remove, (helixer_cluster_start[k] : helixer_cluster_end[k]))
    }
  }
  if(!is.null(genes_to_remove)) {
    helixer_filtered <- sequence_helixer[-genes_to_remove, ]
  } else {
    helixer_filtered <- sequence_helixer
  }
  
  
  helixer_filtered_total = rbind(helixer_filtered_total, helixer_filtered)
  
  
}
nrow(helixer_filtered_total)
nrow(helixer)

helixer <- helixer_filtered_total
export_gff(annotations.data.frame = helixer, output = ".", 
           file.name = paste0(assembly_name, "_helixer_filtered"), seqid = 1, source = 2, type = 3, 
           start = 4, end = 5, strand = 7, attributes = 9, attribute.names = "")

write.csv(helixer, paste0(assembly_name, "_helixer_filtered.csv"))







