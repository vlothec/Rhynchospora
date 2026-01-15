#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

suppressMessages(library(seqinr))
suppressMessages(library(msa))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

replace_existing_analysis = FALSE
plot = FALSE

data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

templates = read.fasta(file = "/home/pwlodzimierz/Rhynchospora/template1.fasta")
templates = c(templates, read.fasta(file = "/home/pwlodzimierz/Rhynchospora/template2.fasta"))
templates = c(templates, read.fasta(file = "/home/pwlodzimierz/Rhynchospora/template3.fasta"))
templates = c(templates, read.fasta(file = "/home/pwlodzimierz/Rhynchospora/templateholos.fasta"))
templates = c(templates, read.fasta(file = "/home/pwlodzimierz/Rhynchospora/templatecorymbosa.fasta"))
templates = c(templates, read.fasta(file = "/home/pwlodzimierz/Rhynchospora/templatebarbata.fasta"))

#TODO: instead of reading the whole assembly, read in summary file with assembly sequences names and sizes

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 30
print(i)

# for(i in seq_along(data_directories)) 
{
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  print(data_directories[i])
  
  repeat_file = list.files(pattern = "_repeats_with_seq.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, " no repeats!")); next}
  
  array_file = list.files(pattern = "_arrays.csv")
  array_file = array_file[!grepl(pattern = "no_repeats_", array_file)]
  if(length(array_file) != 1) {print(paste0(i, " no arrays!")); next}
  
  assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "TRASH_2025/")[[1]][2]
  assembly_file = grep(strsplit(strsplit(repeat_file, split = "_repeats")[[1]][1], split = "/")[[1]][2], assembly_files)
  print(assembly_file)
  print(assembly_files[assembly_file])
  
  if(!replace_existing_analysis) {
    if(file.exists(paste0(assembly_name, "_repeats_filtered.csv"))) {
      cat("file exists, not overwriting and exiting")
      quit(save = "no", status = 1)
      }
  }
  
  if(length(assembly_file) == 1) {
    assembly = read.fasta(assembly_files[assembly_file], seqtype = "DNA", forceDNAtolower = TRUE)
  } else {
    print(" No assembly found, cannot proceed")
    setwd(def_wd)
    next
  }

  repeats = read.csv(file = repeat_file)
  
  arrays = read.csv(file = array_file)
  
  arrays$arrayID = seq_len(nrow(arrays))
  
  
  if(length(grep("\n", repeats$seqID)) > 0) {
    repeats$seqID[grep("\n", repeats$seqID)] <- unlist(lapply(grep("\n", repeats$seqID), function(X) strsplit(repeats$seqID[X], split = "\n")[[1]][1]))
  }
  if(length(grep("\n", arrays$seqID)) > 0) {
    arrays$seqID[grep("\n", arrays$seqID)] <- unlist(lapply(grep("\n", arrays$seqID), function(X) strsplit(arrays$seqID[X], split = "\n")[[1]][1]))
  }
  
  
  
  cat("assembly chromosme names:", unique(names(assembly)), "\n")
  cat("repeats unique chromosome names:", unique(repeats$seqID), "\n")
  
  repeats$arrayID = 0
  for(j in 1 : nrow(arrays)) {
    print(paste0("B: ", i, " ", j, " / ", nrow(arrays)))
    repeats$arrayID[repeats$seqID == arrays$seqID[j] & 
                      ((repeats$start >= arrays$start[j] & 
                          repeats$start <= arrays$end[j]) | 
                         (repeats$end >= arrays$start[j] & 
                            repeats$end <= arrays$end[j])) ] = j
  }
  
  ### Remove repeats (with no class) with score over 60 =====
  nrow(repeats) # 139458
  maxed = 60
  repeats = repeats[repeats$score <= maxed, ]
  
  ### remove arrays and repeats with less than 100 bp of repeats =========
  min_rep_bp_array <- 100
  nrow(repeats) # 136956
  nrow(arrays)  # 17937
  j = 0
  while (j < nrow(arrays)) {
    j = j + 1
    cat("C: ", i, " ", j, " / ", nrow(arrays), "\n", sep = "")
    if(arrays$class[j] %in% names(templates)) next # unless it had a template
    which_repeats = which(repeats$arrayID == arrays$arrayID[j])
    if(length(which_repeats) == 0) {
      arrays = arrays[-j, ]
      j = j - 1
    } else if(sum(repeats$width[which_repeats]) <= min_rep_bp_array) {
      repeats = repeats[-which_repeats, ]
      arrays = arrays[-j, ]
      j = j - 1
    } else {
      arrays$repeats_number[j] = length(which_repeats)
    }
  }
  ### remove arrays and repeats with less than 5 repeats, when average adist score is higher than 30
  nrow(repeats) # 119060
  nrow(arrays)  # 1321
  j = 0
  minr = 5
  maxed = 30
  while (j < nrow(arrays)) {
    j = j + 1
    cat("D: ", i, " ", j, " / ", nrow(arrays), "\n", sep = "")
    if(arrays$class[j] %in% names(templates)) next # unless it had a template
    which_repeats = which(repeats$arrayID == arrays$arrayID[j])
    if(length(which_repeats) < minr) {
      if(mean(repeats$score[which_repeats]) > maxed) {
        arrays = arrays[-j, ]
        repeats = repeats[-which_repeats, ]
        j = j - 1
      }
    } 
  }
  nrow(repeats) # 115285
  nrow(arrays)  # 889
  
  
  ### Summarise classes ====================================
  classes = data.frame(class = unique(repeats$class))
  
  classes$count = 0
  classes$consensus = ""
  classes$median_length = 0
  classes$length_SD = 0
  classes$mean_edit_score = 0 
  classes$total_bp = 0
  
  min_repeats_to_align = 10
  for(j in seq_len(nrow(classes))) {
    repeats_class = repeats[repeats$class == classes$class[j], ]
    
    cat("E:", i, j, "/", nrow(classes), "\n", sep = " ")
    repeats_class$chr = unlist(lapply(repeats_class$seqID, function(X) (seq_along(unique(names(assembly))))[unique(names(assembly)) %in% X]))
    
    classes$count[j] = nrow(repeats_class)
    classes$median_length[j] = round(median(repeats_class$width))
    classes$length_SD[j] = sd(repeats_class$width)
    classes$total_bp[j] = sum(repeats_class$width)
    
    class_name_nchar = round(median(repeats_class$width))
    
    if(nrow(repeats_class) < min_repeats_to_align) {
      repeats_to_align = seq_len(nrow(repeats_class))
    } else {
      repeats_to_align = sample(seq_len(nrow(repeats_class)), (nrow(repeats_class) %/% 100 + min_repeats_to_align), replace = FALSE)
    }
    # sequences_to_align <- unlist(lapply(repeats_to_align, function(X) paste0(assembly[[repeats_class$chr[X]]][repeats_class$start[X] : repeats_class$end[X]], collapse = "")))
    # sequences_to_align[repeats_class$strand[repeats_to_align] == "-"] <- unlist(lapply(sequences_to_align[repeats_class$strand[repeats_to_align] == "-"], revCompString))
    sequences_to_align <- repeats_class$sequence[repeats_to_align]
    if(length(sequences_to_align) == 0) {
      cat("\n\n\n\n\n\n\n\n\n\n")
      stop(paste0(assembly_file, " did not find repeats in one of the classes: ", classes$class[j], ", investigate"))
    } else if(length(sequences_to_align) == 1) {
      classes$consensus[j] <- sequences_to_align
      classes$mean_edit_score[j] = 0
    } else {
      a <- capture.output({alignment_matrix = tolower(as.matrix(msa(sequences_to_align, method = "ClustalOmega", type = "dna")))})
      classes$consensus[j] <- consensus_N(alignment_matrix, class_name_nchar)
      classes$mean_edit_score[j] = mean(adist(classes$consensus[j], sequences_to_align))
    }
  }
  ### save classes
  classes$importance = classes$count * classes$median_length

  write.csv(classes, paste0(assembly_name, "_classes_filtered.csv"), row.names = FALSE)

  classes$new_class = ""
  classes$new_importance = 0
  
  ### remove classes with less than 2000 bp of repeats for merging, to speed it up ========
  if(nrow(classes[classes$total_bp >= 4000, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 4000, ]
    classes = classes[classes$total_bp >= 4000, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else if(nrow(classes[classes$total_bp >= 2000, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 2000, ]
    classes = classes[classes$total_bp >= 2000, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else if(nrow(classes[classes$total_bp >= 1000, ]) >= 20) {
    classes_discarded = classes[classes$total_bp < 1000, ]
    classes = classes[classes$total_bp >= 1000, ]
    if(nrow(classes_discarded) > 0) {
      classes_discarded$new_class = classes_discarded$class
      classes_discarded$score = 0
      classes_discarded$sum_coverage = classes_discarded$importance
    }
  } else {
    classes_discarded = classes[NULL, ]
  }
  
  ### check kmer based similarity between all ==============
  

  classes = classes[order(classes$importance, decreasing = TRUE), ]
  
  similarity_matrix = matrix(nrow = nrow(classes), ncol = nrow(classes))
  for(j in 1 : nrow(similarity_matrix)) {
    print(paste0("E2: ", i, " ", j, " / ", nrow(similarity_matrix)))
    for(k in j : ncol(similarity_matrix)) {
      similarity_matrix[j,k] = kmer_compare(classes$consensus[j], classes$consensus[k])
      similarity_matrix[k,j] = kmer_compare(classes$consensus[k], classes$consensus[j])
    }
  }
  

  for(j in 1 : nrow(classes)) {
    print(paste0("F: ", i, " ", j, " / ", nrow(classes)))
    if(classes$new_class[j] == "") {
      classes$new_class[j] = classes$class[j]
      other_similar = which(similarity_matrix[j, ] > 0.12 & similarity_matrix[, j] > 0.12)
      classes$new_class[other_similar] = classes$class[j]
      classes$new_importance[j] = sum(classes$importance[other_similar])
    }
  }
  
  # unique(classes$new_class)
  classes_new = classes[classes$new_importance != 0, ]
  classes_new$importance = classes_new$new_importance
  classes_new$class = classes_new$new_class
  classes_new$score = classes_new$mean_edit_score / classes_new$median_length
  classes_new$sum_coverage = classes_new$importance
  
  if(nrow(classes_discarded) != 0) {
    classes_new = rbind(classes_new, classes_discarded)
  } 
  
  
  write.csv(classes_new[c("class", "count", "consensus", "median_length", "length_SD", "score", "sum_coverage")], 
            paste0(assembly_name, "_classes_merged_filtered.csv"), row.names = FALSE)
  
  
  ### Assign new classes to repeats and arrays =============
  
  classes_new$class_num_ID = 1 : nrow(classes_new)
  classes$class_num_ID = 0
  for(j in 1 : nrow(classes)) {
    classes$class_num_ID[j] = classes_new$class_num_ID[which(classes_new$new_class == classes$new_class[j])]
  }
  
  repeats$new_class = ""
  repeats$new_class_num_ID = 0
  arrays$new_class_num_ID = 0
  for(j in 1 : nrow(classes)){
    print(paste0("G: ", i, " ", j, " / ", nrow(classes)))
    repeats$new_class[repeats$class == classes$class[j]] = classes$new_class[j]
    repeats$new_class_num_ID[repeats$class == classes$class[j]] = classes$class_num_ID[j]
    arrays$new_class_num_ID[arrays$class == classes$class[j]] = classes$class_num_ID[j]
  }
  repeats$new_class[repeats$new_class == ""] = repeats$class[repeats$new_class == ""]
  
  ### Shrink arrays to only reach over repeats =============
  
  
  for(j in seq_len(nrow(arrays))) {
    print(paste0("H: ", i, " ", j, " / ", nrow(arrays)))
    repeats_array = repeats[repeats$arrayID == arrays$arrayID[j], ]
    arrays$start[j] = min(repeats_array$start)
    arrays$end[j] = max(repeats_array$end)
  }
  
  
  ### Export ===============================================
  repeats$score[repeats$score == 0] <- 0.00000001
  export_gff(annotations.data.frame = repeats, output = ".", 
             file.name = paste0(assembly_name, "_repeats_filtered"), seqid = 1, source = "TRASH", type = "Satellite_DNA", 
             start = 3, end = 4, strand = 5, score = 6, attributes = c(12, 10, 13), 
             # start = 3, end = 4, strand = 5, score = 6, attributes = c(11, 10, 12), 
             attribute.names = c("Name=", "Family_EDS=", "Class_num_ID="))
  
  export_gff(annotations.data.frame = arrays, output = ".", 
             file.name = paste0(assembly_name, "_arrays_filtered"), seqid = 3, source = "TRASH", type = "Satellite_array", 
             start = 1, end = 2, score = 5, attributes = c(9, 10, 14), 
             attribute.names = c("Name=", "Repeat_no=", "Class_num_ID="))
  
  write.csv(repeats, paste0(assembly_name, "_repeats_filtered.csv"), row.names = FALSE)
  write.csv(arrays, paste0(assembly_name, "_arrays_filtered.csv"), row.names = FALSE)
  
  
  ### Plot distance to closest kmer for the whole genome ===
  # And GC
  # Kind of a dot plot along the sequence
  
  if(plot) {
    window_size_GC = 1000
    window_size_repeats = 2000
    max_repeat = 1000
    plotting_window_size = 1000000 # 1 mil
    if(!dir.exists("GC_repNo_repSize_plots")) dir.create("GC_repNo_repSize_plots")
    for(j in seq_along(names(assembly))) {
      print(paste0("I: ", i, " ", j, " / ", length(assembly)))
      plots_number = ceiling(length(assembly[[j]]) / plotting_window_size)
      
      pdf(file = paste0("./GC_repNo_repSize_plots/GC_repNo_repSize_", names(assembly)[j], ".pdf"), width = 15, height = plots_number * 2)
      par(mar = c(1,4,1,1), mfrow = c(plots_number * 2, 1))
      
      for(k in seq_len(plots_number)) {
        print(paste0("J: ", i, " sequence ", j, " / ", length(assembly), " plot ", k, " / ", plots_number))
        plotting_win_start = (k - 1) * plotting_window_size + 1
        plotting_win_end = (k) * plotting_window_size
        if(plotting_win_end > length(assembly[[j]])) plotting_win_end = length(assembly[[j]])
        
        substring = assembly[[j]][plotting_win_start : plotting_win_end]
        # Repeat fraction
        sequence_windows_starts <- genomic_bins_starts(start = plotting_win_start, end = plotting_win_end, bin_size = window_size_repeats)
        sequence_windows_ends <- c(sequence_windows_starts[2 : length(sequence_windows_starts)] - 1, plotting_win_end)
        if(length(sequence_windows_starts) == 1) sequence_windows_ends = plotting_win_end
        repeats_fraction = calculate.repeats.percentage.in.windows(windows.starts = sequence_windows_starts, 
                                                                   repeat.starts = repeats$start[repeats$seqID == names(assembly)[j]], 
                                                                   repeat.lengths = repeats$width[repeats$seqID == names(assembly)[j]], 
                                                                   sequence.length = plotting_win_end)
        plot(NULL, NULL, xlim = c(min(sequence_windows_starts), max(sequence_windows_ends)), ylim = c(0,100),
             xlab = "", ylab = "GC+rep")
        points(sequence_windows_starts, repeats_fraction, type = "h", col = "red")
        
        # GC
        sequence_windows_starts <- genomic_bins_starts(start = 1, end = length(substring), bin_size = window_size_GC)
        sequence_windows_ends <- c(sequence_windows_starts[2 : length(sequence_windows_starts)] - 1, length(substring))
        if(length(sequence_windows_starts) == 1) sequence_windows_ends = length(substring)
        gc_values = calculate.GC.in.windows(windows.starts = sequence_windows_starts, 
                                            windows.ends = sequence_windows_ends,
                                            sequence = substring)
        at_values = 100 - gc_values
        # gc_values = (gc_values - 40) * 2 + 40 # Blow up 2x around 40
        # at_values = (at_values - 60) * 2 + 60 # Blow up 2x around 60
        points((sequence_windows_starts + plotting_win_start - 1), gc_values, type = "l", col = "blue")
        points((sequence_windows_starts + plotting_win_start - 1), at_values, type = "l", col = "green")
        
        
        # Repeat median length
        sequence_windows_starts <- genomic_bins_starts(start = plotting_win_start, end = plotting_win_end, bin_size = window_size_repeats)
        sequence_windows_ends <- c(sequence_windows_starts[2 : length(sequence_windows_starts)] - 1, plotting_win_end)
        if(length(sequence_windows_starts) == 1) sequence_windows_ends = plotting_win_end
        repeats_N_in_win = calculate.repeats.sizes.in.windows(windows.starts = sequence_windows_starts, 
                                                              windows.ends = sequence_windows_ends, 
                                                              repeat.starts = repeats$start[repeats$seqID == names(assembly)[j]], 
                                                              repeat.lengths = repeats$width[repeats$seqID == names(assembly)[j]])
        plot(NULL, NULL, xlim = c(min(sequence_windows_starts), max(sequence_windows_ends)), ylim = c(0, max_repeat),
             xlab = "", ylab = "Rep sizes")
        points(sequence_windows_starts, repeats_N_in_win, type = "h", col = "blue")
        
        # Closest kmer distance
        # closest_distances_in_win = calculate_closest_dist_in_win(windows.starts = sequence_windows_starts - plotting_win_start + 1,
        #                                                          windows.ends = sequence_windows_ends - plotting_win_start + 1,
        #                                                          sequence = substring,
        #                                                          max_repeat = max_repeat)
        # plot(NULL, NULL, xlim = c(min(sequence_windows_starts), max(sequence_windows_ends)), ylim = c(0, max_repeat),
        #      xlab = "", ylab = "Closest N")
        # points(sequence_windows_starts, closest_distances_in_win, type = "h", col = "green")
        
      }
      dev.off()
      
    }
  }
  
  remove(assembly, repeats, arrays)
  gc()
  
}

