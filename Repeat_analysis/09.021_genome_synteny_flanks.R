


input_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
input_id = as.numeric(input_id)# 1 to 37
print(input_id)
input_id <- as.numeric(input_id)
# input_id = 3

keep_chrs_austr <- c(1,4,7)
keep_chrs_austr2 <- c(1,4,7)

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(Biostrings))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(patchwork))


def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

suppressMessages(library(seqinr))

runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv")

repeats_ex1 <- runs_todo$repeats_file_A[input_id]
repeats_ex2 <- runs_todo$repeats_file_B[input_id]

arrays_ex1 <- runs_todo$arrays_file_A[input_id]
arrays_ex2 <- runs_todo$arrays_file_B[input_id]

genome_ex1 <- runs_todo$genome_A[input_id]
genome_ex2 <- runs_todo$genome_B[input_id]

genome_A <- genome_ex1
genome_B <- genome_ex2


# repeats_ex2 <- runs_todo$repeats_file_A[input_id]
# repeats_ex1 <- runs_todo$repeats_file_B[input_id]
# 
# arrays_ex2 <- runs_todo$arrays_file_A[input_id]
# arrays_ex1 <- runs_todo$arrays_file_B[input_id]
# 
# genome_ex2 <- runs_todo$genome_A[input_id]
# genome_ex1 <- runs_todo$genome_B[input_id]

flanking_range = 7500
flanging_skip = 2500




flank_region_matrix <- function(repeats_ex1, repeats_ex2, arrays_ex1, arrays_ex2, genome_ex1, genome_ex2) {
  
  
  flanking_range = 7500
  flanging_skip = 2500
  
  # make chromosome sequences into the /home/pwlodzimierz/Rhynchospora/upload_data/09.02_data if they don't exist
  
  
  query_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_data/query_", flanking_range, "_", flanging_skip, "_", genome_ex1, ".fa")
  subject_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_data/subject_", genome_ex2, ".fa")
  sam_file <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_data/", flanking_range, "_", flanging_skip, "_", genome_ex1, "_", genome_ex2, ".sam")
  
  chr_sizes <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")
  chr_1_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_ex1]
  chr_2_sizes <- chr_sizes$size[chr_sizes$assembly.name == genome_ex2]
  
  chr_1_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_ex1]
  chr_2_names <- chr_sizes$chromosome.name[chr_sizes$assembly.name == genome_ex2]
  
  arrays_1 <- read.csv(arrays_ex1)
  
  arrays_2 <- read.csv(arrays_ex2)
  
  if(genome_ex1 == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    chr_1_sizes <- chr_1_sizes[keep_chrs]
    chr_1_names <- chr_1_names[keep_chrs]
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names,]
  }
  if(genome_ex2 == "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta") {
    keep_chrs <- keep_chrs_austr
    chr_2_sizes <- chr_2_sizes[keep_chrs]
    chr_2_names <- chr_2_names[keep_chrs]
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names,]
  }
  if(genome_ex1 == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    chr_1_sizes <- chr_1_sizes[keep_chrs]
    chr_1_names <- chr_1_names[keep_chrs]
    arrays_1 <- arrays_1[arrays_1$chromosome %in% chr_1_names,]
  }
  if(genome_ex2 == "Rhync_austrobrasiliensis_6344D.hic.hap2.chr.fasta") {
    keep_chrs <- keep_chrs_austr2
    chr_2_sizes <- chr_2_sizes[keep_chrs]
    chr_2_names <- chr_2_names[keep_chrs]
    arrays_2 <- arrays_2[arrays_2$chromosome %in% chr_2_names,]
  }
  
  
  chr_1_starts_adjust <- 0
  for(i in 1 : (length(chr_1_sizes) - 1)) {
    chr_1_starts_adjust <- c(chr_1_starts_adjust, sum(chr_1_sizes[1 : (i)]))
  }
  chr_2_starts_adjust <- 0
  for(i in 1 : ( length(chr_2_sizes) - 1 ) ) {
    chr_2_starts_adjust <- c(chr_2_starts_adjust, sum(chr_2_sizes[1 : (i)]))
  }
  
  
  for(i in 1 : length(chr_1_sizes)) {
    arrays_1$start[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$start[arrays_1$chromosome == chr_1_names[i]] + chr_1_starts_adjust[i]
    arrays_1$end[arrays_1$chromosome == chr_1_names[i]] <- arrays_1$end[arrays_1$chromosome == chr_1_names[i]] + chr_1_starts_adjust[i]
  }
  
  
  for(i in 1 : length(chr_2_sizes)) {
    arrays_2$start[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$start[arrays_2$chromosome == chr_2_names[i]] + chr_2_starts_adjust[i]
    arrays_2$end[arrays_2$chromosome == chr_2_names[i]] <- arrays_2$end[arrays_2$chromosome == chr_2_names[i]] + chr_2_starts_adjust[i]
  }
 
  if(!file.exists(query_file)) { # create the query chromosome file
    print("reading query assembly")
    assemblies <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", recursive = T, full.names = T)
    assemblies <- assemblies[grep(genome_ex1, assemblies)[1]]
    assembly <- readDNAStringSet(assemblies)
    
    merged_seq <- DNAString(paste0(as.character(assembly[chr_1_names]), collapse = ""))
    merged_set <- DNAStringSet(merged_seq)
    names(merged_set) <- paste0("merged_", genome_ex1)
    
    print("saving query sequences")
    
    for(i in 1 : nrow(arrays_1)) {
      if((arrays_1$start[i] - flanking_range - flanging_skip) < 0) {
        write.fasta(sequences = merged_set[[1]][(arrays_1$end[i] + flanging_skip) : (arrays_1$end[i] + flanging_skip + flanking_range)],
                    names = paste0("down_", i), file.out = query_file, open = "a")
      } else if(arrays_1$end[i] + flanging_skip + flanking_range > width(merged_set)) {
        write.fasta(sequences = merged_set[[1]][(arrays_1$start[i] - flanking_range - flanging_skip) : (arrays_1$start[i] - flanging_skip)],
                    names = paste0("up_", i), file.out = query_file, open = "a")
      } else {
        write.fasta(sequences = merged_set[[1]][(arrays_1$start[i] - flanking_range - flanging_skip) : (arrays_1$start[i] - flanging_skip)],
                    names = paste0("up_", i), file.out = query_file, open = "a")
        write.fasta(sequences = merged_set[[1]][(arrays_1$end[i] + flanging_skip) : (arrays_1$end[i] + flanging_skip + flanking_range)],
                    names = paste0("down_", i), file.out = query_file, open = "a")
      }
      
      
    }
    
  }
  
  if(!file.exists(subject_file)) { # create the subject chromosome file
    print("reading subject assembly")
    assemblies <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies", recursive = T, full.names = T)
    assemblies <- assemblies[grep(genome_ex2, assemblies)[1]]
    assembly <- readDNAStringSet(assemblies)
    
    merged_seq <- DNAString(paste0(as.character(assembly[chr_2_names]), collapse = ""))
    merged_set <- DNAStringSet(merged_seq)
    names(merged_set) <- paste0("merged_", genome_ex2)
    
    print("saving subject sequences")
    
    writeXStringSet(merged_set, subject_file)
    
    
    
  }
  
  
  if(!file.exists(sam_file)) { # make the minimap2 alignment
    
    print("making minimap2 alignment")
    
    system2("minimap2", args = c(
      "-a", "-x", "asm20", "--secondary=no", "--eqx",
      "-k15", "-w5",
      "--max-chain-skip", "80",
      "--min-occ-floor", "20",
      "--mask-level", "0.1",
      subject_file, query_file
    ), stdout = sam_file)
    
    print("alignment done")
    
  } 
  
  
  
  
  
  
  # --- PARSE SAM ---
  # Read lines from SAM file and remove headers
  sam <- readLines(sam_file)
  sam <- sam[!grepl("^@", sam)]  # Remove @HD, @SQ, etc.
  
  # Split each line into fields
  split_sam <- strsplit(sam, "\t")
  
  # Filter out clearly invalid lines (fewer than 11 fields)
  split_sam <- split_sam[sapply(split_sam, length) >= 11]
  
  # Normalize length: pad rows with NA to make all rows equal length
  max_len <- max(sapply(split_sam, length))
  split_sam <- lapply(split_sam, function(x) c(x, rep(NA, max_len - length(x))))
  
  # Combine into data.frame
  sam_df <- as.data.frame(do.call(rbind, split_sam), stringsAsFactors = FALSE)
  
  colnames(sam_df)[1:11] <- c("query", "flag", "subject", "subject_start", "mapq", "cigar",
                              "rnext", "pnext", "tlen", "seq", "qual")
  
  sam_df$subject_start <- as.integer(sam_df$subject_start)
  sam_df$flag <- as.integer(sam_df$flag)
  sam_df$mapq <- as.integer(sam_df$mapq)
  sam_df$strand <- ifelse(bitwAnd(sam_df$flag, 16) != 0, "-", "+")
  
  sam_df$alignment_score <- as.numeric(stringr::str_extract(
    apply(sam_df[ , 12:ncol(sam_df)], 1, paste, collapse = "\t"),
    "(?<=AS:i:)\\d+"
  ))
  
  # filter hits
  
  # sam_df <- sam_df[!is.na(sam_df$alignment_score), ]
  
  # sam_df <- sam_df[sam_df$mapq >= 30, ]
  
  
  # --- Compute Alignment Length from CIGAR ---
  parse_cigar <- function(cigar) {
    ops <- unlist(regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar)))
    ops_df <- do.call(rbind, lapply(ops, function(x) {
      val <- as.integer(gsub("[A-Z=]", "", x))
      op  <- gsub("[0-9]", "", x)
      c(op = op, val = val)
    }))
    ops_df <- as.data.frame(ops_df, stringsAsFactors = FALSE)
    ops_df$val <- as.integer(ops_df$val)
    ops_df
  }
  
  get_query_length <- function(cigar) {
    ops <- parse_cigar(cigar)
    sum(ops$val[ops$op %in% c("M", "I", "S", "=", "X")])
  }
  
  get_subject_length <- function(cigar) {
    ops <- parse_cigar(cigar)
    sum(ops$val[ops$op %in% c("M", "D", "=", "X")])  # drop "N"
  }
  
  # Add alignment ends
  sam_df$query_aln_len   <- sapply(sam_df$cigar, get_query_length)
  sam_df$subject_aln_len <- sapply(sam_df$cigar, get_subject_length)
  
  # Compute query and subject end positions
  sam_df$query_start <- 1  # minimap2 always maps full query unless soft clipped
  sam_df$query_end <- sam_df$query_start + sam_df$query_aln_len - 1
  sam_df$subject_end <- sam_df$subject_start + sam_df$subject_aln_len - 1
  
  # --- Keep top n hits per query ---
  top_hits <- sam_df %>%
    group_by(query) %>%
    arrange(desc(alignment_score), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
  # Prepare match_data from top_hits
  match_data <- top_hits %>%
    mutate(
      type = str_extract(query, "^(up|down|both)"),
      index = as.integer(str_extract(query, "(?<=_)[0-9]+"))
    ) %>%
    left_join(arrays_1 %>% mutate(index = row_number()), by = "index")
  
  
  
  
  # Step 1: Extract match type and index from query names
  match_data <- match_data %>%
    mutate(query_flank_origin = case_when(
      type == "up"   ~ start - flanking_range - flanging_skip,
      type == "down" ~ end   + flanging_skip,
      type == "both" ~ start - flanking_range - flanging_skip  # merged = from upstream start
    ))
  
  match_data <- match_data %>%
    mutate(
      abs_query_start = query_flank_origin + query_start - 1,
      abs_query_end   = query_flank_origin + query_end - 1
    )
  
  match_data <- match_data %>%
    filter(!is.na(subject_start), subject_start > 0, !is.na(abs_query_start), !is.na(abs_query_end))
  
  # calculate normalised score that accounts for alignment length: score_frac = (alignment_score / alignment_length) * (alignment_length / query_length)
  
  match_data <- match_data %>%
    mutate(
      alignment_length = sapply(cigar, function(c) {
        matches <- as.numeric(unlist(str_extract_all(c, "\\d+(?=M|=|X)")))
        sum(matches, na.rm = TRUE)
      }),
      query_length = nchar(seq),
      score_frac = (alignment_score / alignment_length) * (alignment_length / query_length)
    )
  
  
  
  # save the data
  match_data <- match_data[, -which(names(match_data) == "seq")]
  write.csv(x = match_data, file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.021_out_data/", flanking_range, "_", flanging_skip, "_", genome_ex1, "_", genome_ex2, ".csv"))
  
  
  # pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/09.02_synteny_plots/hist", flanking_range, "_", flanging_skip, "_", genome_ex1, "_", genome_ex2, ".pdf"), 
  #     width = 8, height = 4, onefile = T)
  # hist(match_data$score_frac, breaks = 100, main = "Normalized Alignment Scores", xlab = "score_frac")
  # hist(match_data$subject_start, breaks = 100, main = "subject_start")
  # hist(match_data$subject_end, breaks = 100, main = "subject_end")
  # hist(match_data$abs_query_start, breaks = 100, main = "abs_query_start")
  # hist(match_data$abs_query_end, breaks = 100, main = "abs_query_end")
  # dev.off()
  
  return(match_data)
  
  
}




match_data_A = flank_region_matrix(repeats_ex1, repeats_ex2, arrays_ex1, arrays_ex2, genome_ex1, genome_ex2)

match_data_B = flank_region_matrix(repeats_ex2, repeats_ex1, arrays_ex2, arrays_ex1, genome_ex2, genome_ex1)

