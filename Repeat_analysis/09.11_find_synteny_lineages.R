
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
# library(seqinr)
# library(GenomicRanges)
# library(IRanges)

runs_todo <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora/synteny_comparisons_table_full_genomes.csv")
lineages_list <- list()
current_lineage_ids <- data.frame()
next_lineage_id <- 1

for(input_id in 15 : 33) {
  genome_A <- runs_todo$genome_A[input_id]
  genome_B <- runs_todo$genome_B[input_id]
  
  file_path <- paste0("/home/pwlodzimierz/Rhynchospora/upload_data/9.031_data/synteny_pairs_", genome_A, "_NA_", genome_B, "_NA.csv")
  
  synteny_pairs_df <- read.csv(file = file_path)

  synteny_pairs_df <- synteny_pairs_df[synteny_pairs_df$is_island_A & synteny_pairs_df$is_island_B, ]
  
  # If first iteration or no active lineages
  if(nrow(current_lineage_ids) == 0) {
    synteny_pairs_df$lineage_id <- seq(next_lineage_id, length.out = nrow(synteny_pairs_df))
    next_lineage_id <- next_lineage_id + nrow(synteny_pairs_df)
    
    # Add A side
    lineages_list[[length(lineages_list) + 1]] <- data.frame(
      lineage_id = synteny_pairs_df$lineage_id,
      genome = genome_A,
      array_id = synteny_pairs_df$island_A_ID
    )
    # Add B side
    lineages_list[[length(lineages_list) + 1]] <- data.frame(
      lineage_id = synteny_pairs_df$lineage_id,
      genome = genome_B,
      array_id = synteny_pairs_df$island_B_ID
    )
    
    current_lineage_ids <- data.frame(array_id = synteny_pairs_df$island_B_ID,
                                      lineage_id = synteny_pairs_df$lineage_id)
  } else {
    # Match A of current with B of previous
    matches <- match(synteny_pairs_df$island_A_ID, current_lineage_ids$array_id)
    
    # Assign existing lineage IDs
    synteny_pairs_df$lineage_id <- current_lineage_ids$lineage_id[matches]
    
    # Handle new lineages (NA matches)
    new_indices <- which(is.na(synteny_pairs_df$lineage_id))
    if(length(new_indices) > 0) {
      synteny_pairs_df$lineage_id[new_indices] <- seq(next_lineage_id, length.out = length(new_indices))
      next_lineage_id <- next_lineage_id + length(new_indices)
      
      # Add A side for NEW lineages
      lineages_list[[length(lineages_list) + 1]] <- data.frame(
        lineage_id = synteny_pairs_df$lineage_id[new_indices],
        genome = genome_A,
        array_id = synteny_pairs_df$island_A_ID[new_indices]
      )
    }
    
    # Add B side for ALL lineages
    lineages_list[[length(lineages_list) + 1]] <- data.frame(
      lineage_id = synteny_pairs_df$lineage_id,
      genome = genome_B,
      array_id = synteny_pairs_df$island_B_ID
    )
    
    # Update current_lineage_ids
    current_lineage_ids <- data.frame(array_id = synteny_pairs_df$island_B_ID,
                                      lineage_id = synteny_pairs_df$lineage_id)
  }
}

lineages_df <- do.call(rbind, lineages_list)

lineages_summary <- table(lineages_df$lineage_id)

lineages_summary <- lineages_summary[lineages_summary > 2]

lineages_summary <- data.frame(lineages_summary)

# Add random colours avoiding greyscale
set.seed(123)
lineages_summary$colour <- hsv(h = runif(nrow(lineages_summary)), 
                               s = runif(nrow(lineages_summary), 0.5, 1), 
                               v = runif(nrow(lineages_summary), 0.5, 1))

write.csv(file = "/home/pwlodzimierz/Rhynchospora/upload_data/9.11_data/synteny_lineages_full_genomes.csv",
          x = lineages_df[lineages_df$lineage_id %in% lineages_summary$Var1, ],
          row.names = FALSE)
write.csv(file = "/home/pwlodzimierz/Rhynchospora/upload_data/9.11_data/synteny_lineages_summary_full_genomes.csv",
          x = lineages_summary,
          row.names = FALSE)

pdf(file = "/home/pwlodzimierz/Rhynchospora/upload_data/9.11_data/synteny_lineages_summary_histogram.pdf")
histo <- hist(lineages_summary$Freq, breaks = 2:8, 
     main = "", ylab = "Count of array lineages", xlab = "Genomes number", xaxt="n")
axis(side = 1, at = histo$mids,labels = 3:8)
dev.off()

