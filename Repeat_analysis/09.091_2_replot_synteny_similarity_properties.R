


setwd("/home/pwlodzimierz/Rhynchospora/upload_data/9.091_mirrored_boxplots")

data_files <- list.files(pattern = "csv")

species <- c('rugosa', 'alba', 'cephalotes', 'ridleyi', 'gaudichaudii', 
             'watsonii', 'radicans', 'pubera', 'breviuscula', 'nervosa', 
             'ciliata', 'colorata', 'tenerrima', 'filiformis', 'austrobrasiliensis', 
             'tenuis', 'riparia', 'barbata', 'holoschoenoides', 'corymbosa')

all_data <- data.frame()
for(i in seq_along(data_files)) {
  data <- read.csv(file = data_files[i])
  
  data$species_A_id <- which(unlist(lapply(species, function(x) {grepl(x, data$species_A)})))
  data$species_B_id <- which(unlist(lapply(species, function(x) {grepl(x, data$species_B)})))
  
  
  
  all_data <- rbind(all_data, data)
}


colors_vector <- c("#E69F00", "#D55E00", 
                   "#FC8D62", "#E5C494", "#D95F02", 
                   "#F0E442", "#FFD92F", 
                   "#A6D854", "#66C2A5", 
                   "#009E73", "#1B9E77", 
                   "#66A61E", "#CC79A7", 
                   "#E78AC3", "#E7298A", 
                   "#B3B3CC", "#8DA0CB", "#7570B3", 
                   "#56B4E9", "#0072B2")

colors_vector_transp <- c("#E69F0095", "#D55E0095", 
                          "#FC8D6295", "#E5C49495", "#D95F0295", "#F0E44295", 
                          "#FFD92F95", "#A6D85495", "#66C2A595", "#009E7395", "#1B9E7795", 
                          "#66A61E95", "#CC79A795", "#E78AC395", "#E7298A95", 
                          "#B3B3CC95", "#8DA0CB95", "#7570B395", "#56B4E995", "#0072B295")

pch_vector <- c(15, 16, 18,
                15, 17, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18)



between_species_data <- all_data[!all_data$haplotype_comparison, ]
between_species_data$species_A_id <- as.numeric(between_species_data$species_A_id)
between_species_data <- between_species_data[order(between_species_data$species_A_id),]

yadjust = 100
max_y_value <- max(c(all_data$arr_num_A, all_data$arr_num_B)) + yadjust + yadjust
points_cex = 1.5

pdf(file = "species_barplot.pdf", 
    width = 8, height = 8)
par(mar = c(2,4,2,4))

plot(x = NULL, y = NULL, xlim = c(0,19), ylim = c(-max_y_value-yadjust,max_y_value+yadjust), 
     xaxt = "n", xlab = "", yaxt = "n", ylab = "")

for(i in 1 : nrow(between_species_data)) {
  
  # lines(x = c(i-0.3,i+0.3), y = rep(max_y_value * between_species_data$similarity_values_A[i]/100, 2), 
  #       col = "green")
  points(x = i, y = max_y_value * between_species_data$similarity_values_A[i]/100, 
         col = "green", pch = 16)
  
  points(x = i, y = 50, pch = pch_vector[i], col = colors_vector[i])
  points(x = i, y = -50, pch = pch_vector[i+1], col = colors_vector[i+1])
  
  
  rect(xleft = i - 0.3, ybottom = yadjust + between_species_data$arr_A_found_as_array_in_B[i], xright = i + 0.3, 
       ytop = yadjust + between_species_data$arr_A_found_as_missing_in_B[i] + between_species_data$arr_A_found_as_array_in_B[i], 
       border = "grey", col = "grey", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = yadjust, xright = i + 0.3, 
       ytop = yadjust + between_species_data$arr_A_found_as_array_in_B[i], 
       border = "orange", col = "orange", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = yadjust, xright = i + 0.3, 
       ytop = yadjust + between_species_data$arr_num_A[i], 
       border = "black")
  
  
  rect(xleft = i - 0.3, ybottom = -yadjust - between_species_data$arr_B_found_as_array_in_A[i], xright = i + 0.3, 
       ytop = -yadjust - between_species_data$arr_B_found_as_missing_in_A[i] - between_species_data$arr_B_found_as_array_in_A[i], 
       border = "grey", col = "grey", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = -yadjust, xright = i + 0.3, 
       ytop = -yadjust - between_species_data$arr_B_found_as_array_in_A[i], 
       border = "orange", col = "orange", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = -yadjust, xright = i + 0.3, 
       ytop = -yadjust - between_species_data$arr_num_B[i], 
       border = "black")
  
  
  
}
text(x = 9, y = 3150, "All arrays", adj = 0)
text(x = 9, y = 2950, "Arrays with a syntenic match", adj = 0)
text(x = 9, y = 2750, "Arrays with syntenic array match", adj = 0)
rect(xleft = 8, ybottom = 3100, xright = 8.9, ytop = 3200, density = 30, col = "white", border = "black")
rect(xleft = 8, ybottom = 2900, xright = 8.9, ytop = 3000, density = 30, col = "grey", border = "grey")
rect(xleft = 8, ybottom = 2700, xright = 8.9, ytop = 2800, density = 30, col = "orange", border = "orange")

axis(4, at = 100 + (seq(0, 100, by = 20) / 100) * (max_y_value - 100), labels = seq(0, 100, by = 20), las = 1)
mtext("Syntenic arrays % similarity", side = 4, line = 3, col = "green", adj = 0.8)

axis(2, at = 100 + seq(0, 3000, by = 500) , labels = seq(0, 3000, by = 500), las = 1)
axis(2, at = -100 - seq(0, 3000, by = 500) , labels = seq(0, 3000, by = 500), las = 1)
mtext("Count of arrays", side = 2, line = 3, col = "black", adj = 0.5)

dev.off()



between_haplotypes_data <- all_data[all_data$haplotype_comparison, ]
between_haplotypes_data$X = 1:nrow(between_haplotypes_data)
between_haplotypes_data <- between_haplotypes_data[c(4,11,9,3,10,7,5,6,13,12,1,2,8),]
yadjust = 100
max_y_value <- max(c(all_data$arr_num_A, all_data$arr_num_B)) + yadjust + yadjust
points_cex = 1.5

pdf(file = "haplotypes_barplot.pdf", 
    width = 8, height = 8)
par(mar = c(2,4,2,4))

plot(x = NULL, y = NULL, xlim = c(0,12), ylim = c(-max_y_value-yadjust,max_y_value+yadjust), 
     xaxt = "n", xlab = "", yaxt = "n", ylab = "")

for(i in 1 : nrow(between_haplotypes_data)) {
  
  # lines(x = c(i-0.3,i+0.3), y = rep(max_y_value * between_haplotypes_data$similarity_values_A[i]/100, 2), 
  #       col = "green")
  points(x = i, y = max_y_value * between_haplotypes_data$similarity_values_A[i]/100, 
        col = "green", pch = 16)
  
  which_species <- which(unlist(lapply(1 : length(species), function(X) grepl(pattern = species[X], x = between_haplotypes_data$species_A[i]))))
  
  points(x = i, y = 0, pch = pch_vector[which_species], col = colors_vector[which_species])
  
  
  rect(xleft = i - 0.3, ybottom = yadjust + between_haplotypes_data$arr_A_found_as_array_in_B[i], xright = i + 0.3, 
       ytop = yadjust + between_haplotypes_data$arr_A_found_as_missing_in_B[i] + between_haplotypes_data$arr_A_found_as_array_in_B[i], 
       border = "grey", col = "grey", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = yadjust, xright = i + 0.3, 
       ytop = yadjust + between_haplotypes_data$arr_A_found_as_array_in_B[i], 
       border = "orange", col = "orange", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = yadjust, xright = i + 0.3, 
       ytop = yadjust + between_haplotypes_data$arr_num_A[i], 
       border = "black")
  
  
  rect(xleft = i - 0.3, ybottom = -yadjust - between_haplotypes_data$arr_B_found_as_array_in_A[i], xright = i + 0.3, 
       ytop = -yadjust - between_haplotypes_data$arr_B_found_as_missing_in_A[i] - between_haplotypes_data$arr_B_found_as_array_in_A[i], 
       border = "grey", col = "grey", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = -yadjust, xright = i + 0.3, 
       ytop = -yadjust - between_haplotypes_data$arr_B_found_as_array_in_A[i], 
       border = "orange", col = "orange", lty = 1, density = 30, lwd = 1)
  rect(xleft = i - 0.3, ybottom = -yadjust, xright = i + 0.3, 
       ytop = -yadjust - between_haplotypes_data$arr_num_B[i], 
       border = "black")
  
  
  
}
text(x = 7, y = 3150, "All arrays", adj = 0)
text(x = 7, y = 2950, "Arrays with a syntenic match", adj = 0)
text(x = 7, y = 2750, "Arrays with syntenic array match", adj = 0)
rect(xleft = 6, ybottom = 3100, xright = 6.9, ytop = 3200, density = 30, col = "white", border = "black")
rect(xleft = 6, ybottom = 2900, xright = 6.9, ytop = 3000, density = 30, col = "grey", border = "grey")
rect(xleft = 6, ybottom = 2700, xright = 6.9, ytop = 2800, density = 30, col = "orange", border = "orange")

axis(4, at = 100 + (seq(0, 100, by = 20) / 100) * (max_y_value - 100), labels = seq(0, 100, by = 20), las = 1)
mtext("Syntenic arrays % similarity", side = 4, line = 3, col = "green", adj = 0.8)

axis(2, at = 100 + seq(0, 3000, by = 500) , labels = seq(0, 3000, by = 500), las = 1)
axis(2, at = -100 - seq(0, 3000, by = 500) , labels = seq(0, 3000, by = 500), las = 1)
mtext("Count of arrays", side = 2, line = 3, col = "black", adj = 0.5)


dev.off()












