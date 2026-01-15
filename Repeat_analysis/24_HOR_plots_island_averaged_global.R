.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

library(ape)
library(nlme)

setwd("/home/pwlodzimierz/Rhynchospora/git_Rhynchospora")

tree_gene_based <- ape::read.tree(file = "SpeciesTree_rooted_node_labels.newick")

score_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/upload_data/23_hor_csvs_averaged_island_per_genome/", full.names = T)

windows = 25
score_per_window <- rep(0, windows)
repeats_per_window <- rep(0, windows)

for(i in seq_along(score_files)) {
  scores <- read.csv(file = score_files[i])
  for(j in seq_len(nrow(scores))) {
    score_per_window[j] = score_per_window[j] + scores$hor_novel_score_per_bin[j]
    repeats_per_window[j] = repeats_per_window[j] + scores$repeat_number_in_bin[j]
  }
}

averaged_scores <- score_per_window / repeats_per_window


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_plot_HOR_island_averaged_values_global.pdf"), 
    width = 7, height = 4)
plot(x = 1 : windows, y = averaged_scores, 
     main = paste0("All species HOR scores per island"), 
     xlab = paste0(windows, " windows averaged scores of ", sum(repeats_per_window), " repeats"),
     ylab = "HOR score: % of other repeats a repeat forms a HOR with", 
     cex.lab = 0.7, cex.axis = 0.7, cex.main = 1, cex.sub = 0.5, 
     type = "b")
dev.off()


###

order <- c("rugosa", "alba", "gaudichaudii", "cephalotes", 
           "ridleyi", 'watsonii', 'radicans', 'pubera', 
           'breviuscula', 'nervosa', 'ciliata', 'colorata', 
           'tenerrima', 'filiformis', 'austrobrasiliensis', 'tenuis', 
           'riparia', 'barbata', 'holoschoenoides', 'corymbosa')


repeat_chr_data <- read.csv("/home/pwlodzimierz/Rhynchospora/upload_data/16_scatter_plots/chromosomes_summary.csv")
repeat_chr_data <- repeat_chr_data[-c(49,50,52,53,55,56),]


windows = 25


chrsize_save <- NULL # in Mbp
te_save <- NULL
horscore_save <- NULL
species_save <- NULL
tyba_save <- NULL

for(i in seq_along(score_files)) {
  scores <- read.csv(file = score_files[i])
  score_per_window <- rep(0, windows)
  repeats_per_window <- rep(0, windows)
  for(j in seq_len(nrow(scores))) {
    score_per_window[j] = score_per_window[j] + scores$hor_novel_score_per_bin[j]
    repeats_per_window[j] = repeats_per_window[j] + scores$repeat_number_in_bin[j]
  }
  averaged_scores <- score_per_window / repeats_per_window
  averaged_score <- mean(averaged_scores)
  
  assembly_name_short = strsplit(strsplit(score_files[i], split = "_HOR_island_averaged_")[[1]][1], split = "island_per_genome//")[[1]][2]
  ordername <- order[unlist(lapply(order, function(X) grepl(X, score_files[i], ignore.case = T)))]
  
  if(grepl("tenuis", assembly_name_short)) assembly_name_short <- "Rtenuis.hap1.chr"
  if(grepl("radicans", assembly_name_short)) assembly_name_short <- "Rhync_radicans.asm.bp.p_ctg.FINAL.chr"
  
  chr_te_total_size <- sum(repeat_chr_data$TE_bp[repeat_chr_data$genome == assembly_name_short & repeat_chr_data$TE_bp != 0])
  chr_mean_size <- mean(repeat_chr_data$chromosome_size[repeat_chr_data$genome == assembly_name_short & repeat_chr_data$TE_bp != 0])
  chr_total_size <- sum(repeat_chr_data$chromosome_size[repeat_chr_data$genome == assembly_name_short & repeat_chr_data$TE_bp != 0])
  chr_tyba_total_size <- sum(repeat_chr_data$tyba_total_bp[repeat_chr_data$genome == assembly_name_short])
  
  te_percentage = 100 * chr_te_total_size / chr_total_size
  tyba_percentage = 100 * chr_tyba_total_size / chr_total_size
  
  te_save <- c(te_save, te_percentage)
  horscore_save <- c(horscore_save, averaged_score)
  species_save <- c(species_save, paste0("R",ordername))
  chrsize_save <- c(chrsize_save, chr_mean_size/1000000)
  tyba_save <- c(tyba_save, tyba_percentage)
}




# Create data frame
pgls_data <- data.frame(
  horscore = horscore_save,
  te = te_save,
  chrsize = chrsize_save,
  tyba = tyba_save,
  species = species_save
)
rownames(pgls_data) <- species_save

pgls_data$log10_te <- log10(pgls_data$te)
pgls_data$log10_chrsize <- log10(pgls_data$chrsize)
pgls_data$log10_tyba <- log10(pgls_data$tyba)
pgls_data$log10_horscore = log10(pgls_data$horscore)

# pgls_data <- pgls_data[-which(pgls_data$species == "Rbarbata"),]

# Fit individual PGLS models
model_te <- gls(log10_horscore ~ log10_te,
                correlation = corBrownian(phy = tree_gene_based, form = ~species),
                data = pgls_data)

model_chrsize <- gls(log10_horscore ~ log10_chrsize,
                     correlation = corBrownian(phy = tree_gene_based, form = ~species),
                     data = pgls_data)

model_tyba <- gls(log10_horscore ~ log10_tyba,
                  correlation = corBrownian(phy = tree_gene_based, form = ~species),
                  data = pgls_data)

# Extract p-values
p_values <- c(
  te = summary(model_te)$tTable["log10_te", "p-value"],
  chrsize = summary(model_chrsize)$tTable["log10_chrsize", "p-value"],
  tyba = summary(model_tyba)$tTable["log10_tyba", "p-value"]
)

# Fit individual PGLS models
model_te <- gls(horscore ~ te,
                correlation = corBrownian(phy = tree_gene_based, form = ~species),
                data = pgls_data)

model_chrsize <- gls(horscore ~ chrsize,
                     correlation = corBrownian(phy = tree_gene_based, form = ~species),
                     data = pgls_data)

model_tyba <- gls(horscore ~ tyba,
                  correlation = corBrownian(phy = tree_gene_based, form = ~species),
                  data = pgls_data)

# Extract p-values
p_values2 <- c(
  te = summary(model_te)$tTable["te", "p-value"],
  chrsize = summary(model_chrsize)$tTable["chrsize", "p-value"],
  tyba = summary(model_tyba)$tTable["tyba", "p-value"]
)

# 1 log corrected
model_te <- gls(horscore ~ log10_te,
                correlation = corBrownian(phy = tree_gene_based, form = ~species),
                data = pgls_data)

model_chrsize <- gls(horscore ~ log10_chrsize,
                     correlation = corBrownian(phy = tree_gene_based, form = ~species),
                     data = pgls_data)

model_tyba <- gls(horscore ~ log10_tyba,
                  correlation = corBrownian(phy = tree_gene_based, form = ~species),
                  data = pgls_data)

# Extract p-values
p_values6 <- c(
  te = summary(model_te)$tTable["log10_te", "p-value"],
  chrsize = summary(model_chrsize)$tTable["log10_chrsize", "p-value"],
  tyba = summary(model_tyba)$tTable["log10_tyba", "p-value"]
)

# non PGLS
p_values3 <- c(
  te = summary(lm(horscore ~ te, data = pgls_data))$coefficients[2, 4],
  chrsize = summary(lm(horscore ~ chrsize, data = pgls_data))$coefficients[2, 4],
  tyba = summary(lm(horscore ~ tyba, data = pgls_data))$coefficients[2, 4]
)
# non PGLS log log
p_values4 <- c(
  te = summary(lm(log10_horscore ~ log10_te, data = pgls_data))$coefficients[2, 4],
  chrsize = summary(lm(log10_horscore ~ log10_chrsize, data = pgls_data))$coefficients[2, 4],
  tyba = summary(lm(log10_horscore ~ log10_tyba, data = pgls_data))$coefficients[2, 4]
)
# non PGLS log
p_values5 <- c(
  te = summary(lm(log10_horscore ~ te, data = pgls_data))$coefficients[2, 4],
  chrsize = summary(lm(log10_horscore ~ chrsize, data = pgls_data))$coefficients[2, 4],
  tyba = summary(lm(log10_horscore ~ tyba, data = pgls_data))$coefficients[2, 4]
)

options(scipen = 999)
round(p_values2, 5) # non log
round(p_values, 5) # log-log
round(p_values6, 5) # log
round(p_values3, 5) # non log non corrected
round(p_values4, 5) # log-log non corrected
round(p_values5, 5) # log non corrected

if(F) {
  # Apply Bonferroni correction
  n_tests <- length(p_values)
  bonferroni_threshold <- 0.05 / n_tests
  p_adjusted <- p.adjust(p_values, method = "bonferroni")
  
  # Create results summary
  results_summary <- data.frame(
    predictor = names(p_values),
    beta = c(coef(model_te)["te"], 
             coef(model_chrsize)["chrsize"], 
             coef(model_tyba)["tyba"]),
    p_value = p_values,
    p_adjusted = p_adjusted,
    significant = p_adjusted < 0.05
  )
  
  print(results_summary)
  
  # Print individual model summaries
  cat("\n=== Model 1: horscore ~ te ===\n")
  summary(model_te)
  
  cat("\n=== Model 2: horscore ~ chrsize ===\n")
  summary(model_chrsize)
  
  cat("\n=== Model 3: horscore ~ tyba ===\n")
  summary(model_tyba)
  
  cat("\n=== Bonferroni corrected threshold: ", bonferroni_threshold, " ===\n")
  
}

# see if should transform

setwd("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots")
if(F) {
  library(ggplot2)
  library(gridExtra)
  
  # 1. Check residuals for all three models
  par(mfrow = c(3, 3))
  
  # Model 1: te
  plot(fitted(model_te), residuals(model_te, type = "normalized"),
       main = "TE: Residuals vs Fitted", xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, lty = 2)
  plot(pgls_data$te, pgls_data$horscore, main = "TE: Raw Data", 
       xlab = "TE", ylab = "horscore")
  abline(model_te, col = "red")
  qqnorm(residuals(model_te, type = "normalized"), main = "TE: Q-Q Plot")
  qqline(residuals(model_te, type = "normalized"))
  
  # Model 2: chrsize
  plot(fitted(model_chrsize), residuals(model_chrsize, type = "normalized"),
       main = "Chrsize: Residuals vs Fitted", xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, lty = 2)
  plot(pgls_data$chrsize, pgls_data$horscore, main = "Chrsize: Raw Data",
       xlab = "Chrsize", ylab = "horscore")
  abline(model_chrsize, col = "red")
  qqnorm(residuals(model_chrsize, type = "normalized"), main = "Chrsize: Q-Q Plot")
  qqline(residuals(model_chrsize, type = "normalized"))
  
  # Model 3: tyba
  plot(fitted(model_tyba), residuals(model_tyba, type = "normalized"),
       main = "Tyba: Residuals vs Fitted", xlab = "Fitted", ylab = "Residuals")
  abline(h = 0, lty = 2)
  plot(pgls_data$tyba, pgls_data$horscore, main = "Tyba: Raw Data",
       xlab = "Tyba", ylab = "horscore")
  abline(model_tyba, col = "red")
  qqnorm(residuals(model_tyba, type = "normalized"), main = "Tyba: Q-Q Plot")
  qqline(residuals(model_tyba, type = "normalized"))
  
  par(mfrow = c(1, 1))
  
  # 2. Try log transformations for models that show non-linearity
  # Log-transform predictors (add small constant if any zeros exist)
  pgls_data$log_te <- log(pgls_data$te + 0.01)
  pgls_data$log_chrsize <- log(pgls_data$chrsize + 0.01)
  pgls_data$log_tyba <- log(pgls_data$tyba + 0.01)
  
  # Fit models with log-transformed predictors
  model_te_log <- gls(horscore ~ log_te,
                      correlation = corBrownian(phy = tree_gene_based, form = ~species),
                      data = pgls_data)
  
  model_chrsize_log <- gls(horscore ~ log_chrsize,
                           correlation = corBrownian(phy = tree_gene_based, form = ~species),
                           data = pgls_data)
  
  model_tyba_log <- gls(horscore ~ log_tyba,
                        correlation = corBrownian(phy = tree_gene_based, form = ~species),
                        data = pgls_data)
  
  # Compare AIC values (lower is better)
  aic_comparison <- data.frame(
    model = c("te_linear", "te_log", "chrsize_linear", "chrsize_log", 
              "tyba_linear", "tyba_log"),
    AIC = c(AIC(model_te), AIC(model_te_log),
            AIC(model_chrsize), AIC(model_chrsize_log),
            AIC(model_tyba), AIC(model_tyba_log))
  )
  
  print("=== AIC Comparison (lower is better) ===")
  print(aic_comparison)
  
  # Summaries of log models
  cat("\n=== Log-transformed Models ===\n")
  summary(model_te_log)
  summary(model_chrsize_log)
  summary(model_tyba_log)
  
  # Extract p-values for log models and compare
  p_values_log <- c(
    te_log = summary(model_te_log)$tTable["log_te", "p-value"],
    chrsize_log = summary(model_chrsize_log)$tTable["log_chrsize", "p-value"],
    tyba_log = summary(model_tyba_log)$tTable["log_tyba", "p-value"]
  )
  
  p_fdr_log <- p.adjust(p_values_log, method = "fdr")
  
  results_log <- data.frame(
    predictor = names(p_values_log),
    beta = c(coef(model_te_log)["log_te"],
             coef(model_chrsize_log)["log_chrsize"],
             coef(model_tyba_log)["log_tyba"]),
    p_raw = p_values_log,
    p_fdr = p_fdr_log,
    sig_fdr = p_fdr_log < 0.05
  )
  
  print("\n=== Log-transformed Results with FDR ===")
  print(results_log)
}



### replot the 3 plots with linear regression


order2 <- c("Rrugosa", "Ralba", "Rgaudichaudii", "Rcephalotes", 
            "Rridleyi", 'Rwatsonii', 'Rradicans', 'Rpubera', 
            'Rbreviuscula', 'Rnervosa', 'Rciliata', 'Rcolorata', 
            'Rtenerrima', 'Raustrobrasiliensis', 'Rtenuis', 
            'Rriparia', 'Rbarbata', 'Rholoschoenoides', 'Rcorymbosa')

pgls_data <- pgls_data[match(order2, pgls_data$species), ]

pgls_data$pch =  c(15, 16, 17, 18,
                   15, 16, 17, 18,
                   15, 16, 17, 18,
                   15, 17, 18,
                   15, 16, 17, 18)

pgls_data$col = c("#E69F00", 
                  "#D55E00", 
                  "#D95F02", 
                  "#FC8D62", 
                  "#E5C494", 
                  "#F0E442", 
                  "#FFD92F", 
                  "#A6D854", 
                  "#66C2A5", 
                  "#009E73", 
                  "#1B9E77", 
                  "#66A61E", 
                  "#CC79A7", 
                  "#E7298A", 
                  "#B3B3CC", 
                  "#8DA0CB", 
                  "#7570B3", 
                  "#56B4E9", 
                  "#0072B2")


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_hor_score_per_tyba_percentage_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$tyba, y = pgls_data$horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per Tyba percentage"), 
     xlab = paste0("Genome wide Tyba percentage"),
     ylab = "Mean HOR score", 
     type = "b", ylim = c(0,60), cex = 2)
# abline(gls(horscore ~ tyba,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(horscore ~ tyba, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_hor_score_per_mean_chr_size_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$chrsize, y = pgls_data$horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per mean chromosome size"), 
     xlab = paste0("Mean chromosome size, Mbp"),
     ylab = "Mean HOR score", 
     type = "b", ylim = c(0,60), cex = 2)
# abline(gls(horscore ~ chrsize,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(horscore ~ chrsize, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_hor_score_per_TE_percentage_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$te, y = pgls_data$horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per TE percentage"), 
     xlab = paste0("Genome wide TE percentage"),
     ylab = "Mean HOR score", 
     type = "b", ylim = c(0,60), cex = 2)
# abline(gls(horscore ~ te,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(horscore ~ te, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()


### now log transformed 

pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_LOG10_hor_score_per_tyba_percentage_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$tyba, y = pgls_data$log10_horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per Tyba percentage"), 
     xlab = paste0("Genome wide Tyba percentage"),
     ylab = "LOG10 Mean HOR score", 
     type = "b", cex = 2)
# abline(gls(log10_horscore ~ tyba,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(log10_horscore ~ tyba, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_LOG10_hor_score_per_mean_chr_size_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$chrsize, y = pgls_data$log10_horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per mean chromosome size"), 
     xlab = paste0("Mean chromosome size, Mbp"),
     ylab = "LOG10 Mean HOR score", 
     type = "b", cex = 2)
# abline(gls(log10_horscore ~ chrsize,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(log10_horscore ~ chrsize, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()


pdf(file = paste0("/home/pwlodzimierz/Rhynchospora/upload_data/24_residuals_correlation_plots/24_plot_HOR_island_mean_LOG10_hor_score_per_TE_percentage_regline.pdf"), 
    width = 7, height = 7)
plot(x = pgls_data$te, y = pgls_data$log10_horscore, pch = pgls_data$pch, col = pgls_data$col, lty = 0,
     main = paste0("All species HOR scores per per TE percentage"), 
     xlab = paste0("Genome wide TE percentage"),
     ylab = "LOG10 Mean HOR score", 
     type = "b", cex = 2)
# abline(gls(log10_horscore ~ te,
#            correlation = corBrownian(phy = tree_gene_based, form = ~species),
#            data = pgls_data), col = "gray", lwd = 2, lty = 1)
abline(lm(log10_horscore ~ te, data = pgls_data ), col = "black", lwd = 2, lty = 1)
dev.off()

