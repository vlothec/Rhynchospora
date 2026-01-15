.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

def_wd = "/home/pwlodzimierz/Rhynchospora/git_Rhynchospora"
setwd(def_wd)
source("./aux_fun.R")

library(msa)
library(ggplot2)
library(GenomicRanges)

setwd("/home/pwlodzimierz/Rhynchospora/local_analysis")





data_directories <- list.dirs(path = "/home/pwlodzimierz/Rhynchospora/TRASH_2025", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[!grepl(pattern = "asiliensis_6344D.hic.hap1.chr.1.4.7.fast", data_directories)]
data_directories <- data_directories[!grepl(pattern = "Rhync_tenuis_6", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2024", recursive = FALSE, full.names = TRUE)
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2022", recursive = FALSE, full.names = TRUE))
assembly_files <- c(assembly_files, list.files(path = "/home/pwlodzimierz/Rhynchospora/assemblies/fastas_Andre_2025", recursive = FALSE, full.names = TRUE))

data_directories <- data_directories[c(1,3,5,6,7,10,11,13,15,17,19,24,27,29,33,34,26,36,37,38)]
#

islands <- NULL
for(i in 1 : length(data_directories)) {
  print(i)
  
  setwd(data_directories[i])
  
  
  islands_file = list.files(pattern = "islands_genome", full.names = TRUE)
  if(length(islands_file) > 1) {
    islands_file <- islands_file[grep("islands_genome.csv", islands_file)]
  }
  
  islands <- rbind(islands, read.csv(file = islands_file))
}
#   "R_barbata_hap1_chrs.fa"
#   "R_holoschoenoides_hap1_chr.fasta"
#   "R_pubera_ref_2n10_v2.chr.fasta"
#   "Ralba.chr.fasta"
#   "Rbreviuscula.hap1.chr.fasta"
#   "Rcephalotes.hap2.chr.fasta"
#   "Rcolorata.hap1.chr.fasta"
#   "Rhync_austrobrasiliensis_6344D.hic.hap1.chr.fasta"
#   "Rhync_ciliata_6041A.hic.hap1.out_JBAT.FINAL.chr.fasta"
#   "Rhync_corymbosa_6179D.asm.hic.FINAL.hap1.chr.fasta"
#   "Rhync_nervosa_6321B.asm.hic.hap1.p_ctg.FINAL.chr.fasta"
#   "Rhync_ridleyi.asm.hic.hap1.p_ctg.FINAL.chr.fasta"
#   "Rhync_riparia_5519C.hap1.chr.fasta"
#   "Rhync_tenuis_ref.hap1.chr.fasta"
#   "Rhync_watsonii_6179B.asm.hic.p_ctg.FINAL_chr.fasta"
#   "Rhynchospora_gaudichaudii_6228B.asm.hic.p_ctg.FINAL.chr.fasta"
#   "Rradicans.016.FINAL.chr.fasta"
#   "Rrugosa.chr.fasta"
#   "Rtenerrima.chr.fasta"

islands$genome_class <- ""

islands$genome_class[grep("alba", islands$genome, ignore.case = TRUE)] = "alba"
islands$genome_class[grep("austrobrasiliensis", islands$genome, ignore.case = TRUE)] = "austrobrasiliensis"
islands$genome_class[grep("brevi", islands$genome, ignore.case = TRUE)] = "breviuscula"
islands$genome_class[grep("barbata", islands$genome, ignore.case = TRUE)] = "barbata"
islands$genome_class[grep("cephalotes", islands$genome, ignore.case = TRUE)] = "cephalotes"
islands$genome_class[grep("corymbosa", islands$genome, ignore.case = TRUE)] = "corymbosa"
islands$genome_class[grep("ciliata", islands$genome, ignore.case = TRUE)] = "ciliata"
islands$genome_class[grep("colorata", islands$genome, ignore.case = TRUE)] = "colorata"
islands$genome_class[grep("filiformis", islands$genome, ignore.case = TRUE)] = "filiformis"
islands$genome_class[grep("gaudichaudii", islands$genome, ignore.case = TRUE)] = "gaudichaudii"
islands$genome_class[grep("holoschoenoides", islands$genome, ignore.case = TRUE)] = "holoschoenoides"
islands$genome_class[grep("nervosa", islands$genome, ignore.case = TRUE)] = "nervosa"
islands$genome_class[grep("pubera", islands$genome, ignore.case = TRUE)] = "pubera"
islands$genome_class[grep("ridleyi", islands$genome, ignore.case = TRUE)] = "ridleyi"
islands$genome_class[grep("rugosa", islands$genome, ignore.case = TRUE)] = "rugosa"
islands$genome_class[grep("radicans", islands$genome, ignore.case = TRUE)] = "radicans"
islands$genome_class[grep("riparia", islands$genome, ignore.case = TRUE)] = "riparia"
islands$genome_class[grep("tenuis", islands$genome, ignore.case = TRUE)] = "tenuis"
islands$genome_class[grep("tenerrima", islands$genome, ignore.case = TRUE)] = "tenerrima"
islands$genome_class[grep("watsonii", islands$genome, ignore.case = TRUE)] = "watsonii"

islands_backup = islands

order <- c("rugosa", "alba", "gaudichaudii", "cephalotes", "ridleyi",
           'watsonii', 'radicans', 'pubera', 'breviuscula', 'nervosa', 
           'ciliata', 'colorata', 'tenerrima', 'filiformis', 'austrobrasiliensis', 
           'tenuis', 'riparia', 'barbata', 'holoschoenoides', 'corymbosa')

islands$genome_class <- factor(islands$genome_class, levels = order)

setwd("/home/pwlodzimierz/Rhynchospora/upload_data/16_scatter_plots")


islands$chromosome_size <- 0

chr_sizes_table <- read.csv(file = "/home/pwlodzimierz/Rhynchospora/metadata/chr.no.and.sizes.full.csv")

for(i in 1 : length(data_directories)) {
  print(i)
  assembly <- strsplit(strsplit(data_directories[i], split = "/TRASH_2025/")[[1]][2], split = "/")[[1]][1]
  
  assembly2 <- strsplit(assembly, split = ".fa")[[1]][1]
  
  sizes_genome <- chr_sizes_table[chr_sizes_table$assembly.name == assembly,]
  chromosomes <- unique(sizes_genome$chromosome.name)
  for(j in seq_along(chromosomes)) {
    islands$chromosome_size[islands$genome == assembly2 & islands$chromosome == chromosomes[j]] <- sizes_genome$size[j]
  }
}

write.csv(x = islands, file = "islands_with_chr_size.csv", row.names = FALSE)
# islands <- read.csv(file = "islands_with_chr_size.csv")


### summarise genomes ###

chromosomes_summary <- data.frame(genome = vector(mode = "character", length = 0),
                                  chromosome = vector(mode = "character", length = 0),
                                  no_of_islands = vector(mode = "numeric", length = 0),
                                  total_bp_of_islands = vector(mode = "numeric", length = 0),
                                  chromosome_size = vector(mode = "numeric", length = 0),
                                  mean_island_size = vector(mode = "numeric", length = 0),
                                  mean_spacing_size = vector(mode = "numeric", length = 0),
                                  total_spacing_size = vector(mode = "numeric", length = 0),
                                  mean_island_ED = vector(mode = "numeric", length = 0),
                                  tyba_total_bp = vector(mode = "numeric", length = 0),
                                  tyba_total_no = vector(mode = "numeric", length = 0),
                                  genome_class = vector(mode = "character", length = 0),
                                  exons_bp = vector(mode = "character", length = 0), # TODO
                                  TE_bp = vector(mode = "character", length = 0), # TODO
                                  LTR_bp = vector(mode = "character", length = 0), # TODO
                                  TE_in_islands_bp = vector(mode = "character", length = 0), # TODO
                                  LTR_in_islands_bp = vector(mode = "character", length = 0)) # TODO

genomes <- unique(islands$genome)

for(i in seq_along(genomes)) {
  print(i)
  islands_genome <- islands[islands$genome == genomes[i],]
  
  chromosomes <- unique(islands_genome$chromosome)
  
  has_genes <- FALSE
  path <- grep(genomes[i], data_directories)
  genes <- list.files(path = data_directories[path], pattern = "helixer_filtered.csv", full.names = T)
  if(length(genes) == 1) {
    genes <- read.csv(file = genes)
    genes <- genes[genes$V3 == "exon",]
    has_genes <- TRUE
  }
  
  has_edta <- FALSE
  path <- grep(genomes[i], data_directories)
  edta <- list.files(path = data_directories[path], pattern = "edta_filtered.csv", full.names = T)
  if(length(edta) == 1) {
    edta <- read.csv(file = edta)
    has_edta <- TRUE
  }
  
  
  for(j in seq_along(chromosomes)) {
    exons_bp = 0
    TE_bp = 0
    LTR_bp = 0
    TE_in_islands_bp = 0
    LTR_in_islands_bp = 0
    
    if(has_genes) {
      genes_chr <- genes[genes$V1 == chromosomes[j], ]
      gr <- GRanges(seqnames = genes_chr$V1, 
                    ranges = IRanges(start = genes_chr$V4, end = genes_chr$V5))
      exons_bp <- sum(width(reduce(gr)))
    } 
    if(has_edta) {
      edta_chr <- edta[edta$V1 == chromosomes[j], ]
      gr <- GRanges(seqnames = edta_chr$V1, 
                    ranges = IRanges(start = edta_chr$V4, end = edta_chr$V5))
      TE_bp <- sum(width(reduce(gr)))
      
      LTR_bp = 0 # TODO
      TE_in_islands_bp = 0 # TODO
      LTR_in_islands_bp = 0 # TODO
    }
    
    islands_chromosome <- islands_genome[islands_genome$chromosome == chromosomes[j],]
    
    spacings = c(islands_chromosome$start, islands_chromosome$chromosome_size[1]) - c(1, islands_chromosome$end)
    
    chromosomes_summary <- rbind(chromosomes_summary, data.frame(
      genome = genomes[i],
      chromosome = chromosomes[j],
      no_of_islands = nrow(islands_chromosome),
      total_bp_of_islands = sum(islands_chromosome$total_bp),
      chromosome_size = islands_chromosome$chromosome_size[1],
      mean_island_size = mean(islands_chromosome$total_bp),
      mean_spacing_size = mean(spacings),
      total_spacing_size = sum(spacings),
      mean_island_ED = sum(islands_chromosome$tyba_internal_ED_to_consensus_size_normalised * islands_chromosome$tyba_no) / sum(islands_chromosome$tyba_no),
      tyba_total_bp = sum(islands_chromosome$tyba_total_bp),
      tyba_total_no = sum(islands_chromosome$tyba_no),
      genome_class = islands_chromosome$genome_class[1],
      exons_bp = exons_bp,
      TE_bp = TE_bp,
      LTR_bp = LTR_bp,
      TE_in_islands_bp = TE_in_islands_bp,
      LTR_in_islands_bp = LTR_in_islands_bp))
  }
  
}
chromosomes_summary$island_fraction = 100 * chromosomes_summary$total_bp_of_islands / chromosomes_summary$chromosome_size

write.csv(x = chromosomes_summary, file = "chromosomes_summary.csv", row.names = FALSE)
chromosomes_summary <- read.csv("chromosomes_summary.csv")

# remove austro extra haplotypes

chromosomes_summary <- chromosomes_summary[-c(49,50,52,53,55,56),]


### ABSTRACT figures plots

order <- c("rugosa", "alba", "gaudichaudii","cephalotes", 
           "ridleyi", "watsonii",
           "radicans",  "pubera", "breviuscula", 
           "nervosa", "ciliata", "colorata","tenerrima", "filiformis", "austrobrasiliensis",
           "tenuis", "riparia", "barbata",
           "holoschoenoides", "corymbosa")


### scatters vs chromosome size
colors_vector <- c("#E69F00bb", 
                   "#D55E00bb", 
                   "#D95F02bb", 
                   "#FC8D62bb", 
                   "#E5C494bb", 
                   "#F0E442bb", 
                   "#FFD92Fbb", 
                   "#A6D854bb", 
                   "#66C2A5bb", 
                   "#009E73bb", 
                   "#1B9E77bb", 
                   "#66A61Ebb", 
                   "#CC79A7bb", 
                   "#E78AC3bb", 
                   "#E7298Abb", 
                   "#B3B3CCbb", 
                   "#8DA0CBbb", 
                   "#7570B3bb", 
                   "#56B4E9bb", 
                   "#0072B2bb")

# Define point shapes (pch values) to improve distinction
pch_vector <- c(15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18,
                15, 16, 17, 18)



# randomise order for more even plotting

chromosomes_summary <- chromosomes_summary[sample(nrow(chromosomes_summary)),]
# Ensure factor levels are used correctly

genome_factor <- factor(chromosomes_summary$genome_class)
genome_indices <- NULL
for(ki in 1 : nrow(chromosomes_summary)) {
  
  genome_indices <- c(genome_indices, which(order == chromosomes_summary$genome_class[ki]))
  
}


################
# Correlations #
################

# get the phylogeny

setwd("/home/pwlodzimierz/Rhynchospora/git_Rhynchospora")

tree_gene_based <- ape::read.tree(file = "SpeciesTree_rooted_node_labels.newick")

library(ape)
library(nlme)

setwd("/home/pwlodzimierz/Rhynchospora/upload_data/16_scatter_plots")

chromosomes_summary$species <- ""
for(i in 1 : nrow(chromosomes_summary)) {
  chromosomes_summary$species[i] <- paste0("R", chromosomes_summary$genome_class[i])
}

chromosomes_summary$chromosome_size_mbp = chromosomes_summary$chromosome_size / 1000000
chromosomes_summary$exons_ratio = 100 * chromosomes_summary$exons_bp / chromosomes_summary$chromosome_size
chromosomes_summary$tyba_total_ratio = 100 * chromosomes_summary$tyba_total_bp / chromosomes_summary$chromosome_size
chromosomes_summary$mean_spacing_kbp <- chromosomes_summary$mean_spacing_size / 1000
chromosomes_summary$mean_island_kbp = chromosomes_summary$mean_island_size / 1000
chromosomes_summary$TE_ratio <- 100 * chromosomes_summary$TE_bp / chromosomes_summary$chromosome_size

model_gene <- gls(exons_ratio ~ chromosome_size_mbp,
                  correlation = corBrownian(phy = tree_gene_based, form = ~species),
                  data = chromosomes_summary ) # it doesn't matter in the end, cannot correct chromosomes with a species tree
model_te <- gls(TE_ratio ~ chromosome_size_mbp,
                correlation = corBrownian(phy = tree_gene_based, form = ~species),
                data = chromosomes_summary )
model_tyba <- gls(tyba_total_ratio ~ chromosome_size_mbp,
              correlation = corBrownian(phy = tree_gene_based, form = ~species),
              data = chromosomes_summary[chromosomes_summary$genome_class != "filiformis",] )
model_spacing <- gls(mean_spacing_kbp ~ chromosome_size_mbp,
              correlation = corBrownian(phy = tree_gene_based, form = ~species),
              data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),] )

summary(gls(mean_spacing_kbp ~ chromosome_size_mbp,
            correlation = corBrownian(phy = tree_gene_based, form = ~species),
            data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata", "tenuis", "pubera")),] ))

model_islandbp <- gls(mean_island_kbp ~ chromosome_size_mbp,
              correlation = corBrownian(phy = tree_gene_based, form = ~species),
              data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),] )
model_islandno <- gls(no_of_islands ~ chromosome_size_mbp,
              correlation = corBrownian(phy = tree_gene_based, form = ~species),
              data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),] )

p_values <- c(
  exons_ratio = summary(model_gene)$tTable["chromosome_size_mbp", "p-value"],
  TE_ratio = summary(model_te)$tTable["chromosome_size_mbp", "p-value"],
  tyba_total_ratio = summary(model_tyba)$tTable["chromosome_size_mbp", "p-value"],
  mean_spacing_kbp = summary(model_spacing)$tTable["chromosome_size_mbp", "p-value"],
  mean_island_kbp = summary(model_islandbp)$tTable["chromosome_size_mbp", "p-value"],
  no_of_islands = summary(model_islandno)$tTable["chromosome_size_mbp", "p-value"]
)


p_values2 <- c(
  exons_ratio = summary(lm(exons_ratio ~ chromosome_size_mbp, data = chromosomes_summary))$coefficients[2, 4],
  TE_ratio = summary(lm(TE_ratio ~ chromosome_size_mbp, data = chromosomes_summary))$coefficients[2, 4],
  tyba_total_ratio = summary(lm(tyba_total_ratio ~ chromosome_size_mbp, data = chromosomes_summary[chromosomes_summary$genome_class != "filiformis",]))$coefficients[2, 4],
  mean_spacing_kbp = summary(lm(mean_spacing_kbp ~ chromosome_size_mbp, data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),]))$coefficients[2, 4],
  mean_island_kbp = summary(lm(mean_island_kbp ~ chromosome_size_mbp, data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),]))$coefficients[2, 4],
  no_of_islands = summary(lm(no_of_islands ~ chromosome_size_mbp, data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),]))$coefficients[2, 4]
)
options(scipen = 999)
round(p_values2, 5)
round(p_values, 5) # it doesn't matter in the end, cannot correct chromosomes with a species tree


if(F) {

  ### collapsed

  chromosomes_summary_collapsed <- data.frame()

  for(i in 1 : length(unique(chromosomes_summary$species))) {
    chromosomes_summary_t <- chromosomes_summary[chromosomes_summary$species == unique(chromosomes_summary$species)[i],]
    chromosomes_summary_collapsed <- rbind(chromosomes_summary_collapsed, data.frame(chromosome_size_mbp = mean(chromosomes_summary_t$chromosome_size_mbp),
                                                                                     exons_ratio = mean(chromosomes_summary_t$exons_ratio),
                                                                                     TE_ratio = mean(chromosomes_summary_t$TE_ratio),
                                                                                     tyba_total_ratio = mean(chromosomes_summary_t$tyba_total_ratio),
                                                                                     mean_spacing_kbp = mean(chromosomes_summary_t$mean_spacing_kbp),
                                                                                     mean_island_kbp = mean(chromosomes_summary_t$mean_island_kbp),
                                                                                     no_of_islands = mean(chromosomes_summary_t$no_of_islands),
                                                                                     species = unique(chromosomes_summary$species)[i]))
  }

  model_gene <- gls(exons_ratio ~ chromosome_size_mbp,
                    correlation = corBrownian(phy = tree_gene_based, form = ~species),
                    data = chromosomes_summary_collapsed )
  model_te <- gls(TE_ratio ~ chromosome_size_mbp,
                  correlation = corBrownian(phy = tree_gene_based, form = ~species),
                  data = chromosomes_summary_collapsed )
  model_tyba <- gls(tyba_total_ratio ~ chromosome_size_mbp,
                    correlation = corBrownian(phy = tree_gene_based, form = ~species),
                    data = chromosomes_summary_collapsed[chromosomes_summary_collapsed$species != "Rfiliformis",] )
  model_spacing <- gls(mean_spacing_kbp ~ chromosome_size_mbp,
                       correlation = corBrownian(phy = tree_gene_based, form = ~species),
                       data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("Rfiliformis", "Rbarbata")),] )
  model_islandbp <- gls(mean_island_kbp ~ chromosome_size_mbp,
                        correlation = corBrownian(phy = tree_gene_based, form = ~species),
                        data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("Rfiliformis", "Rbarbata")),] )
  model_islandno <- gls(no_of_islands ~ chromosome_size_mbp,
                        correlation = corBrownian(phy = tree_gene_based, form = ~species),
                        data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("Rfiliformis", "Rbarbata")),] )

  p_values <- c(
    exons_ratio = summary(model_gene)$tTable["chromosome_size_mbp", "p-value"],
    TE_ratio = summary(model_te)$tTable["chromosome_size_mbp", "p-value"],
    tyba_total_ratio = summary(model_tyba)$tTable["chromosome_size_mbp", "p-value"],
    mean_spacing_kbp = summary(model_spacing)$tTable["chromosome_size_mbp", "p-value"],
    mean_island_kbp = summary(model_islandbp)$tTable["chromosome_size_mbp", "p-value"],
    no_of_islands = summary(model_islandno)$tTable["chromosome_size_mbp", "p-value"]
  )


  p_values2 <- c(
    exons_ratio = summary(lm(exons_ratio ~ chromosome_size_mbp, data = chromosomes_summary_collapsed))$coefficients[2, 4],
    TE_ratio = summary(lm(TE_ratio ~ chromosome_size_mbp, data = chromosomes_summary_collapsed))$coefficients[2, 4],
    tyba_total_ratio = summary(lm(tyba_total_ratio ~ chromosome_size_mbp, data = chromosomes_summary_collapsed[chromosomes_summary_collapsed$species != "filiformis",]))$coefficients[2, 4],
    mean_spacing_kbp = summary(lm(mean_spacing_kbp ~ chromosome_size_mbp, data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("filiformis", "barbata")),]))$coefficients[2, 4],
    mean_island_kbp = summary(lm(mean_island_kbp ~ chromosome_size_mbp, data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("filiformis", "barbata")),]))$coefficients[2, 4],
    no_of_islands = summary(lm(no_of_islands ~ chromosome_size_mbp, data = chromosomes_summary_collapsed[!(chromosomes_summary_collapsed$species %in% c("filiformis", "barbata")),]))$coefficients[2, 4]
  )
  round(p_values2, 5)
  round(p_values, 5)
  
}




### x = chromosome size y = exons %
pdf(file = "x = chromosome size y = exons perc, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp, 
     chromosomes_summary$exons_ratio, 
     pch = pch_vector[genome_indices] , 
     col = colors_vector[genome_indices] ,
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Genes (exons) %", cex = 2)
# abline(model_gene, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(exons_ratio ~ chromosome_size_mbp, data = chromosomes_summary ), col = "black", lwd = 2, lty = 1)
dev.off()



### x = chromosome size y = TEs %
pdf(file = "x = chromosome size y = TEs perc, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp, 
     chromosomes_summary$TE_ratio, 
     pch = pch_vector[genome_indices], 
     col = colors_vector[genome_indices],
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Transposons %", cex = 2)
# abline(model_te, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(TE_ratio ~ chromosome_size_mbp, data = chromosomes_summary[chromosomes_summary$TE_bp > 0,] ), col = "black", lwd = 2, lty = 1)
dev.off()


### x = chromosome size y = repeat %
pdf(file = "x = chromosome size y = tyba perc, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp, 
     chromosomes_summary$tyba_total_ratio, 
     pch = pch_vector[genome_indices], 
     col = colors_vector[genome_indices],
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Tyba %", cex = 2)
# abline(model_tyba, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(tyba_total_ratio ~ chromosome_size_mbp, data = chromosomes_summary[chromosomes_summary$tyba_total_bp > 0,] ), col = "black", lwd = 2, lty = 1)
dev.off()

### x = chromosome size y = mean spacing without barbata
pdf(file = "x = chromosome size y = spacing Kbp, no barbata, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata"))], 
     chromosomes_summary$mean_spacing_kbp[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata"))], 
     pch = pch_vector[genome_indices][!(chromosomes_summary$genome_class %in% c("filiformis", "barbata"))], 
     col = colors_vector[genome_indices][!(chromosomes_summary$genome_class %in% c("filiformis", "barbata"))],
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Mean spacing size, Kbp", cex = 2)
# abline(model_spacing, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(mean_spacing_kbp ~ chromosome_size_mbp, data = chromosomes_summary[!(chromosomes_summary$genome_class %in% c("filiformis", "barbata")),] ), col = "black", lwd = 2, lty = 1)
dev.off()

### x = chromosome size y = mean island
pdf(file = "x = chromosome size y = island Kbp, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp[chromosomes_summary$genome_class != "filiformis"], 
     chromosomes_summary$mean_island_kbp[chromosomes_summary$genome_class != "filiformis"], 
     pch = pch_vector[genome_indices][chromosomes_summary$genome_class != "filiformis"], 
     col = colors_vector[genome_indices][chromosomes_summary$genome_class != "filiformis"],
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Mean centromere array size, Kbp", cex = 2)
# abline(model_islandbp, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(mean_island_kbp ~ chromosome_size_mbp, data = chromosomes_summary[chromosomes_summary$genome_class != "filiformis",] ), col = "black", lwd = 2, lty = 1)
dev.off()


### x = chromosome size y = islands number
pdf(file = "x = chromosome size y = islands number, chromosome scatters.pdf", width = 7, height = 7, onefile = TRUE)
par(mar = c(4,4,2,2), oma = c(1,1,1,1))
plot(chromosomes_summary$chromosome_size_mbp[chromosomes_summary$genome_class != "filiformis"], 
     chromosomes_summary$no_of_islands[chromosomes_summary$genome_class != "filiformis"], 
     pch = pch_vector[genome_indices][chromosomes_summary$genome_class != "filiformis"], 
     col = colors_vector[genome_indices][chromosomes_summary$genome_class != "filiformis"],
     xlim = c(0,400),
     xlab = "Chromosome size, Mbp",
     ylab = "Number of centromeric arrays", cex = 2)
# abline(model_islandno, col = "gray", lwd = 2, lty = 1)# it doesn't matter in the end, cannot correct chromosomes with a species tree
abline(lm(no_of_islands ~ chromosome_size_mbp, data = chromosomes_summary[chromosomes_summary$genome_class != "filiformis",] ), col = "black", lwd = 2, lty = 1)
dev.off()




