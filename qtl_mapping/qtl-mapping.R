## clear workspace
rm(list = ls())

## Load libraries
library(qtl)
library(snow)
library(formatR)
#library(qtl2)

# ========================================================================
# Step 1. Read in cross data
# ========================================================================

cross.data <- read.cross("csvs", genfile="~/Dropbox/Costus/costus-genetic-mapping/linkage_map/mapthis_LG.csv", 
                         phefile="~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping/costus_pheno_2024Aug16.csv", 
                         estimate.map=FALSE, genotypes=c("AA","AB","BB"))
summary(cross.data)

# ========================================================================
# Step 2. Calculate genotype probabilities
# ========================================================================
## Calculate conditional genotype probabilities, conditional on the available marker data
## step = step size (in cM) at which the probabilities are calculated
## error.prob = Assumed genotyping error rate
data_prob <- calc.genoprob(cross.data, step = 2, error.prob=0)

# ========================================================================
# Step 3. Run permutations to get genome-wide LOD significance thresholds
# ========================================================================
## Run scanone permutations 
data_prob_perm <- data_prob
data_prob_perm$pheno <- data_prob_perm$pheno[,2:24]
operm.hk <- scanone(data_prob_perm, method="hk", n.perm=1000, n.cluster = 6)
saveRDS(operm.hk, "~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping/scanone_1000perm.rds")
summary(operm.hk)
# LOD thresholds (1002 permutations)
# lod
# 5%  3.72
# 10% 3.31

# ========================================================================
# Step 4. QTL mapping with scanone
# ========================================================================
## Define functions
## Perform qtl scan
perform_qtl_scan <- function(data_prob, trait_col, trait_name, f1_covar, perm_file) {
  # QTL scan using Haley-Knott regression
  scan_result <- scanone(data_prob, method = "hk", pheno.col = trait_col, addcovar = f1_covar)
  
  # Load permutation data
  perm_data <- readRDS(perm_file)
  cutoff_lod_0.1 <- summary(perm_data)[c("10%"),]
  cutoff_lod_0.05 <- summary(perm_data)[c("5%"),]
  
  # Save results
  write.table(scan_result, file = paste0(trait_name, "-F1parent.scanone.hk.tsv"), sep = "\t", quote = FALSE)
  saveRDS(scan_result, file = paste0(trait_name, "-F1parent.scanone.hk.rds"))
  
  return(list(scan_result = scan_result, cutoff_lod_0.1 = cutoff_lod_0.1, cutoff_lod_0.05 = cutoff_lod_0.05))
}

## Print summary of scan
qtl_summary <- function(scan_result, perm_data, alpha_levels = c(0.05, 0.1, 0.15)) {
  summaries <- lapply(alpha_levels, function(alpha) {
    summary(scan_result, perms = perm_data, alpha = alpha, pvalues = TRUE)
  })
  names(summaries) <- paste0("alpha_", alpha_levels)
  return(summaries)
}

## Plot QTL peaks
plot_qtl_peaks <- function(scan_result, cutoff_lod, chromosomes, positions, trait_name) {
  # Save plot to PDF:
  pdf(paste0(trait_name, "_qtl_peaks.pdf"))
  par(mfrow = c(1, length(chromosomes)))
  for (i in 1:length(chromosomes)) {
    plot(scan_result, chr = chromosomes[i], main = paste0("Chr", chromosomes[i], " - LOD cutoff ", round(cutoff_lod, 2)))
    abline(h = cutoff_lod, col = "red", lwd = 3)
    abline(v = positions[chromosomes[i]], col = "green", lwd = 3)
  }
  dev.off()
  
  # Plot to screen:
  par(mfrow = c(1, length(chromosomes)))
  for (i in 1:length(chromosomes)) {
    plot(scan_result, chr = chromosomes[i], main = paste0("Chr", chromosomes[i], " - LOD cutoff ", round(cutoff_lod, 2)))
    abline(h = cutoff_lod, col = "red", lwd = 3)
    abline(v = positions[chromosomes[i]], col = "green", lwd = 3)
  }
}

## Extract the positions of significant peaks and save them based on LOD thresholds
save_peak_positions <- function(scan_result, cutoff_lod_0.1, cutoff_lod_0.05, trait_name) {
  write.table(scan_result[scan_result$lod > cutoff_lod_0.1,], 
              file = paste0(trait_name, "_scanone-hk_peaks-above-0.1LODcutoff.tsv"), sep = "\t", quote = FALSE)
  
  write.table(scan_result[scan_result$lod > cutoff_lod_0.05,], 
              file = paste0(trait_name, "_scanone-hk_peaks-above-0.05LODcutoff.tsv"), sep = "\t", quote = FALSE)
}

## Plot QTL effects
plot_qtl_effects <- function(data_prob, chromosomes, positions, trait_name) {
  ## Save plot to PDF:
  pdf(paste0(trait_name, "_qtl_effect_plot.pdf"), width = 12)
  par(mfrow = c(1, length(chromosomes)))
  
  for (i in 1:length(chromosomes)) {
    marker <- find.marker(data_prob, chr = chromosomes[i], pos = positions[i])
    data_sim <- sim.geno(data_prob, step=1, n.draws=16)
    effectplot(data_sim, pheno.col = trait_name, mname1 = marker, 
               main = paste0(chromosomes[i], "@", positions[i]))
  }
  
  dev.off()
  
  ## Plot to screen:
  par(mfrow = c(1, length(chromosomes)))
  for (i in 1:length(chromosomes)) {
    marker <- find.marker(data_prob, chr = chromosomes[i], pos = positions[i])
    data_sim <- sim.geno(data_prob, step=1, n.draws=16)
    effectplot(data_sim, pheno.col = trait_name, mname1 = marker, 
               main = paste0(chromosomes[i], "@", positions[i]))
  }
}

## Calculate LOD and Bayesian credible intervals for multiple chromosomes
calc_lod_intervals <- function(scan_result, chromosomes) {
  intervals <- list()
  for (chr in chromosomes) {
    intervals[[paste0("chr", chr)]] <- list(
      lod_1.5 = lodint(scan_result, chr = chr, drop = 1.5, expandtomarkers = TRUE),
      lod_2 = lodint(scan_result, chr = chr, drop = 2, expandtomarkers = TRUE),
      bayes = bayesint(scan_result, chr = chr, prob = 0.95, expandtomarkers = TRUE)
    )
  }
  return(intervals)
}

## Function: prepares data frame with chr,pos,lod columns to create manhattan plot
#  make sure 'cutoff_lod' is set to the LOD threshold
mh_plot_multi <-function(chr_pos_lod){
  cpl.df <- as.data.frame(chr_pos_lod)
  data_trim<-{}
  # if using non-integer chrom names (i.e. not group numbers),
  # use this instead of following line
  # data_trim$CHR<-chr_to_group(as.character.factor(cpl.df$chr))
  data_trim$CHR<-as.numeric(cpl.df$chr)
  data_trim$BP<-(cpl.df$pos)*1e6
  data_trim$P<-cpl.df$lod
  data_trim<-as.data.frame(data_trim)
  data_trim<-na.omit(data_trim)
  # data_trim$CHR <- as.character.factor(data_trim$CHR)
  manhattan_multi(data_trim,
                  genomewideline=cutoff_lod,
                  suggestiveline=-1, 
                  cex.axis=1.5, 
                  cex.lab=4)
  
}

## Generate Manhattan plots based on QTL scan results
plot_manhattan <- function(scanone_results, cutoff_lod, trait_name) {
  ## load the qqman manhattan plot function
  source("~/Dropbox/Costus/genetic_mapping/R:qtl/adapted_qqman.R")
  pdf(paste0(trait_name, "_manhattan_plot.pdf"), width = 10, height = 6)
  mh_plot_multi(scanone_results)
  title(main = paste0(trait_name, " QTL: QTL map"))
  dev.off()
  
  mh_plot_multi(scanone_results)
  title(main = paste0(trait_name, " QTL: QTL map"))
  
}

## Refine the QTL positions and fit the model with interactions
refine_qtl_model <- function(data_prob, trait_col, chromosomes, positions, f1_covar) {
  qtl <- makeqtl(data_prob, chr = chromosomes, pos = positions, what = "prob")
  
  rqtl <- refineqtl(data_prob, pheno.col = trait_col, qtl = qtl, method = "hk")
  refined_positions <- rqtl$pos
  
  # Fit additive model
  out_fq <- fitqtl(data_prob, pheno.col = trait_col, qtl = rqtl, method = "hk", 
                   get.ests = TRUE, covar = f1_covar, formula = y ~ Q1 + Q2 + Q3 + Q4 + F1parent)
  
  # Fit interaction model
  out_fqi <- fitqtl(data_prob, pheno.col = trait_col, qtl = rqtl, method = "hk", 
                    get.ests = TRUE, covar = f1_covar, 
                    formula = y ~ Q1 + Q2 + Q3 + Q4 + Q1:Q2 + Q1:Q3 + Q1:Q4 + Q2:Q3 + Q2:Q4 + Q3:Q4 + F1parent)
  
  return(list(additive_model = out_fq, interaction_model = out_fqi, refined_positions = refined_positions))
}








