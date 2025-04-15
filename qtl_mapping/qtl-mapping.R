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
                         phefile="~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping/costus_pheno_rqtl_2025Jan24.csv", 
                         estimate.map=FALSE, genotypes=c("AA","AB","BB"))
summary(cross.data)

## Rename the chromosomes based on their dominant mappings to C. lasius genome
## But remember, chr. 4 is a composite with 2 and 9, and chr. 7 is a composite with 5
new_names <- c("2", "4", "5", "7", "3", "6", "8", "1", "9")
names(cross.data$geno)[1:length(new_names)] <- new_names

## Check that linkage groups were renamed appropriately
chrnames(cross.data)

# Reorder the geno list by sorting the names numerically
cross.data$geno <- cross.data$geno[order(as.numeric(names(cross.data$geno)))]

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
# Step 4. Define functions for QTL mapping
# ========================================================================

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 1. PERFORM QTL SCAN
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 2. PRINT SUMMARY OF SCAN
qtl_summary <- function(scan_result, perm_data, alpha_levels = c(0.05, 0.1, 0.15)) {
  summaries <- lapply(alpha_levels, function(alpha) {
    summary(scan_result, perms = perm_data, alpha = alpha, pvalues = TRUE)
  })
  names(summaries) <- paste0("alpha_", alpha_levels)
  return(summaries)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 3. PLOT QTL PEAKS
plot_qtl_peaks <- function(scan_result, cutoff_lod, chromosomes, positions, trait_name) {
  
  pdf_width <- max(c(length(chromosomes) * 3, 7))
  
  # Save plot to PDF:
  pdf(paste0(trait_name, "_qtl_peaks.pdf"), width = pdf_width)
  par(mfrow = c(1, length(chromosomes)))
  for (i in 1:length(chromosomes)) {
    plot(scan_result, chr = chromosomes[i], main = paste0("Chr", chromosomes[i], " - LOD cutoff ", round(cutoff_lod, 2)))
    abline(h = cutoff_lod, col = "red", lwd = 3)
    abline(v = positions[i], col = "green", lwd = 3)
  }
  dev.off()
  
  # Plot to screen:
  par(mfrow = c(1, length(chromosomes)))
  for (i in 1:length(chromosomes)) {
    plot(scan_result, chr = chromosomes[i], main = paste0("Chr", chromosomes[i], " - LOD cutoff ", round(cutoff_lod, 2)))
    abline(h = cutoff_lod, col = "red", lwd = 3)
    abline(v = positions[i], col = "green", lwd = 3)
  }
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 4. EXTRACT THE POSITIONS OF SIGNIFICANT PEAKS
save_peak_positions <- function(scan_result, cutoff_lod_0.1, cutoff_lod_0.05, trait_name) {
  write.table(scan_result[scan_result$lod > cutoff_lod_0.1,], 
              file = paste0(trait_name, "_scanone-hk_peaks-above-0.1LODcutoff.tsv"), sep = "\t", quote = FALSE)
  
  write.table(scan_result[scan_result$lod > cutoff_lod_0.05,], 
              file = paste0(trait_name, "_scanone-hk_peaks-above-0.05LODcutoff.tsv"), sep = "\t", quote = FALSE)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 5. PLOT QTL EFFECTS
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 6. CALCULATE LOD AND BAYESIAN CREDIBLE INTERVALS
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 7. PREPARE DATA TO CREATE MANHATTAN PLOT
mh_plot_multi <-function(chr_pos_lod, cutoff_lod){
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

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 8. GENERATE MANHATTAN PLOTS
plot_manhattan <- function(scanone_results, cutoff_lod, trait_name) {
  ## load the qqman manhattan plot function
  source("~/Dropbox/Costus/genetic_mapping/R:qtl/adapted_qqman.R")
  pdf(paste0(trait_name, "_manhattan_plot.pdf"), width = 10, height = 6)
  mh_plot_multi(scanone_results, cutoff_lod)
  title(main = paste0(trait_name, " QTL: QTL map"))
  dev.off()
  
  mh_plot_multi(scanone_results, cutoff_lod)
  title(main = paste0(trait_name, " QTL: QTL map"))
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 9. REFINE THE QTL POSITIONS AND FIT THE MODEL WITH INTERACTIONS
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

# ========================================================================
# Step 5. QTL mapping with scanone
# ========================================================================

# Load data
F1parent <- pull.pheno(data_prob, "F1parent")
perm_file <- "scanone_1000perm.rds"

for (trait_name in names(data_prob$pheno)[25:50]) {

  print(paste("Now running analyses for", trait_name))
  
  # QTL scan and save results
  scan_result <- perform_qtl_scan(data_prob, trait_name, trait_name, F1parent, 
                                  perm_file)
  
  # Summarize the peaks at different alpha levels
  qtl_peaks <- qtl_summary(scan_result$scan_result, readRDS(perm_file))
  print(qtl_peaks)
  
  # Skip to next trait if there aren't any significant LOD peaks
  if (nrow(qtl_peaks$alpha_0.1) == 0) {
    print(paste("There were no LOD peaks above the threshold for", trait_name))
    next
  }
  
  # Plot QTL peaks
  plot_qtl_peaks(scan_result$scan_result, scan_result$cutoff_lod_0.1, 
                 qtl_peaks$alpha_0.1$chr, qtl_peaks$alpha_0.1$pos, trait_name)
  
  # Extract the positions of significant peaks & save them based on LOD thresholds
  save_peak_positions(scan_result$scan_result, scan_result$cutoff_lod_0.1, 
                      scan_result$cutoff_lod_0.05, trait_name)
  
  # Calculate LOD intervals
  lod_intervals <- calc_lod_intervals(scan_result$scan_result, chromosomes = 
                                        qtl_peaks$alpha_0.1$chr)
  print(lod_intervals)
  
  # Plot QTL effects
  plot_qtl_effects(data_prob, chromosomes = qtl_peaks$alpha_0.1$chr, 
                   positions = qtl_peaks$alpha_0.1$pos, trait_name = trait_name)
  
  # Pull out the genotype probabilities at the nearest pseudomarkers
  qtl <- makeqtl(data_prob, chr = qtl_peaks$alpha_0.1$chr, 
                 pos = qtl_peaks$alpha_0.1$pos, what = "prob")
  
  # Fit a multiple-QTL model using Haley-Knott regression (method = "hk")
  if (length(qtl$name) == 1) {formula = "y ~ Q1 + F1parent"}
  if (length(qtl$name) == 2) {formula = "y ~ Q1 + Q2 + F1parent"}
  if (length(qtl$name) == 3) {formula = "y ~ Q1 + Q2 + Q3 + F1parent"}
  if (length(qtl$name) == 4) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + F1parent"}
  if (length(qtl$name) == 5) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + F1parent"}
  if (length(qtl$name) == 6) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + F1parent"}
  if (length(qtl$name) == 7) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + F1parent"}
  if (length(qtl$name) == 8) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + F1parent"}
  
  out.fq <- fitqtl(data_prob, pheno.col=trait_name, qtl=qtl, method="hk", 
                   get.ests=TRUE, covar=as.data.frame(F1parent), 
                   formula = formula)
  print(summary(out.fq))
  
  # Determine whether there is an interaction between the two QTL by fitting the 
  # model with the interaction
  #out.fqi <- fitqtl(data_prob, pheno.col=trait_name, qtl=qtl, method="hk", 
  #                  get.ests=TRUE, covar=as.data.frame(F1parent), 
  #                  formula = y ~ Q1 + Q2 + Q3 + Q1:Q2 + Q1:Q3 + Q2:Q3 + F1parent)
  #summary(out.fqi)
  
  # Also assess interaction with addint, which adds one interaction at a time
  if (length(qtl$name) > 1) {addint(data_prob, pheno.col = trait_name, qtl=qtl, method="hk")}
  
  #################
  # Refine the qtl positions
  # each QTL is moved to the position giving the highest likelihood,
  # and the entire process is repeated until no further improvement in likelihood
  # can be obtained
  rqtl <- refineqtl(data_prob, pheno.col = trait_name, qtl=qtl, method="hk")
  rqtl
  
  if (length(rqtl$name) == 1) {formula = "y ~ Q1 + Q2 + F1parent"}
  if (length(rqtl$name) == 2) {formula = "y ~ Q1 + Q2 + Q3 + F1parent"}
  if (length(rqtl$name) == 3) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + F1parent"}
  if (length(rqtl$name) == 4) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + Q5 + F1parent"}
  
  # Look for additional qtl, using the refined qtl positions
  out.aq <- addqtl(data_prob, pheno.col = trait_name, qtl=rqtl, method="hk", 
                   covar=as.data.frame(F1parent), 
                   formula = formula)
  
  # Additional QTL with LOD that exceed 0.1 cutoff
  add.qtl <- max(out.aq)
  print(add.qtl)
  
  # Pull out the genotype probabilities at the nearest pseudomarkers for the 
  # additional QTL
  qtl <- makeqtl(data_prob, chr = as.numeric(c(rqtl$chr, as.character(add.qtl$chr))), 
                 pos = as.numeric(c(rqtl$pos, as.character(add.qtl$pos))), what = "prob")
  
  if (length(qtl$name) == 1) {formula = "y ~ Q1 + F1parent"}
  if (length(qtl$name) == 2) {formula = "y ~ Q1 + Q2 + F1parent"}
  if (length(qtl$name) == 3) {formula = "y ~ Q1 + Q2 + Q3 + F1parent"}
  if (length(qtl$name) == 4) {formula = "y ~ Q1 + Q2 + Q3 + Q4 + F1parent"}
  
  # Fit a multiple-QTL model with these new QTL added
  out.fq <- fitqtl(data_prob, pheno.col=trait_name, qtl=qtl, method="hk", 
                   get.ests=TRUE, covar=as.data.frame(F1parent), 
                   formula = formula)
  summary(out.fq)
  
  # Check for interactions between new and old QTL
  addint(data_prob, pheno.col = trait_name, qtl=qtl, method="hk")
  
  # Perform and plot Manhattan
  plot_manhattan(scan_result$scan_result, scan_result$cutoff_lod_0.1, trait_name)
}








