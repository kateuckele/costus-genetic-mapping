## Calculate correlations among traits
## K. Uckele Nov. 28, 2023

##libraries
library(GGally)
library(ggcorrplot)
library(dplyr)
library(readr)


## read in phenotypic data
pheno <- read_csv(
  "~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/processed_data/costus_pheno_rqtl_2025Jan24.csv",
  col_names = TRUE,
  na = c("", "NA", "-")  # Specify additional strings to be treated as NA
)

## separate traits into two distinct sets  
morphometric <- pheno[,2:26]
color <- pheno[,27:49]

## change column order
#pheno <- pheno %>% select(ANL, ANW, CLL, COL, COLL, LABL, LABW, RALA, RAST,
#                          STAE, STAL, STATL, STAW, STIW, STYL, TUA, CAL, INFA,
#                          VNG, EFNSC, FNSC, VEFN, VFN)

##Correlations
cor(pheno, method = "pearson", use = "pairwise.complete.obs")

##Check correlations (as scatterplots), distribution and print correlation coefficient 
ggpairs(pheno, title="correlogram with ggpairs()") 

################################################################################
# Compute a correlation matrix for morphometric traits
corr <- round(cor(morphometric, method = "pearson", use = "pairwise.complete.obs"), 1)
head(corr)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(morphometric)
p.mat[lower.tri(p.mat)] <- NA
p.mat.corrected <- apply(p.mat, 2, function(x) p.adjust(x, method = "holm"))
p.mat.corrected[lower.tri(p.mat.corrected)] <- p.mat.corrected[upper.tri(p.mat.corrected)]

# Visualize the correlation matrix
# --------------------------------
# method = "square" (default)
ggcorrplot(corr, type = "lower", p.mat = p.mat.corrected, lab = TRUE, insig = "blank")


# Reordering the correlation matrix
# --------------------------------
# using hierarchical clustering
pdf("~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/figures/morphometric_traits_corr_plot_p_corrected.pdf", height = 10, width = 10)
ggcorrplot(corr, hc.order = TRUE, insig = "blank", 
           type = "lower", p.mat = p.mat.corrected, lab = TRUE)
dev.off()

################################################################################
# Compute a correlation matrix for color traits
corr <- round(cor(color, method = "pearson", use = "pairwise.complete.obs"), 1)
head(corr)

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(color)
p.mat[lower.tri(p.mat)] <- NA
p.mat.corrected <- apply(p.mat, 2, function(x) p.adjust(x, method = "holm"))
p.mat.corrected[lower.tri(p.mat.corrected)] <- p.mat.corrected[upper.tri(p.mat.corrected)]

# Visualize the correlation matrix
# --------------------------------
# method = "square" (default)
ggcorrplot(corr, type = "lower", p.mat = p.mat.corrected, lab = TRUE, insig = "blank")


# Reordering the correlation matrix
# --------------------------------
# using hierarchical clustering
pdf("~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/figures/colormetric_traits_corr_plot_p_corrected.pdf", height = 10, width = 10)
ggcorrplot(corr, hc.order = FALSE, insig = "blank", 
           type = "lower", p.mat = p.mat.corrected, lab = TRUE)
dev.off()
