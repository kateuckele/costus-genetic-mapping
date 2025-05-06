###############################################################################
## Costus plotting workflow
##  • Reads cross object (R/qtl v1)
##  • Renames linkage groups → C.lasius chromosomes
##  • Computes segregation‑distortion p‑values
##  • Runs scanone (Haley–Knott)
##  • Produces per‑chromosome plots: LOD  +  –log10(p)
###############################################################################

## 0.  House‑keeping ----------------------------------------------------------
rm(list = ls())
library(qtl)

setwd("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping")

# ========================================================================
# Step 1. Read in cross data
# ========================================================================

## 1.  Read cross -------------------------------------------------------------
cross <- read.cross(
  format   = "csvs",
  genfile  = "../linkage_map/mapthis_LG.csv",
  phefile  = "costus_pheno_rqtl_2025Jan24.csv",
  estimate.map = FALSE,
  genotypes    = c("AA","AB","BB")
)


## 2.  Rename linkage groups --------------------------------------------------
# cross$geno is a named list: LG1, LG2, …
new_chr_names <- c("2","4","5","7","3","6","8","1","9")
stopifnot(length(new_chr_names) == length(cross$geno))
names(cross$geno) <- new_chr_names

# Re‑order geno list numerically
cross$geno <- cross$geno[ order(as.numeric(names(cross$geno))) ]
message("Chromosomes after renaming: ", paste(chrnames(cross), collapse = ", "))

## 3.  Segregation‑distortion p‑values ---------------------------------------
gt_tab            <- geno.table(cross, scanone.output=TRUE) # χ² test per marker

## 4.  Genotype probabilities (for completeness) -----------------------------
cross <- calc.genoprob(cross, step = 2, error.prob = 0.001)

## 5.  Scanone + plotting loop ------------------------------------------------
traits <- names(cross$pheno)[-1]     # drop the ID column
perm      <- readRDS("scanone_1000perm.rds") # load permutation thresholds
lod_thr   <- summary(perm, alpha = 0.05)[1]   # derive the 5% threshold

for (trait in traits) {
  cat("\n=== Processing trait:", trait, "===\n")
  
  # 5a. run scanone (Haley–Knott)
  lod <- scanone(cross, pheno.col = trait, method = "hk")
  peaks     <- summary(lod, perms = perm, alpha = 0.05)
  keep_chr  <- peaks$chr                       # chromosomes with sig QTL
  
  # Skip to next trait if there aren't any significant LOD peaks
  if (length(keep_chr) == 0) {
    print(paste("There were no LOD peaks above the threshold for", trait))
    next
  }
  
  # 5b. open PDF for this trait
  pdf(paste0("QTL_segdist_", trait, ".pdf"), width = 7, height = length(keep_chr)*4.5)
  par(mfrow = c(length(keep_chr), 1),
      mar   = c(5, 4, 2.5, 1),
      cex   = 0.85)
  
  for (chr in keep_chr) {
    ## Pull distortion for this chromosome
    idx <- gt_tab$chr == chr
    x   <- gt_tab$pos[idx]
    y   <- gt_tab$neglog10P[idx]
    
    ## Set a common Y‑limit for overlay
    lod_chr <- lod[lod$chr == chr, "lod"]
    ylim_max <- max(lod_chr, y, na.rm = TRUE)
    
    ## Distortion trace (upper panel style)
    plot(x, y, type = "l", lty = 2, col = "darkgreen",
         ylim = c(0, ylim_max),
         xlab = paste0("Chr ", chr, " position (cM)"),
         ylab = "LOD / –log10(p)",
         main = paste("Trait:", trait, "| Chr", chr))
    abline(h = -log10(0.01), lty = 3, col = "grey60")
    
    ## Overlay LOD curve
    plot(lod, chr = chr, add = TRUE, col = "steelblue", lwd = 2)
    
    legend("topright",
           legend = c("LOD", "Seg. distortion"),
           col    = c("steelblue", "darkgreen"),
           lty    = c(1, 2),
           lwd    = 2,
           bg     = "white", cex = 0.8)
  }
  dev.off()
}






