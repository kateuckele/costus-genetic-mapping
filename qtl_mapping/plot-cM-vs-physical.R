## 0.  House‑keeping ----------------------------------------------------------
rm(list = ls())
library(qtl)

setwd("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping")

## ---- libraries -------------------------------------------------------------
library(purrr)      # map/imap helpers
library(stringr)    # tidy string handling
library(dplyr)      # piping / data wrangling
library(ggplot2)    # plotting

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

## ---- 2.  Extract marker information ---------------------------------------

marker_df <- imap_dfr(cross$geno, ~{
  map_vec  <- .x$map
  tibble(
    chr      = .y,                               # linkage‑group label
    cM       = as.numeric(map_vec),              # genetic position
    bp       = as.numeric(str_extract(names(map_vec), "(?<=-)[0-9]+")),
    chr_num  = as.numeric(str_extract(names(map_vec), "(?<=Chrom)\\d+"))
  )
})

## ---- 3.  quick sanity check (optional) ------------------------------------
head(marker_df)

## ---- 4.  Plot --------------------------------------------------------------
pdf("~/Dropbox/Costus/costus-genetic-mapping/linkage_map/plot_cM_vs_physical.pdf")
ggplot(marker_df, aes(x = bp/1e6, y = cM,
                      colour = factor(chr_num))) +
  geom_point(size = 0.4, alpha = 0.7) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ chr, scales = "free_x") +
  labs(
    x = "Genomic position (Mb)",
    y = "Genetic position (cM)",
    colour = "Chromosome",
    title = "Genetic (cM) vs. physical (bp) positions, coloured by chromosome"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
dev.off()



