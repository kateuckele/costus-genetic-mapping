########################################################################################
# GC CONTENT ANALYSIS FOR TRANSLOCATIONS IN COSTUS
# Author: K. Uckele
# Date: 2/27/2025
# Description:
# This script analyzes GC content across specific genomic regions of interest (ROIs)
# and compares GC content between different regions within each chromosome.
#
# Main Steps:
# 1. Load GC content data from "gc_content.txt" (generated from bedtools nuc).
# 2. Read a BED file ("roi.bed") containing regions of interest.
# 3. Assign each GC content window to a specific ROI or mark it as "Background."
# 4. Compute summary statistics (mean, median, standard deviation) for each ROI.
# 5. Generate boxplots comparing GC content distributions across regions.
# 6. Perform pairwise statistical comparisons (Wilcoxon test and t-test) for ROIs.
# 7. Apply Holm-Bonferroni correction to adjust for multiple testing.
# 8. Determine which ROI has significantly higher GC content in each pairwise test.
# 9. Save summary statistics and statistical test results to CSV files.
#
# Required R Packages:
# - dplyr (data manipulation)
# - ggplot2 (visualization)
# - readr (reading files)
# - tidyr (data tidying)
#
# Inputs:
# - "gc_content.txt" (output from bedtools nuc, contains GC content per genomic window)
# - "roi.bed" (BED file specifying regions of interest)
#
# Outputs:
# - "gc_summary_results.csv" (summary statistics for each ROI)
# - "pairwise_gc_tests_holm_intra_chromosomal.csv" (statistical test results)
# - Boxplots comparing GC content distributions across regions
#
# Notes:
# - The analysis is performed separately for chromosomes 2, 5, 7, and 9.
# - Statistical significance is tested using both Wilcoxon and t-tests.
# - Holm-Bonferroni correction is applied to reduce false positives in multiple testing.
# - The script automatically determines which ROI has significantly higher GC content.
########################################################################################

# set the working directory to where the GC content data and BED files are stored
setwd("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/translocations/gc_content_translocations/")

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Read GC content data from bedtools output
gc_data <- read_tsv("gc_content.txt", col_names = TRUE)

# Assign column names based on bedtools nuc output
colnames(gc_data) <- c("chr", "start", "end", "pct_at", "pct_gc", 
                       "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len")

# Read Region of Interest (ROI) BED file
roi <- read_tsv("roi.bed", col_names = FALSE)

# Assign column names to the ROI data
colnames(roi) <- c("chr", "roi_start", "roi_end", "roi_name")

# Function to determine which ROI a given GC content window belongs to
assign_roi <- function(chr, start, end, roi_df) {
  matches <- roi_df %>%
    filter(chr == !!chr & roi_start <= end & roi_end >= start)
  
  if (nrow(matches) > 0) {
    return(matches$roi_name[1])  # Assign the first matching ROI name
  } else {
    return("Background")  # Assign as background if no match
  }
}

# Apply the function to classify each GC content window into an ROI or background
gc_data$region <- mapply(assign_roi, gc_data$chr, gc_data$start, gc_data$end, MoreArgs = list(roi_df = roi))

# Compute summary statistics (mean, standard deviation, median, and count) for each ROI
summary_gc <- gc_data %>%
  group_by(region) %>%
  summarise(mean_gc = round(mean(pct_gc),4), sd_gc = round(sd(pct_gc),4), median_gc = round(median(pct_gc),4), n = n())

# Print summary statistics
print(summary_gc)

# Save summary to a CSV file
write.csv(summary_gc, "gc_summary_results.csv", row.names = FALSE)

# Create boxplots to compare GC content across ROIs for different chromosomes

ggplot(gc_data[gc_data$chr == "Chrom2",], aes(x = region, y = pct_gc, fill = region)) +
  geom_boxplot() +
  labs(title = "GC Content in ROIs vs. Background",
       x = "Region",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gc_data[gc_data$chr == "Chrom5",], aes(x = region, y = pct_gc, fill = region)) +
  geom_boxplot() +
  labs(title = "GC Content in ROIs vs. Background",
       x = "Region",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gc_data[gc_data$chr == "Chrom7",], aes(x = region, y = pct_gc, fill = region)) +
  geom_boxplot() +
  labs(title = "GC Content in ROIs vs. Background",
       x = "Region",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(gc_data[gc_data$chr == "Chrom9",], aes(x = region, y = pct_gc, fill = region)) +
  geom_boxplot() +
  labs(title = "GC Content in ROIs vs. Background",
       x = "Region",
       y = "GC Content (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Define chromosomes of interest
chromosomes <- c("Chrom2", "Chrom5", "Chrom7", "Chrom9")

# Create empty lists to store results for statistical tests
wilcox_pvalues <- c()
t_test_pvalues <- c()
comparisons <- c()
gc_comparison_wilcoxon <- c()  
gc_comparison_ttest <- c()  

# Perform pairwise comparisons within each chromosome
for (chr in chromosomes) {
  # Filter GC content data for the current chromosome
  chr_gc_data <- gc_data %>% filter(chr == !!chr)
  
  # Get unique ROIs on this chromosome
  roi_names <- unique(chr_gc_data$region)
  roi_names <- roi_names[roi_names != "Background"]  # Exclude Background if not needed
  
  # Pairwise comparisons within the same chromosome
  for (i in 1:(length(roi_names) - 1)) {
    for (j in (i + 1):length(roi_names)) {
      roi1 <- roi_names[i]
      roi2 <- roi_names[j]
      
      # Extract GC content for each ROI
      roi1_gc <- chr_gc_data$pct_gc[chr_gc_data$region == roi1]
      roi2_gc <- chr_gc_data$pct_gc[chr_gc_data$region == roi2]
      
      # Wilcoxon test (non-parametric)
      wilcox_test <- wilcox.test(roi1_gc, roi2_gc, conf.int=TRUE)
      wilcox_pvalues <- c(wilcox_pvalues, wilcox_test$p.value)
      
      # T-test (parametric)
      t_test <- t.test(roi1_gc, roi2_gc)
      t_test_pvalues <- c(t_test_pvalues, t_test$p.value)
      
      # Store the comparison name
      comparisons <- c(comparisons, paste(chr, roi1, "vs", roi2))
      
      # Determine which ROI has higher GC content based on Wilcoxon test
      gc_direction_wilcoxon <- "Not significant"
      if (wilcox_test$p.value < 0.01) {
        if (wilcox_test$conf.int[1] > 0) {
          gc_direction_wilcoxon <- paste("Higher in", roi1, "than", roi2)
        } else {
          gc_direction_wilcoxon <- paste("Higher in", roi2, "than", roi1)
        }
      }
      
      # Determine which ROI has higher GC content based on t-test
      gc_direction_ttest <- "Not significant"
      if (t_test$p.value < 0.01) {
        if (mean(roi1_gc, na.rm = TRUE) > mean(roi2_gc, na.rm = TRUE)) {
          gc_direction_ttest <- paste("Higher in", roi1, "than", roi2)
        } else {
          gc_direction_ttest <- paste("Higher in", roi2, "than", roi1)
        }
      }
      
      # Store results
      gc_comparison_wilcoxon <- c(gc_comparison_wilcoxon, gc_direction_wilcoxon)
      gc_comparison_ttest <- c(gc_comparison_ttest, gc_direction_ttest)
    }
  }
}

# Apply Holm-Bonferroni correction for multiple comparisons
wilcox_adj_p <- p.adjust(wilcox_pvalues, method = "holm")
t_test_adj_p <- p.adjust(t_test_pvalues, method = "holm")

# Create a dataframe to store test results
pairwise_results <- data.frame(
  Comparison = comparisons,
  Wilcoxon_p = wilcox_pvalues,
  Wilcoxon_adj_p = wilcox_adj_p,
  Ttest_p = t_test_pvalues,
  Ttest_adj_p = t_test_adj_p,
  gc_comparison_wilcoxon = gc_comparison_wilcoxon, 
  gc_comparison_ttest = gc_comparison_ttest 
)

# Print results
print(pairwise_results)

# Save results to CSV file
write.csv(pairwise_results, "pairwise_gc_tests_holm_intra_chromosomal.csv", row.names = FALSE)

