# -------------------------------------------------------------------
# Generate R/QTL Phenotype Input File
# Author: K. Uckele
# Date: July 24, 2024
# -------------------------------------------------------------------

# ----------------------------
# 1. Setup Environment
# ----------------------------

# Clear the workspace to ensure no residual variables affect the script
rm(list = ls())

# Load necessary libraries
# Ensure that the required packages are installed. If not, install them.
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}

# Define a custom operator `%notin%` for convenience
`%notin%` <- Negate(`%in%`)

# ----------------------------
# 2. Load Data
# ----------------------------

# Define file paths (adjust paths as necessary)
pheno_path <- "~/Dropbox/Costus/genetic_mapping/phenotypes/phenotypes_nooutliers.csv"
color_path <- "~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/processed_data/spectral_shape_descriptors.csv"
plant_lab_map_path <- "~/Dropbox/Costus/genetic_mapping/R:QTL/plate_prep_for_Davis.csv"
gen_path <- "~/Dropbox/Costus/genetic_mapping/R:QTL/costus_gen_2024July23.csv"

# Read CSV files into data frames
# It's good practice to set stringsAsFactors = FALSE to prevent automatic factor conversion
pheno <- read.csv(pheno_path, header = TRUE, stringsAsFactors = FALSE)
color <- read.csv(color_path, header = TRUE, stringsAsFactors = FALSE)
plant_lab_map <- read.csv(plant_lab_map_path, header = TRUE, stringsAsFactors = FALSE)
gen <- read.csv(gen_path, stringsAsFactors = FALSE)  # Assuming no header is needed; otherwise, add header = TRUE

# ----------------------------
# 3. Data Preparation
# ----------------------------

# Rename the first column of `pheno` to "id" for consistency
# Ensure that the first column exists to prevent errors
if (ncol(pheno) >= 1) {
  names(pheno)[1] <- "id"
} else {
  stop("The `pheno` data frame does not have any columns.")
}

# Combine the `pheno` and `color` data frames using a full join on "id"
# Ensure that the dplyr package is loaded for `full_join`
joined_pheno_color <- pheno %>%
  full_join(color, by = "id")

# -------------------------------------
# 4. Filter Hybrids in Plant-Lab-ID Map
# -------------------------------------

# Extract hybrids from the Plant-Lab-ID map where PlantID contains "x"
hybrids <- plant_lab_map %>%
  filter(grepl("x", PlantID)) %>%
  select(PlantID, LabID)

# ----------------------------------
# 5. Remove F1 Hybrids (keeping F2s)
# ----------------------------------

# Define the F1 hybrids to remove
f1_hybrids_to_remove <- c("Costus lasius x bracteatus #62", "Costus bracteatus x lasius #39")

# Remove the specified F1 hybrids from the `hybrids` data frame
hybrids <- hybrids[!hybrids$PlantID %in% f1_hybrids_to_remove, ]

# ----------------------------
# 6. Check Sequenced F2 Samples
# ----------------------------

# Identify genotyped F2s not present in the Plant-Lab-ID map
missing_f2_samples <- gen$id[gen$id %notin% hybrids$LabID]
print(missing_f2_samples)
# [1] "" "" --> indicates no genotyped samples are missing

# ----------------------------
# 7. Handle Duplicated LabIDs
# ----------------------------

# Count the number of duplicated LabIDs in the Plant-Lab-ID map
num_duplicated_labids <- sum(duplicated(hybrids$LabID))
print(num_duplicated_labids)  
# [1] 1 --> one labID is duplicated 

# Identify which LabID is duplicated
duplicated_labid <- hybrids$LabID[duplicated(hybrids$LabID)]
print(duplicated_labid)  
# "21_054" is duplicated

# Remove the second occurrence of the duplicated LabID
nrow(hybrids)
# [1] 498
hybrids <- hybrids[!duplicated(hybrids$LabID), ]
# Validate the second occurrence was removed
nrow(hybrids)
# [1] 497

# ----------------------------
# 8. Handle Duplicated PlantIDs
# ----------------------------

# Count the number of duplicated PlantIDs in the Plant-Lab-ID map
num_duplicated_plantids <- sum(duplicated(hybrids$PlantID))
print(num_duplicated_plantids)  
# [1] 8 --> 8 plantIDs are duplicated

# Identify which PlantIDs are duplicated
duplicated_plantids <- hybrids$PlantID[duplicated(hybrids$PlantID)]
print(duplicated_plantids)
# [1] "39 x 39 -114" "39 x 39 -129" "39 x 39 -115" "39 x 39 -49" 
# [5] "39 x 39 -50"  "39 x 39 -51"  "39 x 39 -123" "39 x 39 -117"

# Remove duplicated PlantIDs, keeping the first occurrence
hybrids <- hybrids[!duplicated(hybrids$PlantID), ]

# Verify the number of rows after removal
print(nrow(hybrids))  
# [1] 489

# ----------------------------
# 9. Reformat PlantIDs
# ----------------------------

# The PlantID format is assumed to be "F1_parent x F2_individual"
# Use this format to extract the F1 parent and F2 individual numbers

# Split PlantID by " x " to extract F1 parent
plantID_first <- sapply(strsplit(hybrids$PlantID, " x "), `[`, 1)

# Split PlantID by "-" to extract F2 individual
plantID_second <- sapply(strsplit(hybrids$PlantID, "-"), `[`, 2)

# Check if both parts were successfully extracted
if (length(plantID_first) != nrow(hybrids) | length(plantID_second) != nrow(hybrids)) {
  stop("Error in splitting PlantID into F1 parent and F2 individual.")
}

# Reassign the PlantID column in the format "F1parent_F2individual"
hybrids$PlantID <- paste0(plantID_first, "_", plantID_second)

# --------------------------------------------------------
# 10. Cross-Check Phenotyped Plants with Plant-Lab-ID Map
# --------------------------------------------------------

# Check how many phenotyped plants aren't in the Plant-Lab-ID Map
num_pheno_not_in_plate <- sum(joined_pheno_color$id %notin% hybrids$PlantID)
print(num_pheno_not_in_plate)  
# [1] 12 --> 12 plants were phenotyped but aren't in map

# Identify which phenotyped plants aren't in the Plant-Lab-ID Map
phenos_not_in_plate <- joined_pheno_color$id[joined_pheno_color$id %notin% hybrids$PlantID]
print(phenos_not_in_plate)
# [1] "39_135" "39_60"  "39_8"   "62_0"   "39_56"  "62_354" "62_54"  "62_89"  "125"    "126"    "39"    
# [12] "BRAC"

# Remove phenotyped plants that weren't sequenced (i.e., not in the Plant-Lab-ID Map)
joined_pheno_color <- joined_pheno_color[joined_pheno_color$id %in% hybrids$PlantID, ]

# --------------------------------------------------------
# 11. Cross-Check Plant-Lab-ID Map with Phenotyped Plants
# --------------------------------------------------------

# Count how many plants in the Plant-Lab-ID Map aren't phenotyped
num_plate_not_pheno <- sum(hybrids$PlantID %notin% joined_pheno_color$id)
print(num_plate_not_pheno)  
# [1] 231 --> 231 plants were sequenced but not phenotyped

# Count how many plants in the Plant-Lab-ID Map are sequenced and phenotyped
num_plate_pheno <- sum(hybrids$PlantID %in% joined_pheno_color$id)
print(num_plate_pheno)  
# [1] 258 --> 258 plants were sequenced and phenotyped

# --------------------------------------------------------
# 12. Add F1 Parent Information (covariate in R/QTL models)
# --------------------------------------------------------

# Extract the F1 parent identifier from the first two characters of the 'id'
# Ensure that 'id' has at least two characters to prevent errors
joined_pheno_color$F1parent <- substr(joined_pheno_color$id, 1, 2)

# --------------------------------------------------------
# 13. Map PlantIDs to LabIDs in Phenotype File
# --------------------------------------------------------

# Substitute the LabID for each PlantID in the phenotype file
# This ensures that phenotype ids align with genotype ids

# Create a lookup table from hybrids
hybrid_lookup <- hybrids %>% select(PlantID, LabID)

# Join the phenotype data with the hybrid lookup to replace 'id' with 'LabID'
joined_pheno_color <- joined_pheno_color %>%
  left_join(hybrid_lookup, by = c("id" = "PlantID")) %>%
  # Replace 'id' with 'LabID' from the lookup
  mutate(id = LabID) %>%
  # Remove the temporary 'LabID' column
  select(-LabID)

# ----------------------------
# 14. Verify ID Consistency Between Phenotype and Genotype Data
# ----------------------------

# Check if all IDs in the phenotype file are present in the genotype file
all_ids_match <- sum(joined_pheno_color$id %in% gen$id) == nrow(joined_pheno_color)
print(all_ids_match)  
# [1] TRUE --> all IDs in the phenotype file are present in the genotype file

# ----------------------------
# 15. Final Formatting for R/QTL
# ----------------------------

# Replace any NA values with "-" as required by R/QTL
joined_pheno_color[is.na(joined_pheno_color)] <- "-"

# ----------------------------
# 16. Export the Prepared Data
# ----------------------------

# Define the output file path
output_path <- "~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/processed_data/costus_pheno_rqtl_2025Jan24.csv"

# Write the prepared phenotype data to a CSV file
# Use quote = FALSE to prevent quoting of character fields
# Use row.names = FALSE to omit row numbers in the output
write.csv(joined_pheno_color, output_path, quote = FALSE, row.names = FALSE)
