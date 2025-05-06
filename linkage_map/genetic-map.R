## clear workspace
rm(list = ls())

## Load libraries
library(qtl)
library(stringr)
library(ggplot2)
library(dplyr)

# ========================================================================
# Step 1. Read in map data
# ========================================================================
mapthis <- read.cross("csv", "~/Dropbox/Costus/costus-genetic-mapping/linkage_map/", "costus_gen_2025April1_unthinned.csv", estimate.map=FALSE)
summary(mapthis)
## The warning message indicates that some chromosomes are greater than 1000 cM in length.
## This is because I have 16,000+ markers assigned to a single linkage group. 
## It is assuming every marker is one centimorgan. We can safely ignore this warning. 

## Total markers:      16769

# ========================================================================
# Step 2. Filter cross data
# ========================================================================

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 1. REMOVE IDENTICAL MARKERS
## Find sets of markers that have identical genotype data and that are adjacent to each other
dup <- findDupMarkers(mapthis, exact.only=FALSE, adjacent.only=TRUE)
## Drop adjacent markers that are duplicated
mapthis <- drop.markers(mapthis, unlist(dup))
summary(mapthis)
## Total markers:      4114

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 2. REMOVE INDIVIDUALS AND MARKERS WITH HIGH LEVELS OF MISSING DATA
## Plot missingness across individuals and markers
par(mfrow=c(1,2), las=1)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",
       main="No. genotypes by marker")

## Keep individuals with over 1000 typed markers
mapthis <- subset(mapthis, ind=(ntyped(mapthis)>1000))

## Drop sites typed in less than 450 individuals 
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 450])
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)
## Total markers:      4051

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 3. REMOVE UNUSUALLY SIMILAR INDIVIDUALS 
## Count proportion of matching genotypes between all pairs of individuals
cg <- comparegeno(mapthis, what = "proportion")
## Plot a histogram of the proportions of matching genotypes for all pairs of individuals
par(mfrow=c(1,1), las=1)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

## Identify which pairs of individuals have higher than 90% matching genotypes
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

## Extract the genotype matrix from the cross object 'mapthis' 
## Create a contingency table comparing the genotype calls between individuals 383 and 423 across all markers.
g <- pull.geno(mapthis)
table(g[383,], g[423,])

## Mark positions where both individuals have non-missing genotypes and their genotypes differ
## Set the genotype data for the first individual to missing (NA) at those positions
tozero <- !is.na(g[wh[1],]) & !is.na(g[wh[2],]) & g[wh[1],] != g[wh[2],]
mapthis$geno[[1]]$data[wh[1],tozero] <- NA

## Remove the second individual, keeping the first
mapthis <- subset(mapthis, ind=-wh[2])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 4. INVESTIGATE PATTERNS OF SEGREGATION DISTORTION 

## Create new cross object with linkage groups assigned according to C. lasius genome
newLG <- substr(names(mapthis$geno$`1`$map), 6, 6)
# Extract the original marker data from the single linkage group
orig_data <- mapthis$geno[[1]]$data
orig_map <- mapthis$geno[[1]]$map
new_map <- as.numeric(substr(names(mapthis$geno$`1`$map), 8, nchar(names(mapthis$geno$`1`$map))))
names(new_map) <- names(orig_map)
# Get the unique linkage groups
uniqueLG <- sort(unique(newLG))
# Create a new list for the genotype data by splitting markers according to newLG
newgeno <- lapply(uniqueLG, function(lg) {
  # Identify the marker indices corresponding to this new linkage group
  marker_indices <- which(newLG == lg)
  # Subset the original data and map for these markers.
  list(
    data = orig_data[, marker_indices, drop = FALSE],
    map  = new_map[marker_indices]
  )
})
names(newgeno) <- uniqueLG
# Create a new cross object with the re-organized genotype data.
new_cross <- mapthis
new_cross$geno <- newgeno

## Plot patterns of segregation distortion across the genome
## Create a table of genotype distributions (including P-values from chi-square tests for Mendelian segregation)
gt <- geno.table(new_cross, scanone.output=TRUE)
## Plot p-values (top) and genotype frequencies (bottom)
pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/segregation_distortion_Clasius_groups.pdf")
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")), gap=10, bandcol="gray70") # plot.scanone
plot(gt, chr = 2, lod=3:5, ylab="Genotype frequency", gap=10, bandcol="gray70") # plot.scanone
abline(h=c(0.25, 0.5), lty=2, col="gray")
dev.off()

## Mark translocation region on Chromosome 2
plot(gt, type="p", xlim=c(130e6,163e6), chr = "2",lod=3:5, ylab="Genotype frequency", gap=10, bandcol="gray70", col = c("black", "blue", "red")) # plot.scanone
abline(h=c(0.25, 0.5), lty=2, col="gray")
abline(v=c(145249000, 156891000), lty=2, col="red")
title(main = "Chromosome 2")
legend("bottomleft",                # or "bottomleft", "topleft", etc.
       legend = colnames(gt)[5:7],# text keys; use the same columns you plotted
       col    = c("black", "blue", "red"),          # matching colours
       lty    = 1,                # all lines use linetype 1 in plot.scanone
       lwd    = 2,                # same line width you get from plot()
       bty    = "n")              # drop the legend box (optional)

## Mark translocation region on Chromosome 5
plot(gt, type = "p", xlim=c(78e6, 108e6), chr = "5",lod=3:5, ylab="Genotype frequency", gap=10, bandcol="gray70", col = c("black", "blue", "red")) # plot.scanone
abline(h=c(0.25, 0.5), lty=2, col="gray")
abline(v=c(85210000, 91687000), lty=2, col="red")
title(main = "Chromosome 5")
legend("bottomright",                # or "bottomleft", "topleft", etc.
       legend = colnames(gt)[5:7],# text keys; use the same columns you plotted
       col    = c("black", "blue", "red"),          # matching colours
       lty    = 1,                # all lines use linetype 1 in plot.scanone
       lwd    = 2,                # same line width you get from plot()
       bty    = "n")              # drop the legend box (optional)

## Mark translocation region on Chromosome 9
plot(gt, type = "p", xlim=c(74e6, 82e6), chr = "9",lod=3:5, ylab="Genotype frequency", gap=10, bandcol="gray70", col = c("black", "blue", "red")) # plot.scanone
abline(h=c(0.25, 0.5), lty=2, col="gray")
abline(v=c(75457000, 77619000), lty=2, col="red")
title(main = "Chromosome 9")
legend("bottomright",                # or "bottomleft", "topleft", etc.
       legend = colnames(gt)[5:7],# text keys; use the same columns you plotted
       col    = c("black", "blue", "red"),          # matching colours
       lty    = 1,                # all lines use linetype 1 in plot.scanone
       lwd    = 2,                # same line width you get from plot()
       bty    = "n")              # drop the legend box (optional)

## Is ancestry biased across the genome? 
sum(gt$AA)
# [1] 1042.686 (should be 1013)
sum(gt$BB)
# [1] 918.2241 (should be 1013)
sum(gt$AB)
# [1] 2090.09 (should be 2026)

total <- sum(gt$AA) + sum(gt$AB) + sum(gt$BB)
total/4
# [1] 1012.75

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 5. REMOVE MARKERS WITH SIGNALS OF DISTORTED SEGREGATION 
## We expect the genotypes to appear with the frequencies 1:2:1. 
## Moderate departures from these frequencies are not unusual and may indicate the presence of partially lethal alleles. 
## Gross departures from these frequencies often indicate problematic markers that should be omitted, at least initially.

## Create a table of genotype distributions (including P-values from chi-square tests for Mendelian segregation)
gt <- geno.table(mapthis)
## Subset the markers with significant p-values (after adjusting for multiple tests)
gt_sig <- gt[gt$P.value < 0.05/totmar(mapthis),]
## Plot a histogram of significant p-values
hist(gt_sig$P.value, main="Histogram of significant p-values indicating segregation distortion")
## Omit the worst of these markers
todrop <- rownames(gt[gt$P.value < 4e-6,]) ## 418 markers
mapthis <- drop.markers(mapthis, todrop)
summary(mapthis)
## Total markers:      3633 

### Plot genotype frequencies
## pull out raw genotype data
g <- pull.geno(mapthis)
## for each individual, calculate the proportion of each genotype (AA, AB, BB)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
## for each genotype, plot average individual genotype frequency
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
## One individual has 100% AB genotype, but this is OK because this is the F1

gt <- geno.table(mapthis, scanone.output = TRUE)
sum(gt$AA)
# [1] 927.1093
sum(gt$BB)
# [1] 834.7141
sum(gt$AB)
# [1] 1871.177

total <- sum(gt$AA) + sum(gt$AB) + sum(gt$BB)

total/2
# [1] 908.25

# ========================================================================
# Step 3. Form linkage groups
# ========================================================================

## Estimate pairwise recombination fractions between all pairs of genetic markers
mapthis <- est.rf(mapthis)
## Extract recombination fractions and LOD scores
rf <- pull.rf(mapthis)
lod <- pull.rf(mapthis, what="lod")

## Plot LOD scores by recombination fractions
pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/LODbyRF.pdf")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
dev.off()

## Form linkage groups
## Two markers will be placed in the same linkage groups if they have estimated RF ≤ max.rf and LOD score ≥ min.lod
## Try a range of min LOD values: 
for (i in seq(10, 200, 10)) {
  lg <- formLinkageGroups(mapthis, max.rf=0.4, min.lod=i)
  print(paste("min LOD is", i))
  print(table(lg[,2]))
}

## min LOD of 50 formed 9 large linkage groups and one small linkage group:
## 1   2   3   4   5   6   7   8   9    10 
## 522 470 417 411 410 401 396 318 262  26 

## Try a range of max RF values: 
for (i in seq(0.1, 0.8, 0.1)) {
  lg <- formLinkageGroups(mapthis, max.rf=i, min.lod=50)
  print(paste("min LOD is", i))
  print(table(lg[,2]))
}

## Results essentially unchanged

## Assign min.LOD to 50, since we have 9 chromosomes
## Reorganize markers into these new linkage groups
mapthis_LG <- formLinkageGroups(mapthis, max.rf=0.4, min.lod=50, reorgMarkers=TRUE)

## Investigate the association between C. lasius chromosome numbers and linkage groups
# Create a two-column matrix: first column is the C. lasius chromosome (extracted from row names),
# second column is the linkage group number.
LG <- formLinkageGroups(mapthis, max.rf=0.4, min.lod=50)
print(table(LG[,2]))
lg_before_after <- cbind(as.numeric(substr(rownames(LG), 6, 6)), LG$LG)
colnames(lg_before_after) <- c("lasius_chrom", "LG")

# For each linkage group, create a table of the associated C. lasius chromosomes.
assoc_tables <- lapply(sort(unique(lg_before_after[, "LG"])), function(x) {
  table(lg_before_after[lg_before_after[, "LG"] == x, "lasius_chrom"])
})
names(assoc_tables) <- paste0("LG", sort(unique(lg_before_after[, "LG"])))
assoc_tables  # view the resulting tables

## Rename the linkage groups based on their associated C. lasius chromosome numbers.
new_names <- c("2", "4.2.9", "5", "7.5", "3", "6", "8", "1", "9", "10")
names(mapthis_LG$geno)[1:length(new_names)] <- new_names

## Check that linkage groups were renamed appropriately
chrnames(mapthis_LG)

## Plot a grid showing the recombination fractions (lower triangle) and LOD scores (upper triangle) for all pairs of markers
pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/plotRF.pdf")
plotRF(mapthis_LG, alternate.chrid=TRUE)
dev.off()

## Move markers on linkage group 10 to linkage group 1
c10markers <- markernames(mapthis_LG, chr=10)
for(marker in c10markers) {
  mapthis_LG <- movemarker(mapthis_LG, marker, 1)
}

## Write the cross to a csv file
write.cross(mapthis_LG, format="csv", filestem="~/Dropbox/Costus/costus-genetic-mapping/linkage_map/mapthis_LG")

# ========================================================================
# Step 4. Order markers on chromosomes
# ========================================================================

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 1. ORDER MARKERS ON CHROMOSOME 4.2.9
## Order markers and pull the map for composite chromosome 4.2.9
#mapthis_LG <- orderMarkers(mapthis_LG, chr = "4.2.9", window = 3)
df <- data.frame(
  snp = names(pull.map(mapthis_LG, chr = "4.2.9")[[1]]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    chromosome = as.numeric(str_extract(snp, "(?<=Chrom)\\d+")),
    position   = as.numeric(str_extract(snp, "(?<=-)\\d+")),
    order      = row_number()
  )
# Plot positions sequentially, colored by chromosome
ggplot(df, aes(x = order, y = position, color = factor(chromosome))) +
  geom_point(size = 3) +
  labs(x = "SNP Order", y = "Position", color = "Chromosome") +
  theme_minimal()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 2. ORDER MARKERS ON CHROMOSOME 7.5
# Order markers and pull the map for composite chromosome 7.5
#mapthis_LG <- orderMarkers(mapthis_LG, chr="7.5", window=3)
df <- data.frame(
  snp = names(pull.map(mapthis_LG, chr = "7.5")[[1]]),
  stringsAsFactors = FALSE
) %>%
  mutate(
    chromosome = as.numeric(str_extract(snp, "(?<=Chrom)\\d+")),
    position   = as.numeric(str_extract(snp, "(?<=-)\\d+")),
    order      = row_number()
  )

# Plot positions sequentially, colored by chromosome
ggplot(df, aes(x = order, y = position, color = factor(chromosome))) +
  geom_point(size = 3) +
  labs(x = "SNP Order", y = "Position", color = "Chromosome") +
  theme_minimal()

## Plot a grid showing the recombination fractions (lower triangle) and LOD scores (upper triangle) for all pairs of markers
pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/plotRF.2.pdf")
plotRF(mapthis_LG, alternate.chrid=TRUE)
dev.off()

# ========================================================================
# Step 5. Investigate the number of crossovers per individual
# ========================================================================
## F2 individuals result from two independent meiotic events (one from each parent).
## Meiosis is expected to generate 1-2 crossover per chromosome. 
## --> An F2 individual is expected to have 2-4 crossovers per chromosome. 
## 2-4 crossovers x 9 chromosomes = ~18-36 crossovers

## Count the observed number of crossovers in each individual
plot(countXO(mapthis_LG), ylab="Number of crossovers")
## Calculate the median number of crossovers
quantile(countXO(mapthis_LG))
#   0%  25%  50%  75%  100% 
#   0   14   16   19   30 

# ========================================================================
# Step 6. Estimate the genetic map
# ========================================================================

#newmap <- est.map(mapthis_LG, error.prob=0.005)
#mapthis <- replace.map(mapthis, newmap)
summaryMap(newmap)

newmap_copy <- newmap
mapthis_LG_copy <- mapthis_LG

mapthis_LG <- replace.map(mapthis_LG, newmap)

pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/plotMap.pdf")
plotMap(mapthis_LG, show.marker.names=FALSE)
dev.off()

pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/plotMap.4.2.9.pdf", width = 60)
plotMap(mapthis_LG, show.marker.names=TRUE, chr = "4.2.9", horizontal = TRUE)
dev.off()

pdf("/Users/kathrynuckele/Dropbox/Costus/costus-genetic-mapping/linkage_map/plotMap.7.5.pdf", width = 60)
plotMap(mapthis_LG, show.marker.names=TRUE, chr = "7.5", horizontal = TRUE)
dev.off()

## Write the cross to a csv file
write.cross(mapthis_LG, format="csv", filestem="~/Dropbox/Costus/costus-genetic-mapping/linkage_map/mapthis_LG")

# ========================================================================
# Step 7. Investigate patterns of segregation distortion
# ========================================================================

#If we apply a Bonferroni correction for the 88 tests (88 is the total number of 
# markers we have retained in the data), we would look for P ≥ 0.05/88 which 
## corresponds to −log10 P ≥ 3.25

mapthis_LG <- read.cross("csv", dir = "~/Dropbox/Costus/costus-genetic-mapping/linkage_map/", 
                         file="mapthis_LG.csv",
                         estimate.map=FALSE, genotypes=c("AA","AB","BB"))

gt <- geno.table(mapthis_LG, scanone.output=TRUE)
par(mfrow=c(2,1))
plot(gt, ylab=expression(paste(-log[10], " P-value")), gap=10, bandcol="gray70") # plot.scanone
plot(gt, lod=3:5, ylab="Genotype frequency", gap=10, bandcol="gray70") # plot.scanone
abline(h=c(0.25, 0.5), lty=2, col="gray")
legend("bottomright",                # or "bottomleft", "topleft", etc.
       legend = colnames(gt)[5:7],# text keys; use the same columns you plotted
       col    = c("black", "blue", "red"),          # matching colours
       lty    = 1,                # all lines use linetype 1 in plot.scanone
       lwd    = 2,                # same line width you get from plot()
       bty    = "n")              # drop the legend box (optional)







