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
morphometric <- pheno[,2:24]
color <- pheno[,25:50]

## change column order
#pheno <- pheno %>% select(ANL, ANW, CLL, COL, COLL, LABL, LABW, RALA, RAST,
#                          STAE, STAL, STATL, STAW, STIW, STYL, TUA, CAL, INFA,
#                          VNG, EFNSC, FNSC, VEFN, VFN)

##Correlations
cor(morphometric, method = "pearson", use = "pairwise.complete.obs")

##Check correlations (as scatterplots), distribution and print correlation coefficient 
ggpairs(morphometric, title="correlogram with ggpairs()") 

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

################################################################################
# Building a correlation network
# --------------------------------

# install.packages(c("igraph", "ggraph", "tidygraph", "ggplot2"))

library(igraph)
library(tidygraph)
library(ggraph)
library(ggplot2)

# correlation matrix
cor_mat <- cor(morphometric, use = "pairwise.complete.obs", method = "pearson")

# convert correlation matrix to edge list
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("from", "to", "r")

# make sure names are characters, not factors
cor_df$from <- as.character(cor_df$from)
cor_df$to   <- as.character(cor_df$to)

# remove self-correlations
cor_df <- cor_df[cor_df$from != cor_df$to, ]

# keep only one copy of each pair
cor_df <- cor_df[cor_df$from < cor_df$to, ]

# create node list so all traits are retained
nodes <- data.frame(name = colnames(cor_mat))

# now make a plotting graph with only stronger edges
cor_df_plot <- subset(cor_df, abs(r) >= 0.5)

cor_df_plot$weight <- abs(cor_df_plot$r)
cor_df_plot$sign <- ifelse(cor_df_plot$r > 0, "positive", "negative")

g_plot <- graph_from_data_frame(
  d = cor_df_plot,
  vertices = nodes,
  directed = FALSE
)

tg_plot <- as_tbl_graph(g_plot)

pdf(
  "~/Dropbox/Costus/costus-genetic-mapping/phenotype/results/figures/morphometric_traits_network_large_legend.pdf",
  height = 10,
  width = 10
)

ggraph(tg_plot, layout = "kk") +
  geom_edge_link(aes(width = weight, alpha = weight, linetype = sign)) +
  geom_node_point(color = "grey70", size = 14) +
  geom_node_text(aes(label = name), repel = TRUE, size = 8) +
  scale_edge_width(
    name = expression("|Correlation|"),
    range = c(0.3, 4.5),
    limits = c(0.5, 1),
    breaks = c(0.5, 0.6, 0.7, 0.8, 0.9)
  ) +
  scale_edge_alpha(range = c(0.25, 0.9), guide = "none") +
  scale_edge_linetype_manual(
    name = "Correlation sign",
    values = c("positive" = "solid", "negative" = "dashed"),
    labels = c("positive" = "Positive", "negative" = "Negative")
  ) +
  guides(
    edge_width = guide_legend(order = 1),
    edge_linetype = guide_legend(order = 2)
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.key.width = unit(1.8, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.spacing.y = unit(0.35, "cm"),
    legend.margin = margin(5, 5, 5, 5),
    plot.margin = margin(10, 10, 10, 10)
  )

dev.off()

# identify 'cliques'
max_cliques(g_plot, min = 3) |> lapply(names)
clique_num(g_plot)
max_cliques(g_plot, min = 3) |> lapply(names)
largest_cliques(g_plot) |> lapply(names)








