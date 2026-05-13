## ============================================================
## Generate scanone() permutation RDS (multi-phenotype EM, n.perm)
## Matches qtl-mapping-May2026.R: genoprobs, perm_pheno_cols, method="em", snow cluster.
##
## Command line (run from the qtl_mapping directory):
##
##   cd ~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping
##
##   # Optional: log everything to a file (create logs dir first so tee works)
##   mkdir -p results/logs
##   Rscript scripts/generate_scanone_permutations.R 2>&1 | tee results/logs/generate_scanone_permutations.log
##
##   # Or run without tee (messages only on the terminal)
##   Rscript scripts/generate_scanone_permutations.R
##
## Output path: set `perm_file_out` below (default is timestamped under
## results/processed_data/). Then point qtl-mapping-May2026.R `perm_file_in`
## at that file and keep `perm_pheno_cols` identical in both scripts.
##
## Runtime: multi-trait 1000 permutations often ~10+ minutes (depends on CPU;
## uses snow with n.cluster if the snow package is installed).
## ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(qtl))

if (requireNamespace("snow", quietly = TRUE)) {
  suppressPackageStartupMessages(library(snow))
}

## ----------------------- Configuration ----------------------
qtl_root <- path.expand("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping")
results_dir <- file.path(qtl_root, "results")
processed_dir <- file.path(results_dir, "processed_data")
log_dir <- file.path(results_dir, "logs")
for (d in c(processed_dir, log_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

cross_genfile <- file.path(qtl_root, "data", "mapthis_LG.csv")
cross_phefile <- file.path(qtl_root, "data", "costus_pheno_rqtl_2025Jan24.csv")

genoprob_step <- 0.5
genoprob_error_prob <- 0.001

perm_pheno_cols <- 2:24

perm_n <- 1000L
perm_n_cluster <- 6L

perm_file_out <- file.path(
  processed_dir,
  paste0("scanone_", perm_n, "perm_multitrait_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
)

.ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_msg <- function(...) {
  cat("[", .ts(), "] ", paste0(..., collapse = ""), "\n", sep = "")
}

## --------------------------- Run ----------------------------
log_msg("Output RDS: ", perm_file_out)
if (file.exists(perm_file_out)) {
  stop("Refusing to overwrite existing file. Delete or rename:\n  ", perm_file_out)
}

log_msg("Reading cross …")
cross.data <- read.cross(
  "csvs",
  genfile = cross_genfile,
  phefile = cross_phefile,
  estimate.map = FALSE,
  genotypes = c("AA", "AB", "BB")
)

log_msg("calc.genoprob …")
data_prob <- calc.genoprob(
  cross.data,
  step = genoprob_step,
  error.prob = genoprob_error_prob
)

if (max(perm_pheno_cols) > ncol(data_prob$pheno)) {
  stop(
    "perm_pheno_cols requests column ", max(perm_pheno_cols),
    " but only ", ncol(data_prob$pheno), " phenotype columns exist."
  )
}

data_prob_perm <- data_prob
data_prob_perm$pheno <- data_prob_perm$pheno[, perm_pheno_cols, drop = FALSE]

perm_cols <- colnames(data_prob_perm$pheno)
log_msg(
  "Perm scan phenotypes (n=", length(perm_cols), "): ",
  paste(head(perm_cols, 5), collapse = ", "), ", …"
)

log_msg("scanone() permutations: n.perm=", perm_n, " method=em (this can take a long time) …")
pheno_col_perm <- seq_len(ncol(data_prob_perm$pheno))
perm_args <- list(
  cross = data_prob_perm,
  method = "em",
  n.perm = perm_n,
  pheno.col = pheno_col_perm
)
if ("snow" %in% .packages()) {
  perm_args$n.cluster <- perm_n_cluster
  log_msg("Using snow: n.cluster=", perm_n_cluster)
} else {
  log_msg("snow not attached; single-machine permutations.")
}

t0 <- proc.time()
perm_data <- do.call(scanone, perm_args)
t1 <- proc.time()
log_msg("Elapsed: ", signif((t1 - t0)[3], 4), " s")

log_msg("dim(perm_data): ", paste(dim(perm_data), collapse = " x "))
if (!is.null(colnames(perm_data))) {
  log_msg("colnames(perm_data): ", paste(head(colnames(perm_data), 8), collapse = ", "), ", …")
}

cat("\n--- summary(perm_data) ---\n")
print(summary(perm_data))

saveRDS(perm_data, perm_file_out)
log_msg("Saved: ", perm_file_out)
log_msg("Set perm_file_in in qtl-mapping-May2026.R to this path; keep perm_pheno_cols in sync.")
