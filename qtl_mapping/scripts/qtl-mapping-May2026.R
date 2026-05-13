## ============================================================
## QTL mapping pipeline (R/qtl): scanone, fitqtl, plots, outputs
## Outputs: plots → results/plots; runtime tee logs → results/logs; per-trait transcripts →
##   results/trait_qtl_reports/*.txt; peak interval tables →
##   results/peak_intervals/<YYYYMMDD_alpha*/*.tsv (one subfolder per date×alpha; same day+alpha overwrites)
## Permutations: run scripts/generate_scanone_permutations.R first; set perm_file_in below.
##
## How to run (terminal): from the qtl_mapping directory (the folder that contains scripts/ and results/),
##   mkdir -p results/logs
##   Rscript scripts/qtl-mapping-May2026.R
## Optional — save a full run transcript under results/logs/:
##   Rscript scripts/qtl-mapping-May2026.R 2>&1 | tee results/logs/qtl-mapping-May2026_$(date +%Y%m%d_%H%M%S).log
## Paths inside the script use qtl_root below; adjust qtl_root if your clone is not under ~/Dropbox/...
## ============================================================

## -------------------------- Setup ---------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(qtl)
})

qtl_root <- path.expand("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping")
results_dir <- file.path(qtl_root, "results")
processed_dir <- file.path(results_dir, "processed_data")
plot_dir <- file.path(results_dir, "plots")
log_dir <- file.path(results_dir, "logs")
## Per-trait sink() transcripts (scanone, fitqtl, etc.): plain text, not runtime logs.
trait_report_dir <- file.path(results_dir, "trait_qtl_reports")
peak_intervals_dir <- file.path(results_dir, "peak_intervals")
for (d in c(results_dir, plot_dir, log_dir, trait_report_dir, peak_intervals_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
setwd(results_dir)

## ----------------------- Configuration -----------------------
## NULL = all traits in perm_pheno_cols except covar_name; or c("INFA", ...).
trait_cols <- NULL

alpha <- 0.05

genoprob_step <- 0.5
genoprob_error_prob <- 0.001

cross_genfile <- file.path(qtl_root, "data", "mapthis_LG.csv")
cross_phefile <- file.path(qtl_root, "data", "costus_pheno_rqtl_2025Jan24.csv")

covar_name <- "F1parent"

## Must match generate_scanone_permutations.R (same perm_pheno_cols / RDS column order)
perm_file_in <- file.path(
  processed_dir,
  "scanone_1000perm_multitrait_20260513_141534.rds"
)
perm_pheno_cols <- 2:24

## Interval parameters (reported after scanone)
lodint_drops <- c(`lod_1.5` = 1.5, `lod_2` = 2.0)
bayes_prob <- 0.95

## Plot outputs (results/plots); qqman helper lives under project scripts/
qqman_source <- file.path(qtl_root, "scripts", "adapted_qqman.R")

## -------------------------- Logging --------------------------
.ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

log_msg <- function(...) {
  cat("[", .ts(), "] ", paste0(..., collapse = ""), "\n", sep = "")
}

log_section <- function(title) {
  cat("\n", "---- ", title, " ----", "\n", sep = "")
}

log_print <- function(x, label = NULL) {
  if (!is.null(label)) log_msg(label)
  cat(paste(capture.output(print(x)), collapse = "\n"), "\n")
}

start_trait_log <- function(trait_col) {
  report_file <- file.path(
    trait_report_dir,
    paste0("qtl_trait_report_", trait_col, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  )
  sink(report_file, split = TRUE)
  log_msg("Trait report file: ", report_file)
  list(
    file = report_file,
    close = function() sink()
  )
}

## ------------------------ Helpers ---------------------------
make_add_formula <- function(n_qtl, covar_name = "F1parent") {
  q_terms <- paste0("Q", seq_len(n_qtl))
  rhs <- c(q_terms, covar_name)
  as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

make_formula_with_int <- function(n_qtl, interactions = NULL, covar_name = "F1parent") {
  q_terms <- paste0("Q", seq_len(n_qtl))
  rhs <- c(q_terms, interactions, covar_name)
  as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

safe_has_pheno <- function(cross, pheno_name) {
  !is.null(cross$pheno) && (pheno_name %in% colnames(cross$pheno))
}

## summary(scanone, perms=...) with addcovar can error on multi-column scanoneperm (rqtl `[` bug).
perm_column <- function(perm_data, j) {
  if (!inherits(perm_data, "scanoneperm")) {
    return(perm_data)
  }
  nr <- nrow(perm_data)
  nc <- ncol(perm_data)
  if (nc == 1L) {
    return(perm_data)
  }
  dn <- dimnames(perm_data)
  pm <- matrix(as.numeric(perm_data), nrow = nr, ncol = nc, dimnames = dn)
  m <- pm[, j, drop = FALSE]
  structure(
    m,
    class = c("scanoneperm", "matrix"),
    method = attr(perm_data, "method"),
    model = attr(perm_data, "model"),
    type = attr(perm_data, "type")
  )
}

infer_trait_model <- function(x) {
  ## Determine whether to use a binary model for this trait.
  ## Rule: if all non-missing observed values are in {0, 1}, treat as binary.
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    suppressWarnings(xn <- as.numeric(x))
    if (!all(is.na(xn) == is.na(x))) x <- xn
  }

  x <- x[!is.na(x)]
  if (length(x) == 0) return("normal")

  ux <- sort(unique(x))
  if (all(ux %in% c(0, 1))) return("binary")
  "normal"
}

pick_best_interaction <- function(int_results, alpha = 0.05) {
  ## Goal: return an interaction term string like "Q1:Q2" or NULL.
  sm <- tryCatch(summary(int_results), error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  
  sm_df <- as.data.frame(sm)
  
  ## Identify columns for p-value and LOD, and optionally QTL indices
  pcol <- intersect(
    c("pval", "Pval", "p.value", "P.value", "Pvalue", "pvalue"),
    colnames(sm_df)
  )
  
  lodcol <- intersect(
    c("lod", "LOD"),
    colnames(sm_df)
  )
  
  q1col <- intersect(
    c("q1", "qtl1", "Q1", "i"),
    colnames(sm_df)
  )
  
  q2col <- intersect(
    c("q2", "qtl2", "Q2", "j"),
    colnames(sm_df)
  )
  
  if (length(pcol) < 1 || length(lodcol) < 1) return(NULL)
  
  pcol <- pcol[1]
  lodcol <- lodcol[1]
  
  sig <- sm_df[
    !is.na(sm_df[[pcol]]) & sm_df[[pcol]] <= alpha,
    ,
    drop = FALSE
  ]
  
  if (nrow(sig) == 0) return(NULL)
  
  best <- sig[which.max(sig[[lodcol]]), , drop = FALSE]
  
  ## Prefer explicit q1/q2 columns if they exist; otherwise parse rowname.
  term <- NULL
  
  if (length(q1col) >= 1 && length(q2col) >= 1) {
    q1 <- as.integer(best[[q1col[1]]])
    q2 <- as.integer(best[[q2col[1]]])
    
    if (!is.na(q1) && !is.na(q2)) {
      term <- paste0("Q", q1, ":Q", q2)
    }
  }
  
  if (is.null(term)) {
    rn <- rownames(best)
    
    if (!is.null(rn) && length(rn) >= 1) {
      r <- rn[1]
      
      ## Common formats: "Q1:Q2" or "1:2"
      if (grepl("^Q\\d+:Q\\d+$", r)) {
        term <- r
      }
      
      if (is.null(term) && grepl("^\\d+:\\d+$", r)) {
        parts <- strsplit(r, ":", fixed = TRUE)[[1]]
        term <- paste0("Q", parts[1], ":Q", parts[2])
      }
    }
  }
  
  term
}

summarize_fitqtl <- function(fit, alpha = 0.05) {
  out <- list()

  if (!is.null(fit$result.full)) out$result_full <- fit$result.full
  if (!is.null(fit$result.drop)) out$result_drop <- fit$result.drop
  if (!is.null(fit$ests) && !is.null(fit$ests$ests)) out$ests <- fit$ests$ests

  pve <- NA_real_
  if (!is.null(fit$result.full)) {
    cn <- colnames(fit$result.full)
    pve_col <- intersect(c("perc.var", "PVE", "pve", "%var", "PercentVar"), cn)
    if (length(pve_col) >= 1) pve <- suppressWarnings(as.numeric(fit$result.full[1, pve_col[1]]))
  }
  out$pve <- pve

  ## Significant terms from drop-one table
  sig <- NULL
  if (!is.null(fit$result.drop)) {
    pcol <- intersect(c("Pvalue", "pvalue", "p.val", "pval", "Pval"), colnames(fit$result.drop))
    if (length(pcol) >= 1) {
      pv <- fit$result.drop[[pcol[1]]]
      sig <- fit$result.drop[!is.na(pv) & pv <= alpha, , drop = FALSE]
    }
  }
  out$significant_terms <- sig
  out
}

calculate_covariate_only_pve_fitqtl <- function(cross,
                                                trait_col,
                                                qtl,
                                                covar_df,
                                                covar_name,
                                                method = "hk",
                                                alpha = 0.05) {
  
  covar_only_formula <- as.formula(
    paste("y ~", covar_name)
  )
  
  fit_covar_only <- fitqtl(
    cross,
    pheno.col = trait_col,
    qtl = qtl,
    method = method,
    get.ests = TRUE,
    covar = covar_df,
    formula = covar_only_formula
  )
  
  covar_only_sum <- summarize_fitqtl(
    fit_covar_only,
    alpha = alpha
  )
  
  list(
    fit = fit_covar_only,
    summary = covar_only_sum,
    pve = covar_only_sum$pve
  )
}

compute_intervals_by_chr <- function(scan_result, chrs, lodint_drops, bayes_prob = 0.95) {
  ## Returns a single data.frame of intervals for each chr.
  out <- data.frame(
    chr = character(),
    interval = character(),  # lod_1.5 | lod_2 | bayes
    marker = character(),
    pos = numeric(),
    lod = numeric(),
    stringsAsFactors = FALSE
  )

  for (chr in unique(chrs)) {
    chr <- as.character(chr)

    ## lodint intervals at specified drops
    for (nm in names(lodint_drops)) {
      drop <- as.numeric(lodint_drops[[nm]])
      x <- lodint(scan_result, chr = chr, drop = drop, expandtomarkers = TRUE)
      out <- rbind(out, data.frame(
        chr = paste0("chr", chr),
        interval = nm,
        marker = rownames(x),
        pos = x$pos,
        lod = x$lod,
        stringsAsFactors = FALSE
      ))
    }

    ## Bayesian interval
    bx <- bayesint(scan_result, chr = chr, prob = bayes_prob, expandtomarkers = TRUE)
    out <- rbind(out, data.frame(
      chr = paste0("chr", chr),
      interval = "bayes",
      marker = rownames(bx),
      pos = bx$pos,
      lod = bx$lod,
      stringsAsFactors = FALSE
    ))
  }
  out
}

## Manhattan plot helpers (optional; requires `scripts/adapted_qqman.R`)
mh_plot_multi <- function(chr_pos_lod, cutoff_lod) {
  cpl.df <- as.data.frame(chr_pos_lod)
  data_trim <- list()
  data_trim$CHR <- suppressWarnings(as.numeric(cpl.df$chr))
  data_trim$BP <- (cpl.df$pos) * 1e6
  data_trim$P <- cpl.df$lod
  data_trim <- as.data.frame(data_trim)
  data_trim <- na.omit(data_trim)
  manhattan_multi(
    data_trim,
    genomewideline = cutoff_lod,
    suggestiveline = -1,
    cex.axis = 1.5,
    cex.lab = 1.2
  )
}

plot_manhattan <- function(scanone_results, cutoff_lod, trait_name, output_path = plot_dir) {
  if (!file.exists(qqman_source)) {
    log_msg("Skipping Manhattan plot: missing ", qqman_source)
    return(invisible(FALSE))
  }

  ## load the qqman manhattan plot function
  source(qqman_source)

  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_path)) {
    log_msg("Skipping Manhattan plot: could not create plot dir: ", output_path)
    return(invisible(FALSE))
  }

  out_pdf <- file.path(output_path, paste0(trait_name, "_manhattan_plot.pdf"))
  ok <- tryCatch({
    pdf(out_pdf, width = 10, height = 6)
    TRUE
  }, error = function(e) {
    log_msg("Skipping Manhattan plot: cannot open PDF: ", out_pdf)
    log_msg("pdf() error: ", conditionMessage(e))
    FALSE
  })
  if (!isTRUE(ok)) return(invisible(FALSE))
  on.exit({
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(1, 1))
  mh_plot_multi(scanone_results, cutoff_lod)
  title(main = paste0(trait_name, " QTL: QTL map"))
  log_msg("Wrote Manhattan plot PDF: ", out_pdf)
  invisible(TRUE)
}

safe_dev_off <- function() {
  if (grDevices::dev.cur() > 1) grDevices::dev.off()
}

infer_b_allele_direction <- function(cross, trait_col, chr, pos) {
  ## Uses observed genotype calls at the nearest marker.
  marker <- tryCatch(find.marker(cross, chr = chr, pos = pos), warning = function(w) NA_character_, error = function(e) NA_character_)
  if (is.null(marker) || is.na(marker) || marker == "") {
    return(list(marker = NA_character_, direction = NA_character_, means = NULL))
  }

  g <- pull.geno(cross, chr)
  y <- pull.pheno(cross, trait_col)
  ok <- !is.na(g) & !is.na(y)
  if (!any(ok)) {
    return(list(marker = marker, direction = NA_character_, means = NULL))
  }

  ## Map genotype codes to labels using cross$geno
  gt_codes <- cross$geno[[as.character(chr)]]$data[, marker]
  ## If that fails for any reason, fall back to numeric codes
  if (is.null(gt_codes)) gt_codes <- g

  gt <- gt_codes[ok]
  yy <- y[ok]

  ## Genotype labels as observed (do not rename to AA/AB/BB — order may not match)
  gt_chr <- as.character(gt)
  means <- tapply(yy, gt_chr, mean, na.rm = TRUE)

  dir <- NA_character_
  if (all(c("AA", "BB") %in% names(means))) {
    if (means[["BB"]] > means[["AA"]]) dir <- "B allele increases trait (BB > AA)"
    if (means[["BB"]] < means[["AA"]]) dir <- "B allele decreases trait (BB < AA)"
    if (means[["BB"]] == means[["AA"]]) dir <- "No difference (BB == AA)"
  }

  list(marker = marker, direction = dir, means = means)
}

## ===========================
## Run (multi-trait)
## ===========================

## 1) Read cross
cross.data <- read.cross(
  "csvs",
  genfile = cross_genfile,
  phefile = cross_phefile,
  estimate.map = FALSE,
  genotypes = c("AA", "AB", "BB")
)
cat(paste(capture.output(summary(cross.data)), collapse = "\n"), "\n")

if (!safe_has_pheno(cross.data, covar_name)) {
  stop("Covariate column not found in phenotype data: ", covar_name)
}

## 2) Genotype probabilities
data_prob <- calc.genoprob(
  cross.data,
  step = genoprob_step,
  error.prob = genoprob_error_prob
)

## Covariate data.frame for fitqtl/addqtl/addint (column name must match covar_name)
covar_df <- stats::setNames(
  data.frame(pull.pheno(data_prob, covar_name), stringsAsFactors = FALSE),
  covar_name
)

## 3) Permutations (RDS only; build with scripts/generate_scanone_permutations.R)
if (max(perm_pheno_cols) > ncol(data_prob$pheno)) {
  stop(
    "perm_pheno_cols requests column ", max(perm_pheno_cols),
    " but only ", ncol(data_prob$pheno), " phenotype columns exist."
  )
}
perm_pheno_names <- colnames(data_prob$pheno)[perm_pheno_cols]

if (!file.exists(perm_file_in)) {
  stop("Permutation RDS not found: ", perm_file_in)
}
cat("[", .ts(), "] Loading permutations: ", perm_file_in, "\n", sep = "")
perm_data <- readRDS(perm_file_in)

perm_nc <- ncol(perm_data)
if (perm_nc != length(perm_pheno_names)) {
  stop(
    "Permutation columns (", perm_nc, ") != length(perm_pheno_cols) (",
    length(perm_pheno_names), "). Use an RDS built with the same perm_pheno_cols, ",
    "or regenerate with scripts/generate_scanone_permutations.R."
  )
}
pnc <- colnames(perm_data)
if (is.null(pnc)) {
  stop("perm_data has no colnames(); expected one column name per trait in perm_pheno_cols.")
}
if (!identical(as.character(pnc), as.character(perm_pheno_names))) {
  stop(
    "colnames(perm_data) do not match colnames(data_prob$pheno)[perm_pheno_cols]. ",
    "Regenerate the RDS or fix perm_pheno_cols / phenotype table order."
  )
}

cat("[", .ts(), "] Permutation matrix: ", nrow(perm_data), " x ", ncol(perm_data),
    " (traits=", length(perm_pheno_names), ")\n", sep = "")

cutoff_lod_0.10 <- summary(perm_data)[c("10%"), ]
cutoff_lod_0.05 <- summary(perm_data)[c("5%"), ]

## Determine which traits to run
traits_from_perm_cols <- perm_pheno_names
traits_from_perm_cols <- setdiff(traits_from_perm_cols, covar_name)
if (is.null(trait_cols)) {
  trait_cols <- traits_from_perm_cols
}

## One subfolder per calendar date × alpha (e.g. 20260513_alpha0.05). Re-running the same
## date and alpha replaces that folder; use a different alpha or change the clock only if you need both.
peak_intervals_run_id <- sprintf("%s_alpha%g", format(Sys.Date(), "%Y%m%d"), alpha)
peak_intervals_run_dir <- file.path(peak_intervals_dir, peak_intervals_run_id)
dir.create(peak_intervals_run_dir, recursive = TRUE, showWarnings = FALSE)
cat(
  "[", .ts(), "] Peak interval tables (this run): ", peak_intervals_run_dir, "\n",
  sep = ""
)

## Summary across traits (major results only)
all_trait_summary <- data.frame(
  trait = character(),
  model = character(),
  alpha = numeric(),
  n_qtl = integer(),
  best_chr = character(),
  best_pos = numeric(),
  best_lod = numeric(),
  selected_int = character(),
  pve = numeric(),
  covariate_only_pve = numeric(),
  total_qtl_pve = numeric(),
  trait_report_file = character(),
  peak_intervals_run_id = character(),
  peak_intervals_run_dir = character(),
  stringsAsFactors = FALSE
)

run_trait <- function(trait_col) {
  if (!safe_has_pheno(data_prob, trait_col)) {
    stop("Trait column not found in phenotype data: ", trait_col)
  }

  tl <- start_trait_log(trait_col)
  on.exit(tl$close(), add = TRUE)

  log_section("Trait run")
  trait_model <- infer_trait_model(data_prob$pheno[[trait_col]])
  log_msg("trait=", trait_col, ", model=", trait_model, " (auto), alpha=", alpha)

  ## Missingness summary
  n_total <- nind(data_prob)
  n_miss_trait <- sum(is.na(data_prob$pheno[[trait_col]]))
  n_miss_covar <- sum(is.na(covar_df[[covar_name]]))
  log_msg("Individuals: total=", n_total, ", missing ", trait_col, "=", n_miss_trait, ", missing ", covar_name, "=", n_miss_covar)

  cutoff_vec <- as.numeric(
    if (isTRUE(all.equal(alpha, 0.10))) cutoff_lod_0.10 else cutoff_lod_0.05
  )
  j <- match(as.character(trait_col), perm_pheno_names)
  if (is.na(j)) {
    stop("trait ", trait_col, " is not in perm_pheno_names (perm RDS column set).")
  }
  cutoff_trait <- cutoff_vec[j]
  log_msg("Permutation LOD cutoff at alpha=", alpha, ": ", signif(cutoff_trait, 4))

  ## 4) scanone
  log_section("scanone")
  scan_result <- scanone(
    data_prob,
    method = "em",
    pheno.col = trait_col,
    addcovar = covar_df[[covar_name]],
    model = trait_model
  )
  perm_this <- perm_column(perm_data, j)
  scan_summary <- summary(scan_result, perms = perm_this, alpha = alpha, pvalues = TRUE)
  log_print(scan_summary, "Significant peaks:")

  ## Always report intervals right after scanone (if any peaks)
  intervals_df <- NULL
  if (nrow(scan_summary) > 0) {
    interval_chrs <- unique(scan_summary$chr)
    intervals_df <- compute_intervals_by_chr(
      scan_result = scan_result,
      chrs = interval_chrs,
      lodint_drops = lodint_drops,
      bayes_prob = bayes_prob
    )
    log_print(intervals_df, "Intervals by chromosome (from scanone output):")

    intervals_file <- file.path(peak_intervals_run_dir, paste0(trait_col, "_peak_intervals.tsv"))
    write.table(intervals_df, file = intervals_file, sep = "\t", quote = FALSE, row.names = FALSE)
    log_msg("Wrote intervals table to: ", intervals_file)
  } else {
    log_msg("No significant peaks at alpha=", alpha, ". Ending trait pipeline early.")
  }

  ## Defaults for summary row
  best_chr <- NA_character_
  best_pos <- NA_real_
  best_lod <- NA_real_
  n_qtl <- 0L
  selected_int <- NA_character_
  pve <- NA_real_
  covariate_only_pve <- NA_real_
  total_qtl_pve <- NA_real_

  if (nrow(scan_summary) > 0) {
    ## Best peak from scanone
    best_i <- which.max(scan_summary$lod)
    best_chr <- as.character(scan_summary$chr[best_i])
    best_pos <- as.numeric(scan_summary$pos[best_i])
    best_lod <- as.numeric(scan_summary$lod[best_i])

    ## 5) initial QTL object
    log_section("makeqtl + fitqtl")
    qtl <- makeqtl(data_prob, chr = scan_summary$chr, pos = scan_summary$pos, what = "prob")
    n_qtl <- nqtl(qtl)
    log_msg("Initial QTL count: ", n_qtl)

    ## 6) additive model
    current_formula <- make_add_formula(n_qtl, covar_name = covar_name)
    log_msg("Additive model formula: ", deparse(current_formula))
    fit_add <- fitqtl(
      data_prob, pheno.col = trait_col, qtl = qtl, method = "hk",
      get.ests = TRUE, covar = covar_df, formula = current_formula
    )
    add_sum <- summarize_fitqtl(fit_add, alpha = alpha)
    if (!is.na(add_sum$pve)) {
      pve <- add_sum$pve
      log_msg("Additive model PVE (if available): ", signif(pve, 4))
    }
    if (!is.null(add_sum$significant_terms) && nrow(add_sum$significant_terms) > 0) {
      log_print(add_sum$significant_terms, "Significant terms (drop-one table, if available):")
    }
    if (!is.null(add_sum$result_drop)) log_print(add_sum$result_drop, "Drop-one table:")

    ## 7) test pairwise interactions
    log_section("Interactions (addint)")
    selected_int <- NULL
    if (n_qtl > 1) {
      int_results <- addint(data_prob, pheno.col = trait_col, qtl = qtl, method = "hk", covar = covar_df)
      log_print(int_results, "addint results:")
      selected_int <- pick_best_interaction(int_results, alpha = alpha)
      if (!is.null(selected_int)) log_msg("Selected interaction: ", selected_int)
      if (is.null(selected_int)) log_msg("No significant interaction selected at alpha=", alpha)
    } else {
      log_msg("Only one QTL; skipping addint.")
    }

    ## 8) refit with interaction if selected
    log_section("Refit (+interaction if selected)")
    current_formula <- make_formula_with_int(n_qtl, interactions = selected_int, covar_name = covar_name)
    log_msg("Final model formula: ", deparse(current_formula))
    fit_add2 <- fitqtl(
      data_prob, pheno.col = trait_col, qtl = qtl, method = "hk",
      get.ests = TRUE, covar = covar_df, formula = current_formula
    )
    add2_sum <- summarize_fitqtl(fit_add2, alpha = alpha)
    if (!is.na(add2_sum$pve)) log_msg("Final model PVE (if available): ", signif(add2_sum$pve, 4))
    if (!is.null(add2_sum$significant_terms) && nrow(add2_sum$significant_terms) > 0) {
      log_print(add2_sum$significant_terms, "Significant terms (drop-one table, if available):")
    }

    ## 9) refine positions
    log_section("refineqtl + refit")
    rqtl <- refineqtl(
      data_prob, pheno.col = trait_col, qtl = qtl, method = "hk",
      covar = covar_df, formula = current_formula
    )
    log_print(rqtl, "Refined QTL object:")

    fit_refined <- fitqtl(
      data_prob, pheno.col = trait_col, qtl = rqtl, method = "hk",
      get.ests = TRUE, covar = covar_df, formula = current_formula
    )
    refined_sum <- summarize_fitqtl(fit_refined, alpha = alpha)
    
    if (!is.na(refined_sum$pve)) {
      pve <- refined_sum$pve
      log_msg("Refined full model PVE, covariate + all QTL: ", signif(pve, 4))
    }
    
    ## Calculate covariate-only PVE using fitqtl()
    log_section("Total QTL PVE")
    
    covar_only <- calculate_covariate_only_pve_fitqtl(
      cross = data_prob,
      trait_col = trait_col,
      qtl = rqtl,
      covar_df = covar_df,
      covar_name = covar_name,
      method = "hk",
      alpha = alpha
    )
    
    covariate_only_pve <- covar_only$pve
    
    if (!is.na(covariate_only_pve)) {
      log_msg("Covariate-only model formula: y ~ ", covar_name)
      log_msg("Covariate-only PVE from fitqtl: ", signif(covariate_only_pve, 4))
      
      if (!is.null(covar_only$summary$result_full)) {
        log_print(covar_only$summary$result_full, "Covariate-only fitqtl full model table:")
      }
    }
    
    if (!is.na(pve) && !is.na(covariate_only_pve)) {
      total_qtl_pve <- pve - covariate_only_pve
      
      log_msg(
        "Total phenotypic variance explained by all QTL ",
        "(refined full model PVE - covariate-only PVE): ",
        signif(total_qtl_pve, 4)
      )
    } else {
      log_msg("Could not calculate total QTL PVE because refined full model PVE or covariate-only PVE is NA.")
    }
    
    if (!is.null(refined_sum$significant_terms) && nrow(refined_sum$significant_terms) > 0) {
      log_print(refined_sum$significant_terms, "Significant terms (drop-one table, if available):")
    }
    if (!is.null(refined_sum$result_drop)) log_print(refined_sum$result_drop, "Drop-one table:")
    if (!is.null(refined_sum$ests)) log_print(refined_sum$ests, "Effect estimates:")

    ## 10) add one more QTL
    log_section("addqtl")
    out_aq <- addqtl(
      data_prob, pheno.col = trait_col, qtl = rqtl, method = "hk",
      covar = covar_df, formula = current_formula
    )
    add_summary <- summary(out_aq, perms = perm_this, alpha = alpha, pvalues = TRUE)
    log_print(add_summary, "Candidate added QTL:")

    ## 11) plots -> write PDFs
    log_section("plots")

    ## Manhattan plot (optional)
    plot_manhattan(scan_result, cutoff_lod = cutoff_trait, trait_name = trait_col, output_path = plot_dir)

    ## Combined PDF for scanone-by-chr + effectplots
    trait_plot_pdf <- file.path(plot_dir, paste0(trait_col, "_qtl_plots.pdf"))
    ok_pdf <- tryCatch({
      pdf(trait_plot_pdf, width = 11, height = 8.5)
      TRUE
    }, error = function(e) {
      log_msg("Skipping trait plot PDF: cannot open PDF: ", trait_plot_pdf)
      log_msg("pdf() error: ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok_pdf)) on.exit(safe_dev_off(), add = TRUE)

    available_chrs <- names(data_prob$geno)
    refined_chrs <- unique(as.character(rqtl$chr))
    if (isTRUE(ok_pdf)) {
      for (chr in refined_chrs) {
        if (!(chr %in% available_chrs)) {
          log_msg("Skipping scan plot: chr not found in cross: ", chr)
          next
        }
        plot(scan_result, chr = chr, main = paste0(trait_col, " - Chr", chr, " (alpha=", alpha, ")"))
        abline(h = cutoff_trait, col = "red", lwd = 3)
        idx <- which(as.character(rqtl$chr) == chr)
        if (length(idx) > 0) abline(v = rqtl$pos[idx], col = "green", lwd = 3)
      }
    }

    data_sim <- sim.geno(data_prob, step = 1, n.draws = 16)
    if (isTRUE(ok_pdf)) {
      for (i in seq_along(rqtl$chr)) {
        chr <- as.character(rqtl$chr[i])
        pos <- rqtl$pos[i]
        if (!(chr %in% available_chrs)) {
          log_msg("Skipping effectplot: chr not found in cross: ", chr)
          next
        }
        marker <- tryCatch(find.marker(data_prob, chr = chr, pos = pos), warning = function(w) NA_character_, error = function(e) NA_character_)
        if (is.na(marker) || marker == "") {
          log_msg("Skipping effectplot: could not find marker for chr", chr, "@", pos)
          next
        }
        effectplot(
          data_sim,
          pheno.col = trait_col,
          mname1 = marker,
          main = paste0(trait_col, " effect: chr", chr, "@", pos)
        )
      }

      safe_dev_off()
      log_msg("Wrote trait plots PDF: ", trait_plot_pdf)
    }

    ## Log B-allele direction at refined QTL (nearest marker)
    log_section("B-allele direction (nearest marker)")
    for (i in seq_along(rqtl$chr)) {
      chr <- as.character(rqtl$chr[i])
      pos <- rqtl$pos[i]
      if (!(chr %in% available_chrs)) {
        log_msg("Skipping B-allele direction: chr not found in cross: ", chr)
        next
      }
      dir_info <- infer_b_allele_direction(data_prob, trait_col, chr = chr, pos = pos)
      log_msg("QTL ", i, ": chr", chr, "@", pos, " marker=", dir_info$marker)
      if (!is.null(dir_info$means)) {
        log_print(dir_info$means, "Genotype means:")
      }
      log_msg("Direction: ", dir_info$direction)
    }
  }

  log_msg("Trait completed.")

  list(
    trait = trait_col,
    model = trait_model,
    alpha = alpha,
    n_qtl = n_qtl,
    best_chr = best_chr,
    best_pos = best_pos,
    best_lod = best_lod,
    selected_int = if (is.null(selected_int)) NA_character_ else selected_int,
    pve = pve,
    covariate_only_pve = covariate_only_pve,
    total_qtl_pve = total_qtl_pve,
    trait_report_file = tl$file,
    peak_intervals_run_id = peak_intervals_run_id,
    peak_intervals_run_dir = peak_intervals_run_dir
  )
}

for (tr in trait_cols) {
  cat("[", .ts(), "] Running trait: ", tr, "\n", sep = "")
  row <- run_trait(tr)
  all_trait_summary <- rbind(all_trait_summary, as.data.frame(row, stringsAsFactors = FALSE))
}

summary_file <- file.path(getwd(), paste0("qtl_trait_summary_alpha", alpha, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"))
write.table(all_trait_summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[", .ts(), "] Wrote cross-trait summary to: ", summary_file, "\n", sep = "")
