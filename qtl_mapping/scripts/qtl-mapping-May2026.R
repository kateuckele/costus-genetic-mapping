## ============================================================
## QTL mapping pipeline (R/qtl): scanone, fitqtl, plots, outputs
## Outputs: plots → results/plots; runtime tee logs → results/logs; per-trait transcripts →
##   results/trait_qtl_reports/*.txt; peak interval tables →
##   results/peak_intervals/<YYYYMMDD_alpha*/*.tsv (one subfolder per date×alpha; same day+alpha overwrites);
##   merged per-QTL table → results/qtl_lod15_summary.tsv (from in-memory fitqtl / refineqtl / lodint, not log parsing)
## Permutations: run scripts/generate_scanone_permutations.R first; set perm_dir_in below.
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
## Per-trait sink() transcripts (scanone, fitqtl, etc.)
trait_report_dir <- file.path(results_dir, "trait_qtl_reports")
peak_intervals_dir <- file.path(results_dir, "peak_intervals")
for (d in c(results_dir, plot_dir, log_dir, trait_report_dir, peak_intervals_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}
setwd(results_dir)

## ----------------------- Configuration -----------------------
## NULL = all traits in the permutation manifest except covar_name; or c("INFA", ...).
trait_cols <- NULL

alpha <- 0.05

genoprob_step <- 0.5
genoprob_error_prob <- 0.001

cross_genfile <- file.path(qtl_root, "data", "mapthis_LG.csv")
cross_phefile <- file.path(qtl_root, "data", "costus_pheno_rqtl_2025Jan24.csv")

covar_name <- "F1parent"

## Must match generate_scanone_permutations.R (same perm_pheno_cols, per-trait model,
## and addcovar = F1parent used when building each per-trait RDS).
perm_dir_in <- file.path(
  processed_dir,
  "scanone_1000perm_by_trait_20260608_134657"
)
perm_manifest_in <- file.path(perm_dir_in, "manifest.tsv")
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
make_qtl_formula <- function(n_qtl, interactions = NULL, covar_name = "F1parent") {
  q_terms <- paste0("Q", seq_len(n_qtl))
  rhs <- c(q_terms, interactions, covar_name)
  as.formula(paste("y ~", paste(rhs, collapse = " + ")))
}

safe_has_pheno <- function(cross, pheno_name) {
  !is.null(cross$pheno) && (pheno_name %in% colnames(cross$pheno))
}

resolve_perm_file <- function(file_path, manifest_dir) {
  if (file.exists(file_path)) return(file_path)
  alt <- file.path(manifest_dir, file_path)
  if (file.exists(alt)) return(alt)
  file_path
}

permutation_cutoff <- function(perm_obj, alpha) {
  cutoff_name <- paste0(format(100 * alpha, trim = TRUE, scientific = FALSE), "%")
  sm <- summary(perm_obj)
  if (is.null(dim(sm))) {
    if (!(cutoff_name %in% names(sm))) {
      stop("Permutation summary has no cutoff named ", cutoff_name, ".")
    }
    return(as.numeric(sm[[cutoff_name]]))
  }
  if (!(cutoff_name %in% rownames(sm))) {
    stop("Permutation summary has no cutoff row named ", cutoff_name, ".")
  }
  as.numeric(sm[cutoff_name, 1, drop = TRUE])
}

interaction_formula_term <- function(label, qtl_names = NULL, covar_name = NULL) {
  ## Convert addint row labels such as "2@71.7:4.2.9@20.7" to "Q1:Q2".
  if (is.null(label) || is.na(label) || !nzchar(label)) return(NULL)
  
  if (grepl("^Q\\d+:Q\\d+$", label)) return(label)
  
  if (grepl("^\\d+:\\d+$", label)) {
    parts <- strsplit(label, ":", fixed = TRUE)[[1]]
    return(paste0("Q", parts[1], ":Q", parts[2]))
  }
  
  parts <- strsplit(label, ":", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NULL)
  
  mapped <- parts
  if (!is.null(qtl_names)) {
    qtl_names <- as.character(qtl_names)
    for (i in seq_along(parts)) {
      j <- match(parts[i], qtl_names)
      if (!is.na(j)) mapped[i] <- paste0("Q", j)
    }
  }
  if (!is.null(covar_name)) {
    mapped[mapped == covar_name] <- covar_name
  }
  
  covar_allowed <- rep(FALSE, length(mapped))
  if (!is.null(covar_name)) {
    covar_allowed <- mapped == covar_name
  }
  if (all(grepl("^Q\\d+$", mapped) | covar_allowed)) {
    return(paste(mapped, collapse = ":"))
  }
  
  NULL
}

summarize_addint <- function(int_results, alpha = 0.05, qtl_names = NULL, covar_name = NULL) {
  sm <- tryCatch(summary(int_results), error = function(e) NULL)
  if (is.null(sm)) {
    return(list(
      table = NULL,
      significant_labels = character(),
      formula_terms = character(),
      perc_var = numeric()
    ))
  }
  
  sm_df <- as.data.frame(sm)
  pcol <- "Pvalue(F)"
  vcol <- "%var"
  if (!(pcol %in% colnames(sm_df))) {
    stop("addint summary is missing expected interaction p-value column: ", pcol)
  }
  if (!(vcol %in% colnames(sm_df))) {
    stop("addint summary is missing expected interaction variance column: ", vcol)
  }
  
  sig <- sm_df[
    !is.na(sm_df[[pcol]]) & sm_df[[pcol]] <= alpha,
    ,
    drop = FALSE
  ]
  if (nrow(sig) == 0) {
    return(list(
      table = sm_df,
      significant_labels = character(),
      formula_terms = character(),
      perc_var = numeric()
    ))
  }
  
  labels <- rownames(sig)
  mapped_terms <- character(nrow(sig))
  for (i in seq_len(nrow(sig))) {
    mapped_terms[i] <- interaction_formula_term(
      labels[i],
      qtl_names = qtl_names,
      covar_name = covar_name
    )
  }
  mapped_ok <- !is.na(mapped_terms) & nzchar(mapped_terms)
  formula_terms <- unique(mapped_terms[mapped_ok])
  
  perc_var <- suppressWarnings(as.numeric(sig[[vcol]]))
  
  list(
    table = sm_df,
    significant_labels = labels,
    formula_terms = formula_terms,
    perc_var = perc_var
  )
}

summarize_fitqtl <- function(fit, alpha = 0.05) {
  out <- list()

  if (!is.null(fit$result.full)) out$result_full <- fit$result.full
  if (!is.null(fit$result.drop)) out$result_drop <- fit$result.drop
  if (!is.null(fit$ests) && !is.null(fit$ests$ests)) out$ests <- fit$ests$ests

  pve <- NA_real_
  if (!is.null(fit$result.full) && "%var" %in% colnames(fit$result.full)) {
    pve <- suppressWarnings(as.numeric(fit$result.full[1, "%var"]))
  }
  out$pve <- pve

  sig <- NULL
  if (!is.null(fit$result.drop) && "Pvalue(F)" %in% colnames(fit$result.drop)) {
    pv <- fit$result.drop[, "Pvalue(F)", drop = TRUE]
    sig <- fit$result.drop[!is.na(pv) & pv <= alpha, , drop = FALSE]
  }
  out$significant_terms <- sig
  out
}

log_fit_summary <- function(fit_sum,
                            pve_label = NULL,
                            log_significant = TRUE,
                            log_drop = TRUE,
                            log_ests = FALSE) {
  if (!is.null(pve_label) && !is.na(fit_sum$pve)) {
    log_msg(pve_label, signif(fit_sum$pve, 4))
  }
  if (isTRUE(log_significant) &&
      !is.null(fit_sum$significant_terms) &&
      nrow(fit_sum$significant_terms) > 0) {
    log_print(fit_sum$significant_terms, "Significant terms (drop-one table, if available):")
  }
  if (isTRUE(log_drop) && !is.null(fit_sum$result_drop)) {
    log_print(fit_sum$result_drop, "Drop-one table:")
  }
  if (isTRUE(log_ests) && !is.null(fit_sum$ests)) {
    log_print(fit_sum$ests, "Effect estimates:")
  }
}

## ---- Per-QTL summary table (lod_1.5 row for karyoplote / tables; built during run_trait) ----

empty_qtl_lod15_rows <- function() {
  data.frame(
    trait = character(),
    qtl_id = character(),
    chr = character(),
    interval_type = character(),
    start_marker = character(),
    start_pos = numeric(),
    start_lod = numeric(),
    peak_marker = character(),
    peak_pos = numeric(),
    peak_lod = numeric(),
    end_marker = character(),
    end_pos = numeric(),
    end_lod = numeric(),
    refined_peak_pos = character(),
    additive = character(),
    dominance = character(),
    B_allele_direction = character(),
    qtl_pve_percent = numeric(),
    trait_report_file = character(),
    interval_file = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

## Label matching fitqtl drop-one rownames (e.g. "2@46.1")
qtl_fit_term_label <- function(qobj, i) {
  if (!is.null(qobj$name) && length(qobj$name) >= i && nzchar(as.character(qobj$name[i]))) {
    return(as.character(qobj$name[i]))
  }
  chr <- as.character(qobj$chr[i])
  pos <- as.numeric(qobj$pos[i])
  paste0(chr, "@", sprintf("%.1f", pos))
}

extract_ad_dom_from_fitqtl <- function(fit, chr, pos) {
  add <- ""
  dom <- ""
  if (is.null(fit) || is.null(fit$ests)) {
    return(list(add = add, dom = dom))
  }
  est <- fit$ests$ests
  if (is.null(est)) {
    return(list(add = add, dom = dom))
  }
  chr <- as.character(chr)
  pos <- as.numeric(pos)
  ra <- paste0(chr, "@", sprintf("%.1f", pos), "a")
  rd <- paste0(chr, "@", sprintf("%.1f", pos), "d")
  if (is.matrix(est)) {
    rn <- rownames(est)
    if (!is.null(rn)) {
      if (ra %in% rn) {
        add <- paste(est[ra, , drop = TRUE], collapse = " ")
      }
      if (rd %in% rn) {
        dom <- paste(est[rd, , drop = TRUE], collapse = " ")
      }
    }
  } else if (length(est)) {
    nm <- names(est)
    if (!is.null(nm)) {
      if (ra %in% nm) {
        add <- as.character(est[[ra]])
      }
      if (rd %in% nm) {
        dom <- as.character(est[[rd]])
      }
    }
  }
  list(add = add, dom = dom)
}

extract_pve_from_result_drop <- function(result_drop, term) {
  if (is.null(result_drop)) {
    return(NA_real_)
  }
  if (!nrow(result_drop)) {
    return(NA_real_)
  }
  if (!nzchar(term)) {
    return(NA_real_)
  }
  pc <- "%var"
  if (!(pc %in% colnames(result_drop))) {
    return(NA_real_)
  }
  rn <- rownames(result_drop)
  if (!is.null(rn) && term %in% rn) {
    return(suppressWarnings(as.numeric(result_drop[term, pc, drop = TRUE])))
  }
  NA_real_
}

## One row per refined QTL; interval endpoints from scanone lodint (same interval_type as peak TSVs).
build_qtl_lod15_detail_rows <- function(trait,
                                         interval_type,
                                         intervals_df,
                                         b_allele_info,
                                         rqtl,
                                         fit_refined,
                                         refined_sum,
                                         interval_file,
                                         trait_report_file) {
  out <- empty_qtl_lod15_rows()
  if (is.null(intervals_df) || !nrow(intervals_df)) {
    return(out)
  }
  if (is.null(rqtl) || !length(rqtl$chr)) {
    return(out)
  }
  rd <- refined_sum$result_drop
  rows <- list()
  n <- length(rqtl$chr)
  for (i in seq_len(n)) {
    chr <- as.character(rqtl$chr[i])
    chr_key <- paste0("chr", chr)
    pos <- as.numeric(rqtl$pos[i])
    term <- qtl_fit_term_label(rqtl, i)
    sub <- intervals_df[
      as.character(intervals_df$chr) == chr_key &
        as.character(intervals_df$interval) == interval_type,
      ,
      drop = FALSE
    ]
    notes <- character()
    if (nrow(sub) < 3L) {
      notes <- c(
        notes,
        sprintf("expected_3_%s_rows_chr_%s_got_%d", interval_type, chr_key, nrow(sub))
      )
      next
    }
    o <- order(as.numeric(sub$pos))
    sub <- sub[o, , drop = FALSE]
    i0 <- 1L
    i1 <- nrow(sub)
    ip <- which.max(as.numeric(sub$lod))
    st <- sub[i0, , drop = FALSE]
    ed <- sub[i1, , drop = FALSE]
    pk <- sub[ip, , drop = FALSE]
    ad <- extract_ad_dom_from_fitqtl(fit_refined, chr, pos)
    pve <- extract_pve_from_result_drop(rd, term)
    dir_info <- b_allele_info[[i]]
    b_allele_direction <- if (is.na(dir_info$direction)) "" else dir_info$direction
    rows[[length(rows) + 1L]] <- data.frame(
      trait = trait,
      qtl_id = paste0("Q", i),
      chr = chr_key,
      interval_type = interval_type,
      start_marker = as.character(st$marker),
      start_pos = as.numeric(st$pos),
      start_lod = as.numeric(st$lod),
      peak_marker = as.character(pk$marker),
      peak_pos = as.numeric(pk$pos),
      peak_lod = as.numeric(pk$lod),
      end_marker = as.character(ed$marker),
      end_pos = as.numeric(ed$pos),
      end_lod = as.numeric(ed$lod),
      refined_peak_pos = as.character(pos),
      additive = ad$add,
      dominance = ad$dom,
      B_allele_direction = b_allele_direction,
      qtl_pve_percent = pve,
      trait_report_file = trait_report_file,
      interval_file = interval_file,
      notes = paste(notes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  if (!length(rows)) {
    return(out)
  }
  do.call(rbind, rows)
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

  gt_codes <- cross$geno[[as.character(chr)]]$data[, marker]
  y <- pull.pheno(cross, trait_col)
  ok <- !is.na(gt_codes) & !is.na(y)
  if (!any(ok)) {
    return(list(marker = marker, direction = NA_character_, means = NULL))
  }

  gt <- gt_codes[ok]
  yy <- y[ok]

  gt_chr <- as.character(gt)
  if (all(stats::na.omit(gt_chr) %in% c("1", "2", "3"))) {
    gt_chr <- c(`1` = "AA", `2` = "AB", `3` = "BB")[gt_chr]
  }
  means <- tapply(yy, gt_chr, mean, na.rm = TRUE)

  dir <- NA_character_
  if (all(c("AA", "BB") %in% names(means))) {
    if (means[["BB"]] > means[["AA"]]) dir <- "B allele increases trait (BB > AA)"
    if (means[["BB"]] < means[["AA"]]) dir <- "B allele decreases trait (BB < AA)"
    if (means[["BB"]] == means[["AA"]]) dir <- "No difference (BB == AA)"
  }

  list(marker = marker, direction = dir, means = means)
}

collect_b_allele_info <- function(cross, trait_col, rqtl, available_chrs = names(cross$geno)) {
  out <- vector("list", length(rqtl$chr))
  for (i in seq_along(rqtl$chr)) {
    chr <- as.character(rqtl$chr[i])
    pos <- rqtl$pos[i]
    if (!(chr %in% available_chrs)) {
      out[[i]] <- list(
        chr = chr,
        pos = pos,
        marker = NA_character_,
        direction = NA_character_,
        means = NULL,
        skipped = TRUE
      )
    } else {
      dir_info <- infer_b_allele_direction(cross, trait_col, chr = chr, pos = pos)
      out[[i]] <- c(list(chr = chr, pos = pos, skipped = FALSE), dir_info)
    }
  }
  out
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

## 3) Permutations (one RDS per trait; build with scripts/generate_scanone_permutations.R)
if (max(perm_pheno_cols) > ncol(data_prob$pheno)) {
  stop(
    "perm_pheno_cols requests column ", max(perm_pheno_cols),
    " but only ", ncol(data_prob$pheno), " phenotype columns exist."
  )
}
perm_pheno_names <- colnames(data_prob$pheno)[perm_pheno_cols]

if (!file.exists(perm_manifest_in)) {
  stop("Permutation manifest not found: ", perm_manifest_in)
}
cat("[", .ts(), "] Loading permutation manifest: ", perm_manifest_in, "\n", sep = "")
perm_manifest <- read.table(
  perm_manifest_in,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_perm_cols <- c("trait", "perm_pheno_col", "model", "file")
missing_perm_cols <- setdiff(required_perm_cols, colnames(perm_manifest))
if (length(missing_perm_cols) > 0) {
  stop("Permutation manifest missing columns: ", paste(missing_perm_cols, collapse = ", "))
}
if (!identical(as.character(perm_manifest$trait), as.character(perm_pheno_names))) {
  stop(
    "Permutation manifest traits do not match colnames(data_prob$pheno)[perm_pheno_cols]. ",
    "Regenerate permutations or fix perm_dir_in / perm_pheno_cols."
  )
}
if (!identical(as.integer(perm_manifest$perm_pheno_col), as.integer(perm_pheno_cols))) {
  stop(
    "Permutation manifest perm_pheno_col values do not match perm_pheno_cols. ",
    "Regenerate permutations or fix perm_dir_in / perm_pheno_cols."
  )
}
perm_manifest$file <- vapply(
  perm_manifest$file,
  resolve_perm_file,
  character(1),
  manifest_dir = dirname(perm_manifest_in)
)
missing_perm_files <- perm_manifest$file[!file.exists(perm_manifest$file)]
if (length(missing_perm_files) > 0) {
  stop("Permutation RDS file(s) not found:\n  ", paste(missing_perm_files, collapse = "\n  "))
}
cat("[", .ts(), "] Per-trait permutation files: ", nrow(perm_manifest), "\n", sep = "")

## Determine which traits to run
traits_from_perm_cols <- perm_manifest$trait
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
  significant_ints = character(),
  significant_int_pve = character(),
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

  j <- match(as.character(trait_col), perm_manifest$trait)
  if (is.na(j)) {
    stop("trait ", trait_col, " is not in the permutation manifest.")
  }
  perm_row <- perm_manifest[j, , drop = FALSE]
  trait_model <- as.character(perm_row$model)
  
  log_section("Trait run")
  log_msg("trait=", trait_col, ", model=", trait_model, " (from permutation manifest), alpha=", alpha)

  ## Missingness summary
  n_total <- nind(data_prob)
  n_miss_trait <- sum(is.na(data_prob$pheno[[trait_col]]))
  n_miss_covar <- sum(is.na(covar_df[[covar_name]]))
  log_msg("Individuals: total=", n_total, ", missing ", trait_col, "=", n_miss_trait, ", missing ", covar_name, "=", n_miss_covar)

  log_msg("Loading trait permutations: ", perm_row$file)
  perm_this <- readRDS(perm_row$file)
  cutoff_trait <- permutation_cutoff(perm_this, alpha = alpha)
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
  n_qtl <- 0L
  significant_ints <- NA_character_
  significant_int_pve <- NA_character_
  interaction_terms <- NULL
  pve <- NA_real_
  covariate_only_pve <- NA_real_
  total_qtl_pve <- NA_real_
  qtl_detail_rows <- NULL

  if (nrow(scan_summary) > 0) {
    ## 5) initial QTL object
    log_section("makeqtl + fitqtl")
    qtl <- makeqtl(data_prob, chr = scan_summary$chr, pos = scan_summary$pos, what = "prob")
    n_qtl <- nqtl(qtl)
    log_msg("Initial QTL count: ", n_qtl)

    ## 6) additive model
    current_formula <- make_qtl_formula(n_qtl, covar_name = covar_name)
    log_msg("Additive model formula: ", deparse(current_formula))
    fit_add <- fitqtl(
      data_prob, pheno.col = trait_col, qtl = qtl, method = "hk",
      get.ests = TRUE, covar = covar_df, formula = current_formula
    )
    add_sum <- summarize_fitqtl(fit_add, alpha = alpha)
    if (!is.na(add_sum$pve)) {
      pve <- add_sum$pve
    }
    log_fit_summary(add_sum, pve_label = "Additive model PVE (if available): ")

    ## 7) test pairwise interactions
    log_section("Interactions (addint)")
    int_results <- addint(data_prob, pheno.col = trait_col, qtl = qtl, method = "hk", covar = covar_df)
    log_print(int_results, "addint results:")
    int_sum <- summarize_addint(int_results, alpha = alpha, qtl_names = qtl$name, covar_name = covar_name)
    if (length(int_sum$significant_labels) > 0) {
      significant_ints <- paste(int_sum$significant_labels, collapse = ", ")
      significant_int_pve <- paste(signif(int_sum$perc_var, 4), collapse = ", ")
      interaction_terms <- int_sum$formula_terms
      log_msg("Significant interactions at alpha=", alpha, ": ", significant_ints)
      log_msg("Significant interaction % variance explained: ", significant_int_pve)
      log_msg("Interaction terms included in final model: ", paste(interaction_terms, collapse = ", "))
    } else {
      log_msg("No significant interactions at alpha=", alpha)
    }

    ## 8) refit with all significant interactions
    log_section("Refit (+significant interactions)")
    current_formula <- make_qtl_formula(n_qtl, interactions = interaction_terms, covar_name = covar_name)
    log_msg("Final model formula: ", deparse(current_formula))
    fit_add2 <- fitqtl(
      data_prob, pheno.col = trait_col, qtl = qtl, method = "hk",
      get.ests = TRUE, covar = covar_df, formula = current_formula
    )
    add2_sum <- summarize_fitqtl(fit_add2, alpha = alpha)
    log_fit_summary(add2_sum, pve_label = "Final model PVE (if available): ", log_drop = FALSE)

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
    }
    log_fit_summary(
      refined_sum,
      pve_label = "Refined full model PVE, covariate + all QTL: ",
      log_significant = FALSE,
      log_drop = FALSE
    )
    
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
    
    log_fit_summary(refined_sum, pve_label = NULL, log_ests = TRUE)

    available_chrs <- names(data_prob$geno)
    b_allele_info <- collect_b_allele_info(
      data_prob,
      trait_col,
      rqtl,
      available_chrs = available_chrs
    )

    qtl_detail_rows <- build_qtl_lod15_detail_rows(
      trait = trait_col,
      interval_type = "lod_1.5",
      intervals_df = intervals_df,
      b_allele_info = b_allele_info,
      rqtl = rqtl,
      fit_refined = fit_refined,
      refined_sum = refined_sum,
      interval_file = basename(intervals_file),
      trait_report_file = basename(tl$file)
    )

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
    for (i in seq_along(b_allele_info)) {
      dir_info <- b_allele_info[[i]]
      if (isTRUE(dir_info$skipped)) {
        log_msg("Skipping B-allele direction: chr not found in cross: ", dir_info$chr)
        next
      }
      log_msg("QTL ", i, ": chr", dir_info$chr, "@", dir_info$pos, " marker=", dir_info$marker)
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
    significant_ints = significant_ints,
    significant_int_pve = significant_int_pve,
    pve = pve,
    covariate_only_pve = covariate_only_pve,
    total_qtl_pve = total_qtl_pve,
    trait_report_file = tl$file,
    peak_intervals_run_id = peak_intervals_run_id,
    peak_intervals_run_dir = peak_intervals_run_dir,
    qtl_detail_rows = qtl_detail_rows
  )
}

all_qtl_lod15 <- empty_qtl_lod15_rows()

for (tr in trait_cols) {
  cat("[", .ts(), "] Running trait: ", tr, "\n", sep = "")
  rout <- run_trait(tr)
  all_trait_summary <- rbind(
    all_trait_summary,
    as.data.frame(rout[names(rout) != "qtl_detail_rows"], stringsAsFactors = FALSE)
  )
  if (!is.null(rout$qtl_detail_rows) && nrow(rout$qtl_detail_rows) > 0) {
    all_qtl_lod15 <- rbind(all_qtl_lod15, rout$qtl_detail_rows)
  }
}

summary_file <- file.path(getwd(), paste0("qtl_trait_summary_alpha", alpha, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tsv"))
write.table(all_trait_summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[", .ts(), "] Wrote cross-trait summary to: ", summary_file, "\n", sep = "")

lod15_file <- file.path(getwd(), "qtl_lod15_summary.tsv")
if (nrow(all_qtl_lod15) > 0) {
  write.table(all_qtl_lod15, file = lod15_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("[", .ts(), "] Wrote per-QTL lod summary (", nrow(all_qtl_lod15), " rows): ", lod15_file, "\n", sep = "")
} else {
  cat("[", .ts(), "] No qtl_lod15_summary.tsv (no QTL rows accumulated).\n", sep = "")
}
