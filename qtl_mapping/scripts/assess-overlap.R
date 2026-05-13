## QTL peak-interval overlap (trait x trait x chromosome)
##
## Reads *_peak_intervals.tsv files (from qtl-mapping-May2026.R), filters to one
## interval type (e.g. lod_1.5, lod_2, bayes), and writes pairwise overlap metrics
## plus peak_distance_cM = |peak_pos_trait1 - peak_pos_trait2| (pos = cM from scanone).
##
## Run from the qtl_mapping directory (contains results/):
##   Rscript scripts/assess-overlap.R
##
## Optional: first argument = interval column value (default lod_1.5).
##   Rscript scripts/assess-overlap.R bayes

rm(list = ls())

qtl_root <- path.expand("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping")
setwd(qtl_root)

## Default overlap input: latest results/peak_intervals/<YYYYMMDD_alpha*>/ (May2026 layout),
## or results/peak_intervals/ itself if *_peak_intervals.tsv live there (legacy).
latest_peak_intervals_dir <- function(root = file.path("results", "peak_intervals")) {
  if (!dir.exists(root)) {
    return(root)
  }
  pat <- "_peak_intervals\\.tsv$"
  legacy <- list.files(root, pattern = pat, full.names = TRUE)
  if (length(legacy) > 0) {
    return(root)
  }
  subs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  if (length(subs) == 0) {
    return(root)
  }
  has_tsv <- function(d) {
    length(list.files(d, pattern = pat)) > 0L
  }
  ok <- subs[vapply(subs, has_tsv, logical(1L))]
  if (length(ok) == 0) {
    return(root)
  }
  info <- file.info(ok)
  mt <- suppressWarnings(as.numeric(info$mtime))
  mt[is.na(mt)] <- -Inf
  ok[[which.max(mt)]]
}

trait_from_interval_filename <- function(path) {
  trimws(sub("_peak_intervals\\.tsv$", "", basename(path)))
}

find_qtl_overlaps <- function(dir,
                              pattern    = "_peak_intervals\\.tsv$",
                              interval   = "lod_1.5",
                              out_file   = NULL) {
  suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(purrr)
  })

  if (is.null(out_file)) {
    safe <- gsub("[^A-Za-z0-9._-]+", "_", interval)
    out_file <- file.path(dir, paste0("qtl_overlaps_", safe, ".tsv"))
  }

  if (!dir.exists(dir)) {
    stop("Input directory does not exist: ", dir)
  }

  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) {
    stop(
      "No files matching pattern '", pattern, "' in:\n  ", dir,
      "\nUse latest_peak_intervals_dir() or pass the run folder under results/peak_intervals/."
    )
  }

  dat <- map_dfr(
    files,
    function(fp) {
      read_tsv(fp, col_types = cols(), show_col_types = FALSE) %>%
        mutate(trait = trait_from_interval_filename(fp))
    }
  )

  need <- c("interval", "chr", "pos", "trait", "lod")
  miss <- setdiff(need, names(dat))
  if (length(miss) > 0L) {
    stop("Missing columns in peak interval TSVs: ", paste(miss, collapse = ", "))
  }

  dat_f <- dat %>% filter(.data$interval == .env$interval)
  if (nrow(dat_f) == 0L) {
    stop(
      "No rows with interval == '", interval, "'. ",
      "Available values: ", paste(unique(dat$interval), collapse = ", ")
    )
  }

  ## Peak map position (cM): pos at maximum LOD within this interval × chromosome.
  peaks <- dat_f %>%
    group_by(.data$trait, .data$chr) %>%
    arrange(desc(.data$lod), .data$pos) %>%
    slice(1L) %>%
    ungroup() %>%
    select(trait, chr, peak_cM = pos)

  bounds <- dat_f %>%
    group_by(.data$trait, .data$chr) %>%
    summarize(
      start = min(.data$pos, na.rm = TRUE),
      end   = max(.data$pos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(length = .data$end - .data$start) %>%
    left_join(peaks, by = c("trait", "chr"))

  traits <- unique(bounds$trait)
  if (length(traits) < 2L) {
    overlaps <- tibble::tibble(
      chr = character(),
      trait.1 = character(),
      trait.2 = character(),
      start.1 = numeric(),
      end.1 = numeric(),
      start.2 = numeric(),
      end.2 = numeric(),
      overlap_start = numeric(),
      overlap_end = numeric(),
      overlap_len = numeric(),
      frac_sum = numeric(),
      jaccard = numeric(),
      frac_small = numeric(),
      peak_distance_cM = numeric()
    )
    message("Only ", length(traits), " trait(s) with interval '", interval, "'; writing empty overlaps table.")
  } else {
    trait_pairs <- combn(traits, 2, simplify = FALSE)

    overlaps <- map_dfr(trait_pairs, function(pair) {
      t1 <- pair[1]
      t2 <- pair[2]
      df1 <- filter(bounds, .data$trait == t1)
      df2 <- filter(bounds, .data$trait == t2)

      inner_join(df1, df2, by = "chr", suffix = c(".1", ".2")) %>%
        mutate(
          overlap_start = pmax(.data$start.1, .data$start.2),
          overlap_end   = pmin(.data$end.1, .data$end.2),
          overlap_len   = pmax(0, .data$overlap_end - .data$overlap_start),
          length1       = .data$length.1,
          length2       = .data$length.2,
          frac_sum   = if_else(.data$overlap_len > 0,
            .data$overlap_len / (.data$length1 + .data$length2),
            0
          ),
          jaccard = if_else(
            .data$overlap_len > 0,
            .data$overlap_len / (.data$length1 + .data$length2 - .data$overlap_len),
            0
          ),
          frac_small = if_else(
            .data$overlap_len > 0,
            .data$overlap_len / pmin(.data$length1, .data$length2),
            0
          ),
          peak_distance_cM = abs(.data$peak_cM.1 - .data$peak_cM.2)
        ) %>%
        filter(.data$overlap_len > 0) %>%
        transmute(
          .data$chr,
          trait.1 = .data$trait.1,
          trait.2 = .data$trait.2,
          .data$start.1,
          .data$end.1,
          .data$start.2,
          .data$end.2,
          .data$overlap_start,
          .data$overlap_end,
          .data$overlap_len,
          .data$peak_distance_cM,
          .data$frac_sum,
          .data$jaccard,
          .data$frac_small
        )
    })
  }

  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  write_tsv(overlaps, out_file)
  message("Wrote: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))

  invisible(overlaps)
}

## Run only when this file is the Rscript --file= target (not when source()'d from -e).
is_assess_overlap_main <- function() {
  args <- commandArgs(FALSE)
  i <- grep("^--file=", args)
  if (length(i) < 1L) {
    return(FALSE)
  }
  f <- sub("^--file=", "", args[[i[1L]]])
  grepl("assess-overlap\\.R$", f, ignore.case = TRUE)
}

interval_arg <- commandArgs(trailingOnly = TRUE)
interval_use <- if (length(interval_arg) >= 1L) interval_arg[[1]] else "lod_1.5"

if (is_assess_overlap_main()) {
  in_dir <- latest_peak_intervals_dir()
  find_qtl_overlaps(
    dir      = in_dir,
    interval = interval_use,
    out_file = file.path(
      in_dir,
      paste0("qtl_overlaps_", gsub("[^A-Za-z0-9._-]+", "_", interval_use), ".tsv")
    )
  )
}
