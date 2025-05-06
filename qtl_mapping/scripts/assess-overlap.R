## clear workspace
rm(list = ls())

## Load libraries
library(qtl)

## Set working directory
setwd("~/Dropbox/Costus/costus-genetic-mapping/qtl_mapping/")

find_qtl_overlaps <- function(dir        = ".",
                              pattern    = "*_peak_intervals.tsv",
                              interval   = "lod_1.5",
                              out_file   = "qtl_overlaps.tsv") {
  #— libraries
  library(readr)
  library(dplyr)
  library(purrr)
  
  # 1. list + read
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  dat   <- map_dfr(files,
                   ~ read_tsv(.x, col_types = cols()) %>%
                     mutate(trait = sub(pattern, "", basename(.x)))
  )
  
  # 2. filter by chosen interval
  dat_f <- dat %>% filter(interval == !!interval)
  
  # 3. collapse to one [start,end] per trait × chr
  bounds <- dat_f %>%
    group_by(trait, chr) %>%
    summarize(
      start   = min(pos, na.rm = TRUE),
      end     = max(pos, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(length = end - start)
  
  # 4. build all trait–trait pairs
  trait_pairs <- combn(unique(bounds$trait), 2, simplify = FALSE)
  
  # 5. compute overlaps + metrics
  overlaps <- map_dfr(trait_pairs, function(pair) {
    t1 <- pair[1]; t2 <- pair[2]
    df1 <- filter(bounds, trait == t1)
    df2 <- filter(bounds, trait == t2)
    
    inner_join(df1, df2, by = "chr", suffix = c(".1", ".2")) %>%
      mutate(
        overlap_start = pmax(start.1, start.2),
        overlap_end   = pmin(end.1,   end.2),
        overlap_len   = pmax(0, overlap_end - overlap_start),
        length1       = length.1,
        length2       = length.2,
        # three overlap metrics:
        frac_sum   = if_else(overlap_len > 0,
                             overlap_len / (length1 + length2),
                             0),
        jaccard    = if_else(overlap_len > 0,
                             overlap_len / (length1 + length2 - overlap_len),
                             0),
        frac_small = if_else(overlap_len > 0,
                             overlap_len / pmin(length1, length2),
                             0)
      ) %>%
      filter(overlap_len > 0) %>%
      transmute(
        chr,
        trait.1,
        trait.2,
        start.1,
        end.1,
        start.2,
        end.2,
        overlap_start,
        overlap_end,
        overlap_len,
        frac_sum,
        jaccard,
        frac_small
      )
  })
  
  # 6. write to TSV
  write_tsv(overlaps, out_file)
  
  invisible(NULL)
}


find_qtl_overlaps(
  dir      = "./results/data_frames/",
  interval = "bayes",
  out_file = "./results/data_frames/bayes_qtl_overlaps.tsv"
)

