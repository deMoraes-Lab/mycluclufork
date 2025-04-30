#!/usr/bin/env Rscript

# Globals ----

suppressPackageStartupMessages({
  library(tidyverse)
  library(stringdist)
  library(proxy)
})

HITS_CHAR <- "results/ARCH.tsv"
WEIGHTS <- c(s = 1, i = 1, d = 1, t = 1)

# Helpers ----

damerau <- function(x, y, ...) {
  stringdist::stringdist(x, y, method = "jaccard", q = 7, weight = ...)
}

pr_DB$set_entry(FUN = damerau, names = c("Damerau", "Damerau-Levenshtein", "dl"))
pr_DB$modify_entry(
  names = "dl", description = "The Damerau-Levenshtein distance for strings.",
  formula = "d(i,j) = min (\n
                   0 [i=j=0],\n
                   d(i-1,j) + 1 [i>0],\n
                   d(i,j-1) + 1 [j>0],\n
                   d(i-1, j-1)  + 1(ai!=bj) [i,j>0],\n
                   d(i-2, j-2), + 1(ai!=bj) [i,j>1 & ai=bj-1 & ai-1 = bj])"
)

split_domains <- function(x, pattern = "\\|") {
  map(x, \(i) if (is.na(i)) NA else str_split_1(i, pattern = pattern)) |>
    unlist()
}

# Main ----

RECALCULATE <- TRUE

hits_char <- read_tsv(HITS_CHAR)

if (RECALCULATE) {
  wdl <- dist(as.list(hits_char$char),
    as.list(hits_char$char),
    method = "dl", weight = WEIGHTS
  )

  rownames(wdl) <- (hits_char$neID)
  colnames(wdl) <- (hits_char$neID)
  write_rds(wdl, "jacc_q7.Rds")
} else {
  wdl <- read_rds("jacc_q7.Rds")
}
