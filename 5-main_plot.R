#!/usr/bin/env Rscript

library(tidyverse)
library(dbscan)
library(Rtsne)
library(irlba) # dependendy of Rtsne
library(ggthemes)
library(glue)
library(paletteer)
library(ggrepel)


# Globals ----

CLUSTER <- "results/HDB_wdl.Rds"
DISTANCE_M <- "results/jacc_q7.Rds"
PSEARCH <- "results/jacc_q7_psearch.Rds"

HITS_CHAR <- "results/ARCH.tsv"

SEED <- 424242
set.seed(SEED)

# HDBSCAN

MIN_PTS <- 4

# K-Meams

# Agglomerative

# Affinity

# Spectral

# Helpers ----

# bernoully mask
# generate a mask for a vector using a bernoulli
bmask <- function(x, .prob = 0.64, .length = FALSE) {
  if (.length) lx <- x else lx <- length(x)
  probs <- c(.prob, 1 - .prob)
  sample(c(TRUE, FALSE), lx, replace = TRUE, prob = probs)
}


# delete by NA
# FALSE entries in the mask become NAs
delbyNA <- function(x, keep_mask) {
  x[!keep_mask] <- NA
  x
}

# randomly shuffle
unsort <- function(x) {
  x[sample(1:length(x))]
}

# Read Data ----

hits_char <- read_tsv(HITS_CHAR)
wdl <- read_rds(DISTANCE_M)
perplexity_search <- read_rds(PSEARCH)

# HDBSCAN ----

# HDB_wdl <- hdbscan(as.dist(wdl), minPts = MIN_PTS)
HDB_wdl <- read_rds(CLUSTER)

# k-mediods ----

# k-means ----

# Agglomerative ----

# Main ----

tsne <- perplexity_search$perp_512

plot_data <- as_tibble(tsne$Y)
plot_data$neID <- rownames(wdl)
plot_data$cluster <- factor(HDB_wdl$cluster)
plot_data$outlier <- HDB_wdl$cluster == 0

plot_data <- left_join(plot_data, hits_char, join_by(neID))
plot_data_all <- plot_data


plot_data <- plot_data |> arrange(desc(LARCH))
mask <- bmask(plot_data$V1, .prob = 0.64)
famask <- bmask(plot_data$family, .prob = 0.04)

qsmall <- quantile(plot_data$LARCH, 10 / 100)
qbig <- quantile(plot_data$LARCH, 90 / 100)


max_colors <- max(as.integer(as.character(plot_data$cluster))) + 1
mycolors <- unsort(paletteer_c("grDevices::Dark 3", max_colors))

ggplot(plot_data) +
  geom_text_repel(
    aes(
      x = delbyNA(V1, famask),
      y = delbyNA(V2, famask),
      label = delbyNA(family, famask)
    ),
    max.overlaps = 24, size = 1.5,
    segment.color = NA, alpha = 2 / 3
  ) +
  # Path for all
  geom_path(aes(x = V1, y = V2), alpha = 1 / 8, size = 0.1) +
  # Path for small
  geom_path(aes(x = delbyNA(V1, LARCH <= qsmall), y = delbyNA(V2, LARCH <= qsmall)),
    size = 1, color = "#29C75D"
  ) +
  # Path for big
  geom_path(aes(x = delbyNA(V1, LARCH >= qbig), y = delbyNA(V2, LARCH >= qbig)),
    size = 0.3, alpha = 1 / 2, color = "#C72992"
  ) +
  # Inliers
  geom_jitter(
    aes(
      x = delbyNA(V1, cluster != 0),
      y = delbyNA(V2, cluster != 0)
    ),
    alpha = 1 / 3,
    size = 2
  ) +
  # Outliers
  geom_jitter(
    aes(
      x = delbyNA(V1, cluster == 0),
      y = delbyNA(V2, cluster == 0)
    ),
    color = "red",
    alpha = 1 / 4,
    size = 0.64
  ) +
  # Clusters
  geom_jitter(
    aes(
      x = delbyNA(V1, mask & cluster != 0),
      y = delbyNA(V2, mask & cluster != 0),
      color = delbyNA(cluster, mask & cluster != 0)
    ),
    size = 5,
    shape = 1,
    alpha = 2 / 3
  ) +
  # Theme
  theme_fivethirtyeight(base_size = 18) +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_color_manual(values = mycolors) +
  labs(
    title = "PTTG regions",
    subtitle = "t-SNE visualization, clusters obtained with HDBSCAN over q-gram distances.",
    caption = "author: Becerra-Soto E."
  )
ggsave("results/concept.pdf", width = 11, height = 8.5, units = "in", dpi = 300)
