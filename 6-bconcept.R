#!/usr/bin/env Rscript

library(tidyverse)
library(Rtsne)
library(ggthemes)
library(paletteer)
library(ggrepel)


# Globals ----

CLUSTER <- "results/HDB_wdl.Rds"
DISTANCE_M <- "results/jacc_q7.Rds"
PSEARCH <- "results/jacc_q7_psearch.Rds"

HITS_CHAR <- "results/ARCH.tsv"

SEED <- 455538
set.seed(SEED)

# Helpers ----

# Bernoulli mask
# generate a mask for a vector using a Bernoulli
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
HDB_wdl <- read_rds(CLUSTER)
perplexity_search <- read_rds(PSEARCH)

# Main ----

tsne <- perplexity_search$perp_512

plot_data <- as_tibble(tsne$Y)
plot_data$neID <- rownames(wdl)
plot_data$cluster <- factor(HDB_wdl$cluster)
plot_data$outlier <- HDB_wdl$cluster == 0

plot_data <- left_join(plot_data, hits_char, join_by(neID))
plot_data_all <- plot_data


plot_data <- plot_data |> arrange(desc(LARCH))
mask <- bmask(plot_data$V1, .prob = 1.0)
famask <- bmask(plot_data$family, .prob = 1.0)

qsmall <- quantile(plot_data$LARCH, 10 / 100)
qbig <- quantile(plot_data$LARCH, 90 / 100)

plot_data <- plot_data |>
  arrange(desc(LARCH)) |>
  mutate(Ltop = 1:nrow(plot_data)) |>
  arrange(LARCH) |>
  mutate(Lbot = 1:nrow(plot_data))

max_colors <- max(as.integer(as.character(plot_data$cluster))) + 1
mycolors <- unsort(paletteer_c("grDevices::Dark 3", max_colors))

ggplot(plot_data) +
  # Inliers
  geom_jitter(
    aes(
      x = delbyNA(V1, cluster != 0),
      y = delbyNA(V2, cluster != 0)
    ),
    alpha = 1 / 4,
    size = 0.74
  ) +
  # Outliers
  geom_jitter(
    aes(
      x = delbyNA(V1, cluster == 0),
      y = delbyNA(V2, cluster == 0)
    ),
    color = "red",
    alpha = 1 / 8,
    size = 0.18
  ) +
  # Clusters
  geom_jitter(
    aes(
      x = delbyNA(V1, mask & cluster != 0),
      y = delbyNA(V2, mask & cluster != 0),
      color = delbyNA(cluster, mask & cluster != 0)
    ),
    size = 6,
    shape = 1,
    alpha = 1 / 6
  ) +
  # geom_text_repel(
  #   aes(
  #     x = delbyNA(V1, famask),
  #     y = delbyNA(V2, famask),
  #     label = delbyNA(family, famask)
  #   ),
  #   max.overlaps = 1024, size = 0.32,
  #   segment.color = NA, alpha = 1 / 3,
  #   force = 2, force_pull = 0.5
  # ) +
  # Path for small
    geom_jitter(aes(x = delbyNA(V1, Lbot <= 16), y = delbyNA(V2, Lbot <= 16)), shape = 25, size = 1.6, color = "#29C75D") +
  # Path for big
  geom_jitter(aes(x = delbyNA(V1, Ltop <= 16), y = delbyNA(V2, Ltop <= 16)), shape = 24, size = 1.6, color = "#C72992", alpha) +
  # bad quality
  geom_jitter(
    aes(
      x = delbyNA(V1, LARCH <= 7 | lengtho_ext != 25),
      y = delbyNA(V2, LARCH <= 7 | lengtho_ext != 25)
    ),
    color = "black",
    alpha = 1,
    size = 1.2,
    shape = 4) +
  # Theme
  theme_fivethirtyeight(base_size = 18) +
  theme(legend.position = "none") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) +
  scale_color_manual(values = mycolors) +
  labs(
    title = "PTTG-neighborhoods can be clustered",
    subtitle = "tSNE axes, HDBSCAN clusters, Jaccard q-gram distance.",
    caption = "author: Becerra-Soto E."
  )
ggsave("results/Bcon.pdf", width = 11, height = 8.5, units = "in", dpi = 300)

  

