# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 05: Gene-Set Intersections and Jaccard Similarity
#
# Purpose:
#   - Extract UP and DOWN gene sets by infection stage
#   - Compare differential expression signatures
#   - Quantify gene-set overlap using Jaccard similarity
# ============================================================

library(ggplot2)
library(UpSetR)

# 1. LOAD DIFERENTIAL EXPRESSION RESULTS
de_results <- readRDS("data/processed/differential_expression_results.rds")

# 2. BUILD DIRECTION-SPECIFIC GENE SETS
up_sets <- lapply(de_results, 
    function(df) df$gene_id[df$direction == "UP"]
)

down_sets <- lapply(de_results, 
    function(df) df$gene_id[df$direction == "DOWN"]
)

# 3. Verify that genes IDs are unique within each set
stopifnot(
    all(vapply(up_sets, function(x) !anyDuplicated(x), logical(1))),
    all(vapply(down_sets, function(x) !anyDuplicated(x), logical(1)))
)

cat("\n=== UPREGULATED GENES SET SIZES ===\n")
print(lengths(up_sets))
cat("\n=== DOWNREGULATED GENES SET SIZES ===\n")
print(lengths(down_sets))

# 3. FIRST JACCARD COMPARISION: D7 VS D9
genes_d7 <- up_sets[["D7_vs_D0"]]
genes_d9 <- up_sets[["D9_vs_D0"]]

# Genes shared between D7 and D9
shared_genes <- intersect(genes_d7, genes_d9)
# Genes present in at least one of the two sets
all_genes <- union(genes_d7, genes_d9)
# Calculate Jaccard similarity
jaccard_d7_d9 <- length(shared_genes) / length(all_genes)

cat("\n=== JACCARD SIMILARITY: D7 VS D9 (UP) ===\n")

cat("Genes UP in D7:", length(genes_d7), "\n")
cat("Genes UP in D9:", length(genes_d9), "\n")
cat("Shared genes:", length(shared_genes), "\n")
cat("Union size:", length(all_genes), "\n")
cat("Jaccard index:", round(jaccard_d7_d9, 3), "\n")


# 4. COMPUTE PAIRWISE JACCARD SIMILARITIES

# Calculate Jaccard similarity between two gene sets
jaccard_similarity <- function(set_a, set_b) {

  shared <- intersect(set_a, set_b)
  combined <- union(set_a, set_b)

  if (length(combined) == 0) {
    return(NA_real_)
  }

  return(length(shared) / length(combined))
}

# Calculate a Jaccard similarity matrix
calculate_jaccard_matrix <- function(gene_sets) {

  set_names <- names(gene_sets)

  similarity_matrix <- matrix(
    NA_real_,
    nrow = length(set_names),
    ncol = length(set_names),
    dimnames = list(set_names, set_names)
  )

  for (i in seq_along(set_names)) {
    for (j in seq_along(set_names)) {

      similarity_matrix[i, j] <- jaccard_similarity(
        gene_sets[[i]],
        gene_sets[[j]]
      )

    }
  }

  return(similarity_matrix)
}

# Calculate similarities separately for UP and DOWN genes
jaccard_up <- calculate_jaccard_matrix(up_sets)
jaccard_down <- calculate_jaccard_matrix(down_sets)

cat("\n=== JACCARD SIMILARITY MATRIX: UP ===\n")
print(round(jaccard_up, 3))

cat("\n=== JACCARD SIMILARITY MATRIX: DOWN ===\n")
print(round(jaccard_down, 3))

# Export both similarity matrices
write.csv(jaccard_up, "results/tables/jaccard_up.csv")
write.csv(jaccard_down, "results/tables/jaccard_down.csv")


# 5. VISUALIZE JACCARD SIMILARITY MATRICES

# Convert a similarity matrix into long-format data
prepare_jaccard_data <- function(similarity_matrix, direction) {

  df <- as.data.frame(as.table(similarity_matrix))

  colnames(df) <- c("stage_x", "stage_y", "jaccard")

  # Simplify labels: D7_vs_D0 -> D7
  df$stage_x <- sub("_vs_D0$", "", df$stage_x)
  df$stage_y <- sub("_vs_D0$", "", df$stage_y)

  df$direction <- direction

  return(df)
}

# Prepare UP and DOWN matrices
up_data <- prepare_jaccard_data(jaccard_up, "UP")
down_data <- prepare_jaccard_data(jaccard_down, "DOWN")

# Combine both datasets
jaccard_plot_data <- rbind(up_data, down_data)

# Define chronological order
stage_order <- c("D3", "D5", "D7", "D9", "D11", "D13")

jaccard_plot_data$stage_x <- factor(
  jaccard_plot_data$stage_x,
  levels = stage_order
)

jaccard_plot_data$stage_y <- factor(
  jaccard_plot_data$stage_y,
  levels = rev(stage_order)
)

# Generate heatmaps
p_jaccard <- ggplot(
  jaccard_plot_data,
  aes(x = stage_x, y = stage_y, fill = jaccard)
) +
  geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    aes(
      label = sprintf("%.2f", jaccard),
      color = jaccard >= 0.6
    ),
    size = 3.5
  ) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_gradient(
    low = "#EDF3F4",
    high = "#176A75",
    limits = c(0, 1),
    name = "Jaccard"
  ) +
  scale_color_manual(
    values = c("FALSE" = "#1D2B2B", "TRUE" = "white"),
    guide = "none"
  ) +
  coord_equal() +
  labs(
    title = "Temporal similarity of differential expression signatures",
    subtitle = "Pairwise Jaccard similarity of UP and DOWN gene sets",
    x = "Infection stage",
    y = "Infection stage"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# Export heatmaps
ggsave(
  filename = "results/plots/jaccard_similarity_heatmaps.png",
  plot = p_jaccard,
  width = 13,
  height = 6,
  dpi = 300
)

cat("\nJaccard similarity heatmaps exported successfully.\n")


# 6. VISUALIZE MULTI-SET INTERSECTIONS WITH UpSetR

# Create readable set labels
up_sets_plot <- up_sets
down_sets_plot <- down_sets

names(up_sets_plot) <- sub("_vs_D0$", "", names(up_sets_plot))
names(down_sets_plot) <- sub("_vs_D0$", "", names(down_sets_plot))

# Convert gene sets into binary membership tables
up_membership <- UpSetR::fromList(up_sets_plot)
down_membership <- UpSetR::fromList(down_sets_plot)

# Define the display order
stage_order <- c("D3", "D5", "D7", "D9", "D11", "D13")

# Plot UP gene intersections
png(
  "results/plots/upset_temporal_up.png",
  width = 3600,
  height = 2100,
  res = 300
)

UpSetR::upset(
  up_membership,
  sets = rev(stage_order),
  keep.order = TRUE,
  nsets = 6,
  nintersects = 25,
  order.by = "freq",
  mainbar.y.label = "Genes in exact intersection",
  sets.x.label = "UP genes per stage",
  main.bar.color = "#176A75",
  text.scale = 1.1
)

dev.off()

# Plot DOWN gene intersections
png(
  "results/plots/upset_temporal_down.png",
  width = 3600,
  height = 2100,
  res = 300
)

UpSetR::upset(
  down_membership,
  sets = rev(stage_order),
  keep.order = TRUE,
  nsets = 6,
  nintersects = 25,
  order.by = "freq",
  mainbar.y.label = "Genes in exact intersection",
  sets.x.label = "DOWN genes per stage",
  main.bar.color = "#A66D50",
  text.scale = 1.1
)

dev.off()

# Genes shared across all six infection stages
core_up <- Reduce(intersect, up_sets)
core_down <- Reduce(intersect, down_sets)

cat("\n=== GENES SHARED ACROSS ALL STAGES ===\n")
cat("Shared UP genes:", length(core_up), "\n")
cat("Shared DOWN genes:", length(core_down), "\n")

cat("\nUpSet plots exported successfully.\n")

# 7. EXPORT SHARED GENE SIGNATURES

# Verify that UP and DOWN core signatures do not overlap
stopifnot(
  length(intersect(core_up, core_down)) == 0
)

# Export genes shared across all infection stages
write.csv(
  data.frame(gene_id = core_up),
  file = "results/tables/core_up_six_stages.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = core_down),
  file = "results/tables/core_down_six_stages.csv",
  row.names = FALSE
)

cat("\nCore gene signatures exported successfully.\n")
