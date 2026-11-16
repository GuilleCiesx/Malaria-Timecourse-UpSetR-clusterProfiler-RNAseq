
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 02: Exploratory Quality Control
#
# Purpose:
#   - Assess sample-level RNA-seq count characteristics
#   - Calculate library sizes and detected genes per sample
#   - Summarize the distribution of QC metrics
#   - Identify genes with zero counts across all samples
#   - Export QC metrics and zero-count gene identifiers
# ============================================================

library(ggplot2)

# 1. LOAD DATA

counts <- readRDS("data/raw/01_Bulk_Matrix.Rds")
metadata <- readRDS("data/raw/01_Bulk_SampleTable.Rds")

stopifnot(identical(colnames(counts), as.character(metadata$SampleName)))

# 2. CALCULATE SAMPLE-LEVEL QC METRICS

sample_qc <- data.frame(
  sample_id = colnames(counts),
  group = as.character(metadata$Group),
  library_size = colSums(counts),
  detected_genes = colSums(counts > 0)
)

# 3. INSPECT THE METRICS

cat("\n=== SAMPLE QC METRICS ===\n")
print(head(sample_qc))

cat("\n=== DESCRIPTIVE STATISTICS ===\n")
print(summary(sample_qc[, c("library_size", "detected_genes")]))

# 4. GROUP-LEVEL DESCRIPTIVE STATISTICS

cat("\n=== GROUP-LEVEL QC STATISTICS ===\n")

# Define the chronological order of experimental groups
sample_qc$group <- factor(sample_qc$group,
  levels = c("D0", "D3", "D5", "D7", "D9", "D11", "D13")
)

# Calculate the median of both QC metrics per group
group_medians <- aggregate(
  cbind(library_size, detected_genes) ~ group,
  data = sample_qc,
  FUN = median
)

print(group_medians)

# 5. IDENTIFY GENES WITH ZERO COUNTS

# Calculate total counts per gene across all samples
gene_total_counts <- rowSums(counts)

# Identify genes without counts in any sample
zero_count_gene_ids <- rownames(counts)[gene_total_counts == 0]
stopifnot(length(zero_count_gene_ids) == sum(gene_total_counts == 0))

cat("\n=== ZERO-COUNT GENES ===\n")

cat("Total zero-count genes:",
    length(zero_count_gene_ids), "\n")

# Display the first 20 gene identifiers
print(head(zero_count_gene_ids, 20))

# Export the complete list of gene identifiers
write.csv(
  data.frame(ensembl_gene_id = zero_count_gene_ids),
  file = "results/tables/zero_count_genes.csv",
  row.names = FALSE
)

# 6. EXPORT RESULTS

write.csv(sample_qc,
  file = "results/tables/sample_qc_metrics.csv",
  row.names = FALSE
)

cat("\n=== EXPORTED FILES ===\n")

cat("Sample QC metrics: results/tables/sample_qc_metrics.csv\n")
cat("Zero-count genes: results/tables/zero_count_genes.csv\n")

# 7. VISUALIZE SAMPLE-LEVEL QC METRICS

# Convert library sizes to millions of counts
sample_qc$library_size_millions <- sample_qc$library_size / 1e6

# A. Library sizes by experimental group
p_library <- ggplot(sample_qc,
  aes(x = group, y = library_size_millions)
) +
  geom_jitter(
  position = position_jitter(
    width = 0.12,
    height = 0,
    seed = 42
  ),
  size = 2.5,
  alpha = 0.8
) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 4,
    color = "black",
    fill = "red",
    stroke = 0.8
  ) +
  labs(
    title = "Library sizes across infection stages",
    subtitle = "Individual samples and group medians",
    x = "Days post-infection",
    y = "Total assigned counts (millions)"
  ) +
  theme_bw(base_size = 12)

# B. Detected genes by experimental group
p_genes <- ggplot(sample_qc,
  aes(x = group, y = detected_genes)
) +
  geom_jitter(
  position = position_jitter(
    width = 0.12,
    height = 0,
    seed = 42
  ),
  size = 2.5,
  alpha = 0.8
) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 4,
    color = "black",
    fill = "red",
    stroke = 0.8
  ) +
  labs(
    title = "Detected genes across infection stages",
    subtitle = "Individual samples and group medians",
    x = "Days post-infection",
    y = "Genes with at least one count"
  ) +
  theme_bw(base_size = 12)

# Export figures
ggsave(
  filename = "results/plots/qc_library_sizes.png",
  plot = p_library,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  filename = "results/plots/qc_detected_genes.png",
  plot = p_genes,
  width = 9,
  height = 5,
  dpi = 300
)

cat("\nQC figures exported successfully.\n")

# 8. RELATIONSHIP BETWEEN SAMPLE-LEVEL QC METRICS

# Calculate Spearman's rank correlation coefficient
spearman_rho <- cor(
  sample_qc$library_size,
  sample_qc$detected_genes,
  method = "spearman"
)

cat("\n=== QC METRIC CORRELATION ===\n")
cat("Spearman's rho:", round(spearman_rho, 3), "\n")

# Test the association between library size and detected genes
spearman_test <- cor.test(
  sample_qc$library_size,
  sample_qc$detected_genes,
  method = "spearman",
  exact = FALSE
)

print(spearman_test)

# Format the statistical results for the figure
cor_label <- sprintf(
  "Spearman's rho = %.3f | p = %.2e | n = %d",
  unname(spearman_test$estimate),
  spearman_test$p.value,
  nrow(sample_qc)
)


# Visualize the relationship between QC metrics
p_correlation <- ggplot(
  sample_qc,
  aes(
    x = library_size_millions,
    y = detected_genes,
    color = group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  labs(
    title = "Relationship between library size and detected genes",
    subtitle = cor_label,
    x = "Total assigned counts (millions)",
    y = "Detected genes",
    color = "Infection stage"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title.position = "plot",
    plot.margin = margin(12, 18, 12, 12)
  )

ggsave(
  filename = "results/plots/qc_library_vs_detected_genes.png",
  plot = p_correlation,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\nQC correlation plot exported successfully.\n")
