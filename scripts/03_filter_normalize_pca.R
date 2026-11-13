
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 03: Filtering, Normalization and PCA
#
# Purpose:
#   - Filter lowly expressed genes using edgeR
#   - Calculate TMM normalization factors
#   - Explore global expression patterns using PCA
# ============================================================

library(edgeR)
library(ggplot2)

# 1. LOAD AND VALIDATE INPUT DATA

counts <- readRDS("data/raw/01_Bulk_Matrix.Rds")
metadata <- readRDS("data/raw/01_Bulk_SampleTable.Rds")

stopifnot(
  identical(colnames(counts), as.character(metadata$SampleName))
)

# Define experimental groups in chronological order
groups <- factor(
  metadata$Group,
  levels = c("D0", "D3", "D5", "D7", "D9", "D11", "D13")
)

stopifnot(!anyNA(groups))

# 2. CREATE AN edgeR DGEList OBJECT

dge <- DGEList(counts = counts,group = groups)



# 3. FILTER LOWLY EXPRESSED GENES

# Identify genes with zero counts across all samples
zero_count_genes <- rowSums(counts) == 0

# Identify genes passing the expression filter
keep <- filterByExpr(dge, group = groups)

# Verify that no zero-count gene passes the filter
stopifnot(!any(keep[zero_count_genes]))

# Summarize filtering decisions
cat("\n=== EXPRESSION FILTERING ===\n")

cat("Genes before filtering:", nrow(dge), "\n")
cat("Zero-count genes removed:",
    sum(!keep & zero_count_genes), "\n")
cat("Low-expression genes removed:",
    sum(!keep & !zero_count_genes), "\n")
cat("Genes retained:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")
cat("Retention rate:",
    round(100 * mean(keep), 2), "%\n")

# Apply the filtering criteria
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Verify the number of retained genes
stopifnot(nrow(dge) == sum(keep))

# 4. TMM NORMALIZATION

# Calculate TMM normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# Create a summary of normalization results
normalization_summary <- data.frame(
  sample_id = colnames(dge),
  group = groups,
  library_size = dge$samples$lib.size,
  norm_factor = dge$samples$norm.factors
)

# Calculate effective library sizes
normalization_summary$effective_library_size <-
  normalization_summary$library_size *
  normalization_summary$norm_factor

# Validate normalization factors
stopifnot(
  all(is.finite(normalization_summary$norm_factor)),
  all(normalization_summary$norm_factor > 0)
)

# Display normalization results
cat("\n=== TMM NORMALIZATION ===\n")

print(normalization_summary)

cat("\n=== NORMALIZATION FACTOR SUMMARY ===\n")
print(summary(normalization_summary$norm_factor))

# 5. LOG-CPM TRANSFORMATION

# Calculate TMM-normalized log2 counts per million
log_cpm <- cpm(dge,
  normalized.lib.sizes = TRUE,
  log = TRUE,
  prior.count = 2
)

# Validate the transformed expression matrix
stopifnot(identical(dim(log_cpm), dim(dge$counts)), all(is.finite(log_cpm)))


cat("\n=== LOG-CPM TRANSFORMATION ===\n")
cat("Genes:", nrow(log_cpm), "\n")
cat("Samples:", ncol(log_cpm), "\n")

# Inspect the first five genes and four samples
print(round(log_cpm[1:5, 1:4], 2))

# 6. PRINCIPAL COMPONENT ANALYSIS (PCA)

# Transpose the expression matrix:
# rows = samples, columns = genes
pca_input <- t(log_cpm)

# Perform PCA using centered, unscaled gene expression
pca <- prcomp(pca_input, center = TRUE, scale. = FALSE)

# Calculate explained variance for each principal component
explained_variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)

# Create a data frame with sample coordinates
pca_scores <- data.frame(
  sample_id = rownames(pca$x),
  group = groups,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

# Validate sample correspondence
stopifnot(identical(pca_scores$sample_id, as.character(metadata$SampleName)))

cat("\n=== PRINCIPAL COMPONENT ANALYSIS ===\n")
cat("PC1 explained variance:", round(explained_variance[1], 2), "%\n")
cat("PC2 explained variance:", round(explained_variance[2], 2), "%\n")

print(head(pca_scores))


# 7. VISUALIZE PCA RESULTS

p_pca <- ggplot(pca_scores,
  aes(x = PC1, y = PC2, color = group)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text(
    aes(label = sample_id),
    vjust = -0.8,
    size = 2.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of murine malaria RNA-seq samples",
    subtitle = "TMM-normalized log2-CPM expression",
    x = paste0("PC1 (", round(explained_variance[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_variance[2], 1), "%)"),
    color = "Infection stage"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = "results/plots/pca_infection_stages.png",
  plot = p_pca,
  width = 10,
  height = 7,
  dpi = 300
)

# Export PCA sample coordinates
write.csv(pca_scores,
  file = "results/tables/pca_sample_scores.csv",
  row.names = FALSE
)

cat("\nPCA figure and sample scores exported successfully.\n")

# 8. SAVE PROCESSED DATA

# Save filtered counts and TMM normalization factors
saveRDS(dge, file = "data/processed/dge_filtered_tmm.rds")

# Save normalized log-CPM matrix for downstream exploration
saveRDS(log_cpm, file = "data/processed/log_cpm.rds")

cat("\nProcessed RNA-seq data saved successfully.\n")