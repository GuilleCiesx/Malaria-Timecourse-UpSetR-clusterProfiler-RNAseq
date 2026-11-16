
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 04: Differential Expression Analysis
#
# Purpose:
#   - Define the experimental design matrix
#   - Specify infection-stage contrasts against D0
#   - Fit negative binomial models using edgeR
#   - Identify differentially expressed genes
# ============================================================

library(edgeR)
library(limma)

# 1. LOAD PROCESSED DATA

# Load filtered counts and TMM normalization factors
dge <- readRDS("data/processed/dge_filtered_tmm.rds")

# Load sample metadata
metadata <- readRDS("data/raw/01_Bulk_SampleTable.Rds")

# Validate sample correspondence
stopifnot(identical(colnames(dge), as.character(metadata$SampleName)))

# 2. DEFINE EXPERIMENTAL GROUPS

groups <- factor(metadata$Group, 
    levels = c("D0", "D3", "D5", "D7", "D9", "D11", "D13")
)

stopifnot(!anyNA(groups))

# 3. BUILD DESIGN MATRIX

design <- model.matrix(~ 0 + groups)

# Use experimental group names as column names
colnames(design) <- levels(groups)
rownames(design) <- colnames(dge)

# Validate the design
stopifnot(
  nrow(design) == ncol(dge),
  qr(design)$rank == ncol(design)
)

cat("\n=== EXPERIMENTAL DESIGN ===\n")
print(design)
cat("\nSamples per group:\n")
print(table(groups))

# 4. DEFINE DIFFERENTIAL EXPRESSION CONTRASTS

contrasts <- makeContrasts(
  D3_vs_D0  = D3 - D0,
  D5_vs_D0  = D5 - D0,
  D7_vs_D0  = D7 - D0,
  D9_vs_D0  = D9 - D0,
  D11_vs_D0 = D11 - D0,
  D13_vs_D0 = D13 - D0,
  levels = design
)

cat("\n=== DIFFERENTIAL EXPRESSION CONTRASTS ===\n")
print(contrasts)

# 5. FIT QUASI-LIKELIHOOD NEGATIVE BINOMIAL MODEL

# Fit the model using filtered counts and TMM normalization
fit <- glmQLFit(dge, design = design)

cat("\n=== MODEL FITTING ===\n")
cat("Genes modeled:", nrow(fit$coefficients), "\n")
cat("Samples:", ncol(dge), "\n")
cat("Model coefficients:", ncol(design), "\n")
cat("Design residual degrees of freedom:",
    ncol(dge) - qr(design)$rank, "\n")

cat("\n=== MODEL RESIDUAL DF SUMMARY ===\n")
print(summary(fit$df.residual))

if (!is.null(fit$df.residual.adj)) {
  cat("\n=== ADJUSTED RESIDUAL DF SUMMARY ===\n")
  print(summary(fit$df.residual.adj))
}

#
# 6. DIFFERENTIAL EXPRESSION ACROSS ALL CONTRASTS

# Define the FDR significance threshold
fdr_threshold <- 0.05

# Create directory for differential expression results
dir.create(
  "results/tables/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

# Initialize containers for results
de_results <- list()
de_summary <- data.frame()

# Analyze each infection stage against D0
for (contrast_name in colnames(contrasts)) {

  # Perform quasi-likelihood F-test
  qlf <- glmQLFTest(
    fit,
    contrast = contrasts[, contrast_name]
  )

  # Extract results for all tested genes
  results <- topTags(qlf, n = Inf)$table

  # Preserve Ensembl gene identifiers as a column
  results <- data.frame(
    gene_id = rownames(results),
    results,
    row.names = NULL
  )

  # Classify genes by significance and direction
  results$direction <- ifelse(
    results$FDR < fdr_threshold & results$logFC > 0,
    "UP",
    ifelse(
      results$FDR < fdr_threshold & results$logFC < 0,
      "DOWN",
      "NS"
    )
  )

  # Store results in memory
  de_results[[contrast_name]] <- results

  # Count genes in each category
  summary_row <- data.frame(
    contrast = contrast_name,
    total_genes = nrow(results),
    significant = sum(results$FDR < fdr_threshold),
    upregulated = sum(results$direction == "UP"),
    downregulated = sum(results$direction == "DOWN")
  )

  # Add results to summary table
  de_summary <- rbind(de_summary, summary_row)

# Validate differential expression classification
stopifnot(
  all(
    de_summary$significant ==
      de_summary$upregulated + de_summary$downregulated
  )
)

  # Export complete differential expression table
  write.csv(
    results,
    file = file.path(
      "results/tables/differential_expression",
      paste0(contrast_name, ".csv")
    ),
    row.names = FALSE
  )

  cat("\nCompleted contrast:", contrast_name, "\n")
}

# 7. EXPORT DIFFERENTIAL EXPRESSION SUMMARY

cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
print(de_summary)

write.csv(
  de_summary,
  file = "results/tables/differential_expression_summary.csv",
  row.names = FALSE
)

# Save results for downstream gene-set analyses
saveRDS(
  de_results,
  file = "data/processed/differential_expression_results.rds"
)

cat("\nAll differential expression results exported.\n")

# Save fitted model and contrasts for downstream statistical tests
saveRDS(
  fit,
  file = "data/processed/edgeR_ql_fit.rds"
)

saveRDS(
  contrasts,
  file = "data/processed/edgeR_contrasts.rds"
)