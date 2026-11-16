
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 06: Statistical Robustness and Sensitivity Analysis
#
# Purpose:
#   - Evaluate sensitivity to DEG selection thresholds
#   - Compare UP and DOWN gene-set sizes
#   - Assess shared signatures across all infection stages
# ============================================================

library(ggplot2)
library(edgeR)

# 1. LOAD DIFFERENTIAL EXPRESSION RESULTS

de_results <- readRDS(
  "data/processed/differential_expression_results.rds"
)

# 2. DEFINE SELECTION SCENARIOS

scenarios <- data.frame(
  scenario = c("baseline", "strict_fdr", "large_effect"),
  fdr = c(0.05, 0.01, 0.05),
  min_abs_logFC = c(0, 0, 1)
)

# 3. INITIALIZE RESULT CONTAINERS

robustness_summary <- list()
core_summary <- list()

# 4. COMPARE SELECTION SCENARIOS

for (i in seq_len(nrow(scenarios))) {

  current <- scenarios[i, ]

  up_sets <- list()
  down_sets <- list()

  for (contrast_name in names(de_results)) {

    results <- de_results[[contrast_name]]

    # Select genes meeting the current thresholds
    selected <- (
      !is.na(results$FDR) &
      !is.na(results$logFC) &
      results$FDR < current$fdr &
      abs(results$logFC) >= current$min_abs_logFC
    )

    # Separate UP and DOWN genes
    up_sets[[contrast_name]] <- results$gene_id[
      selected & results$logFC > 0
    ]

    down_sets[[contrast_name]] <- results$gene_id[
      selected & results$logFC < 0
    ]

    # Record DEG counts per contrast
    robustness_summary[[length(robustness_summary) + 1]] <-
      data.frame(
        scenario = current$scenario,
        contrast = contrast_name,
        up = length(up_sets[[contrast_name]]),
        down = length(down_sets[[contrast_name]])
      )
  }

  # Identify genes shared across all six stages
  core_up <- Reduce(intersect, up_sets)
  core_down <- Reduce(intersect, down_sets)

  core_summary[[i]] <- data.frame(
    scenario = current$scenario,
    core_up = length(core_up),
    core_down = length(core_down)
  )
}

# 5. COMBINE AND EXPORT RESULTS

robustness_summary <- do.call(
  rbind, robustness_summary
)

core_summary <- do.call(
  rbind, core_summary
)

cat("\n=== DEG THRESHOLD SENSITIVITY ===\n")
print(robustness_summary)

cat("\n=== CORE SIGNATURE SENSITIVITY ===\n")
print(core_summary)

write.csv(
  robustness_summary,
  "results/tables/robustness_threshold_counts.csv",
  row.names = FALSE
)

write.csv(
  core_summary,
  "results/tables/robustness_core_counts.csv",
  row.names = FALSE
)


# 6. JACCARD SIMILARITY UNDER DIFFERENT THRESHOLDS

# Calculate Jaccard similarity between two gene sets
jaccard_index <- function(set_a, set_b) {
  n_union <- length(union(set_a, set_b))

  if (n_union == 0) {
    return(NA_real_)
  }

  length(intersect(set_a, set_b)) / n_union
}

# Create directory for sensitivity matrices
dir.create(
  "results/tables/jaccard_sensitivity",
  recursive = TRUE,
  showWarnings = FALSE
)

jaccard_summary <- list()

# Recalculate gene sets for every scenario
for (i in seq_len(nrow(scenarios))) {

  current <- scenarios[i, ]

  for (direction in c("UP", "DOWN")) {

    # Select genes according to current criteria
    gene_sets <- lapply(de_results, function(df) {

      selected <- !is.na(df$FDR) &
        !is.na(df$logFC) &
        df$FDR < current$fdr &
        abs(df$logFC) >= current$min_abs_logFC

      if (direction == "UP") {
        return(df$gene_id[selected & df$logFC > 0])
      } else {
        return(df$gene_id[selected & df$logFC < 0])
      }
    })

    # Initialize similarity matrix
    set_names <- names(gene_sets)
    n <- length(gene_sets)

    jaccard_matrix <- matrix(
      NA_real_,
      nrow = n,
      ncol = n,
      dimnames = list(set_names, set_names)
    )

    # Compare every pair of gene sets
    for (a in seq_len(n)) {
      for (b in seq_len(n)) {

        jaccard_matrix[a, b] <- jaccard_index(
          gene_sets[[a]],
          gene_sets[[b]]
        )
      }
    }

    # Export matrix for this scenario and direction
    write.csv(
      jaccard_matrix,
      file = paste0(
        "results/tables/jaccard_sensitivity/jaccard_",
        current$scenario, "_",
        tolower(direction), ".csv"
      )
    )

    # Record two biologically relevant comparisons
    for (pair in list(
      c("D7_vs_D0", "D9_vs_D0"),
      c("D11_vs_D0", "D13_vs_D0")
    )) {

      jaccard_summary[[length(jaccard_summary) + 1]] <-
        data.frame(
          scenario = current$scenario,
          direction = direction,
          comparison = paste(pair, collapse = " / "),
          jaccard = jaccard_matrix[pair[1], pair[2]]
        )
    }
  }
}

# Combine and export comparison results
jaccard_summary <- do.call(rbind, jaccard_summary)

cat("\n=== JACCARD ROBUSTNESS SUMMARY ===\n")
print(jaccard_summary)

write.csv(
  jaccard_summary,
  "results/tables/jaccard_robustness_summary.csv",
  row.names = FALSE
)


# 7. VISUALIZE JACCARD ROBUSTNESS

# Define readable labels for selection scenarios
jaccard_summary$scenario <- factor(
  jaccard_summary$scenario,
  levels = c("baseline", "strict_fdr", "large_effect"),
  labels = c("Baseline", "FDR < 0.01", "|log2FC| >= 1")
)

# Define readable labels for temporal comparisons
jaccard_summary$comparison <- factor(
  jaccard_summary$comparison,
  levels = c(
    "D7_vs_D0 / D9_vs_D0",
    "D11_vs_D0 / D13_vs_D0"
  ),
  labels = c("D7-D9", "D11-D13")
)

# Compare Jaccard similarities across scenarios
p_robustness <- ggplot(
  jaccard_summary,
  aes(x = scenario, y = jaccard, color = comparison)
) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3.5
  ) +
  geom_text(
    aes(label = sprintf("%.3f", jaccard)),
    position = position_dodge(width = 0.4),
    vjust = -0.9,
    size = 3,
    show.legend = FALSE
  ) +
  facet_wrap(~ direction, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Robustness of temporal gene-set similarity",
    subtitle = "Pairwise Jaccard similarity under alternative DEG thresholds",
    x = "Gene selection criterion",
    y = "Jaccard similarity",
    color = "Stage comparison"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  filename = "results/plots/jaccard_robustness.png",
  plot = p_robustness,
  width = 11,
  height = 6,
  dpi = 300
)

cat("\nJaccard robustness figure exported successfully.\n")


# 8. FORMAL FOLD-CHANGE THRESHOLD TESTING

# Load the fitted quasi-likelihood model and contrasts
fit <- readRDS("data/processed/edgeR_ql_fit.rds")
contrasts <- readRDS("data/processed/edgeR_contrasts.rds")

# Define formal testing thresholds
lfc_threshold <- 1
fdr_threshold <- 0.05

# Initialize result containers
treat_up_sets <- list()
treat_down_sets <- list()
treat_summary <- list()

# Test each infection stage against D0
for (contrast_name in colnames(contrasts)) {

  # Test whether the true absolute log2FC exceeds 1
  treat_test <- glmTreat(
    fit,
    contrast = contrasts[, contrast_name],
    lfc = lfc_threshold,
    null = "interval"
  )

  # Extract all results with adjusted p-values
  treat_results <- topTags(treat_test, n = Inf)$table

  # Identify statistically significant genes
  significant <- !is.na(treat_results$FDR) &
    treat_results$FDR < fdr_threshold

  # Separate significant UP and DOWN genes
  treat_up_sets[[contrast_name]] <- rownames(treat_results)[
    significant & treat_results$logFC > 0
  ]

  treat_down_sets[[contrast_name]] <- rownames(treat_results)[
    significant & treat_results$logFC < 0
  ]

  # Summarize formal testing results
  treat_summary[[contrast_name]] <- data.frame(
    contrast = contrast_name,
    up = length(treat_up_sets[[contrast_name]]),
    down = length(treat_down_sets[[contrast_name]])
  )
}

# Combine contrast summaries
treat_summary <- do.call(rbind, treat_summary)
rownames(treat_summary) <- NULL

cat("\n=== FORMAL glmTreat RESULTS ===\n")
print(treat_summary)

# Identify genes significant in all six stages
treat_core_up <- Reduce(intersect, treat_up_sets)
treat_core_down <- Reduce(intersect, treat_down_sets)

cat("\n=== FORMAL CORE SIGNATURES ===\n")
cat("Core UP:", length(treat_core_up), "\n")
cat("Core DOWN:", length(treat_core_down), "\n")

# Export results
write.csv(
  treat_summary,
  "results/tables/glmTreat_summary.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = treat_core_up),
  "results/tables/glmTreat_core_up.csv",
  row.names = FALSE
)

write.csv(
  data.frame(gene_id = treat_core_down),
  "results/tables/glmTreat_core_down.csv",
  row.names = FALSE
)
