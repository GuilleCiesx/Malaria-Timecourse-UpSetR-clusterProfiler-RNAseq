
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 07: Functional Enrichment Analysis
#
# Purpose:
#   - Load shared differential expression signatures
#   - Define an appropriate tested-gene background
#   - Validate gene identifiers
#   - Perform GO Biological Process enrichment
#   - Analyze UP and DOWN core signatures separately
#   - Export functional enrichment results
# ============================================================

library(clusterProfiler)
library(org.Mm.eg.db)

# 1. LOAD GENE SIGNATURES AND BACKGROUND

# Load genes shared across all six infection stages
core_up <- read.csv(
  "results/tables/core_up_six_stages.csv"
)$gene_id

core_down <- read.csv(
  "results/tables/core_down_six_stages.csv"
)$gene_id

# Load filtered and TMM-normalized edgeR object
dge <- readRDS(
  "data/processed/dge_filtered_tmm.rds"
)

# Define the universe of genes tested for differential expression
gene_universe <- rownames(dge$counts)

# 2. VALIDATE INPUT GENE SETS

stopifnot(
  length(core_up) > 0,
  length(core_down) > 0,
  !anyNA(core_up),
  !anyNA(core_down),
  anyDuplicated(core_up) == 0,
  anyDuplicated(core_down) == 0,
  anyDuplicated(gene_universe) == 0,
  all(core_up %in% gene_universe),
  all(core_down %in% gene_universe),
  length(intersect(core_up, core_down)) == 0
)

cat("\n=== FUNCTIONAL ENRICHMENT INPUTS ===\n")

cat("Core UP genes:", length(core_up), "\n")
cat("Core DOWN genes:", length(core_down), "\n")
cat("Tested gene universe:", length(gene_universe), "\n")

# 3. CHECK ENSEMBL ANNOTATION COVERAGE

# Retrieve Ensembl identifiers available in the annotation database
annotated_genes <- AnnotationDbi::keys(
  org.Mm.eg.db,
  keytype = "ENSEMBL"
)

# Count identifiers with available annotations
mapped_up <- sum(core_up %in% annotated_genes)
mapped_down <- sum(core_down %in% annotated_genes)
mapped_universe <- sum(gene_universe %in% annotated_genes)

cat("\n=== ENSEMBL ANNOTATION COVERAGE ===\n")

cat("Annotated core UP genes:", mapped_up, "\n")
cat("Annotated core DOWN genes:", mapped_down, "\n")
cat("Annotated background genes:", mapped_universe, "\n")

cat(
  "Core UP annotation coverage:",
  round(100 * mapped_up / length(core_up), 2),
  "%\n"
)

cat(
  "Core DOWN annotation coverage:",
  round(100 * mapped_down / length(core_down), 2),
  "%\n"
)

cat(
  "Background annotation coverage:",
  round(100 * mapped_universe / length(gene_universe), 2),
  "%\n"
)

# 4. GO BIOLOGICAL PROCESS ENRICHMENT: CORE UP

go_core_up <- enrichGO(
  gene = core_up,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Convert enrichment results into a data frame
go_up_results <- as.data.frame(go_core_up)

cat("\n=== GO ENRICHMENT: CORE UP ===\n")

cat(
  "Enriched biological processes:",
  nrow(go_up_results),
  "\n"
)

if (nrow(go_up_results) > 0) {

  print(
    head(
      go_up_results[
        , c("ID", "Description", "GeneRatio", "p.adjust")
      ],
      10
    )
  )

} else {

  cat("No significant GO Biological Process terms found.\n")

}

# 5. GO BIOLOGICAL PROCESS ENRICHMENT: CORE DOWN

go_core_down <- enrichGO(
  gene = core_down,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = gene_universe,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

# Convert enrichment results into a data frame
go_down_results <- as.data.frame(go_core_down)

cat("\n=== GO ENRICHMENT: CORE DOWN ===\n")

cat(
  "Enriched biological processes:",
  nrow(go_down_results),
  "\n"
)

if (nrow(go_down_results) > 0) {

  print(
    head(
      go_down_results[
        , c("ID", "Description", "GeneRatio", "p.adjust")
      ],
      10
    )
  )

} else {

  cat("No significant GO Biological Process terms found.\n")

}

# 6. EXPORT FUNCTIONAL ENRICHMENT RESULTS

dir.create(
  "results/tables/functional_enrichment",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  go_up_results,
  file = "results/tables/functional_enrichment/go_core_up_bp.csv",
  row.names = FALSE
)

write.csv(
  go_down_results,
  file = "results/tables/functional_enrichment/go_core_down_bp.csv",
  row.names = FALSE
)

cat("\n=== FUNCTIONAL ENRICHMENT SUMMARY ===\n")

cat("GO BP terms enriched in core UP:",
    nrow(go_up_results), "\n")

cat("GO BP terms enriched in core DOWN:",
    nrow(go_down_results), "\n")

cat("\nFunctional enrichment results exported successfully.\n")


# 7. VISUALIZE GO BIOLOGICAL PROCESS ENRICHMENT

# UP core signature
if (nrow(go_up_results) > 0) {

  p_go_up <- enrichplot::dotplot(
    go_core_up,
    showCategory = 15,
    x = "GeneRatio",
    color = "p.adjust",
    label_format = 45
  ) +
    ggplot2::labs(
      title = "GO Biological Process enrichment - Core UP",
      subtitle = "Genes upregulated across all six infection stages"
    ) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(
    filename = "results/plots/go_core_up_bp.png",
    plot = p_go_up,
    width = 12,
    height = 8,
    dpi = 300
  )
}

# DOWN core signature
if (nrow(go_down_results) > 0) {

  p_go_down <- enrichplot::dotplot(
    go_core_down,
    showCategory = 10,
    x = "GeneRatio",
    color = "p.adjust",
    label_format = 45
  ) +
    ggplot2::labs(
      title = "GO Biological Process enrichment - Core DOWN",
      subtitle = "Genes downregulated across all six infection stages"
    ) +
    ggplot2::theme_bw(base_size = 12)

  ggplot2::ggsave(
    filename = "results/plots/go_core_down_bp.png",
    plot = p_go_down,
    width = 12,
    height = 7,
    dpi = 300
  )
}

cat("\nGO enrichment dotplots exported successfully.\n")


# 8. INSPECT GENES CONTRIBUTING TO SELECTED GO TERMS

# Extract gene symbols associated with a GO term
get_go_genes <- function(enrichment_results, go_id) {

  gene_string <- enrichment_results$geneID[
    enrichment_results$ID == go_id
  ]

  stopifnot(length(gene_string) == 1)

  unique(strsplit(gene_string, "/", fixed = TRUE)[[1]])
}

# Investigate selected DOWN biological processes
cilium_genes <- get_go_genes(
  go_down_results,
  "GO:0003341"
)

microtubule_genes <- get_go_genes(
  go_down_results,
  "GO:0007018"
)

sperm_motility_genes <- get_go_genes(
  go_down_results,
  "GO:0097722"
)

# Inspect gene membership and overlap
cat("\n=== GENES ASSOCIATED WITH CILIUM MOVEMENT ===\n")
print(cilium_genes)

cat("\n=== SHARED CILIUM AND SPERM MOTILITY GENES ===\n")
print(intersect(cilium_genes, sperm_motility_genes))

cat("\n=== SHARED CILIUM AND MICROTUBULE GENES ===\n")
print(intersect(cilium_genes, microtubule_genes))

# Export genes associated with cilium movement
write.csv(
  data.frame(gene_symbol = sort(cilium_genes)),
  file = "results/tables/functional_enrichment/cilium_movement_genes.csv",
  row.names = FALSE
)

# 9. EXAMINE CILIARY GENE EXPRESSION ACROSS STAGES

# Load differential expression results
de_results <- readRDS(
  "data/processed/differential_expression_results.rds"
)

# Map ciliary gene symbols to Ensembl identifiers
cilium_mapping <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = cilium_genes,
  columns = "ENSEMBL",
  keytype = "SYMBOL"
)

# Retain valid identifiers from the tested gene universe
cilium_mapping <- unique(
  cilium_mapping[
    !is.na(cilium_mapping$ENSEMBL) &
      cilium_mapping$ENSEMBL %in% gene_universe,
    c("SYMBOL", "ENSEMBL")
  ]
)

# Collect log2FC and FDR values across all stages
cilium_expression <- do.call(
  rbind,
  lapply(names(de_results), function(contrast_name) {

    df <- de_results[[contrast_name]]

    selected <- merge(
      cilium_mapping,
      df[, c("gene_id", "logFC", "FDR")],
      by.x = "ENSEMBL",
      by.y = "gene_id"
    )

    selected$contrast <- contrast_name

    return(selected)
  })
)

# Export ciliary gene differential expression results
write.csv(
  cilium_expression,
  file = "results/tables/functional_enrichment/ciliary_gene_expression.csv",
  row.names = FALSE
)

cat("\n=== CILIARY GENE EXPRESSION SUMMARY ===\n")
print(cilium_expression)


# 10. VISUALIZE CILIARY GENE EXPRESSION PATTERNS

# Extract infection stage from contrast names
cilium_expression$stage <- sub(
  "_vs_D0$",
  "",
  cilium_expression$contrast
)

# Define chronological order
cilium_expression$stage <- factor(
  cilium_expression$stage,
  levels = c("D3", "D5", "D7", "D9", "D11", "D13")
)

# Validate one observation per gene and stage
stopifnot(
  nrow(cilium_expression) == 60,
  anyDuplicated(
    cilium_expression[, c("SYMBOL", "stage")]
  ) == 0
)

# Generate log2 fold-change heatmap
p_ciliary <- ggplot2::ggplot(
  cilium_expression,
  ggplot2::aes(
    x = stage,
    y = SYMBOL,
    fill = logFC
  )
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.7
  ) +
  ggplot2::geom_text(
    ggplot2::aes(
      label = sprintf("%.2f", logFC),
      color = logFC < -1.7
    ),
    size = 3.4
  ) +
  ggplot2::scale_fill_gradient(
    low = "#244E70",
    high = "#E5F0F4",
    limits = c(-3, 0),
    name = "log2FC"
  ) +
  ggplot2::scale_color_manual(
    values = c("FALSE" = "#202020", "TRUE" = "white"),
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Ciliary gene expression across infection stages",
    subtitle = "Core DOWN genes associated with GO:0003341",
    x = "Days post-infection",
    y = "Gene symbol"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

# Export figure
ggplot2::ggsave(
  filename = "results/plots/ciliary_gene_logfc_heatmap.png",
  plot = p_ciliary,
  width = 10,
  height = 7,
  dpi = 300
)

cat("\nCiliary gene log2FC heatmap exported successfully.\n")
