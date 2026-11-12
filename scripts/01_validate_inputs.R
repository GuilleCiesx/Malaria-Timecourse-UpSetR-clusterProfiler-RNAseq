
# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 01: Input Data Validation
#
# Purpose:
#   - Inspect count matrix and sample metadata
#   - Validate sample identifiers and raw counts
#   - Examine experimental group composition
#   - Compare author-provided counts with GEO raw counts
# ============================================================

# 1. LOAD INPUT DATA

# Load the count matrix and metadata provided by the authors
counts <- readRDS("data/raw/01_Bulk_Matrix.Rds")
metadata <- readRDS("data/raw/01_Bulk_SampleTable.Rds")

# Load the original count matrix downloaded from GEO
geo_counts <- read.delim(
  "data/raw/GSE279789_raw_counts.txt",
  row.names = 1,
  check.names = FALSE
)

# Harmonize control sample names in memory
colnames(geo_counts) <- sub(
  "^Naive_",
  "D0_",
  colnames(geo_counts)
)

# 2. INSPECT DATA STRUCTURE

cat("\n=== COUNT MATRIX ===\n")

print(class(counts))
print(dim(counts))
print(head(colnames(counts)))
print(head(rownames(counts)))

cat("\n=== SAMPLE METADATA ===\n")

print(class(metadata))
print(dim(metadata))
print(colnames(metadata))
print(head(metadata))

# 3. VALIDATE IDENTIFIERS

# Convert factor identifiers to character vectors
sample_ids <- as.character(metadata$SampleName)

# Check uniqueness and correspondence
stopifnot(
  !anyNA(colnames(counts)),
  !anyNA(rownames(counts)),
  !anyNA(sample_ids),
  anyDuplicated(colnames(counts)) == 0,
  anyDuplicated(rownames(counts)) == 0,
  anyDuplicated(sample_ids) == 0,
  identical(colnames(counts), rownames(metadata)),
  identical(colnames(counts), sample_ids),
  !anyNA(metadata$Group)
)

cat("\n=== SAMPLE IDENTIFIERS ===\n")

cat("Count matrix samples:", ncol(counts), "\n")
cat("Metadata samples:", nrow(metadata), "\n")
cat("Sample identifiers validated successfully.\n")

# 4. VALIDATE RAW COUNT VALUES

stopifnot(
  is.numeric(counts),
  !anyNA(counts),
  all(is.finite(counts)),
  all(counts >= 0),
  all(counts == round(counts))
)

cat("\n=== COUNT INTEGRITY ===\n")
cat("All raw counts are valid non-negative integers.\n")

# 5. EXPLORE EXPERIMENTAL GROUPS

cat("\n=== EXPERIMENTAL GROUPS ===\n")

print(table(metadata$Group))

cat("\n=== LIBRARY SIZES ===\n")

library_sizes <- colSums(counts)

print(summary(library_sizes))

cat("\n=== GENES WITH ZERO COUNTS ===\n")

zero_count_genes <- sum(rowSums(counts) == 0)

print(zero_count_genes)

# 6. INSPECT REPLICATE IMBALANCE

cat("\n=== D3 AND D9 SAMPLE DETAILS ===\n")

samples_to_check <- metadata[
  as.character(metadata$Group) %in% c("D3", "D9"),
  c("SampleName", "Group", "Run", "Sample_type")
]

print(samples_to_check)

# 7. COMPARE GEO AND AUTHORS' MATRICES

cat("\n=== MATRIX COMPARISON ===\n")

same_dimensions <- identical(
  dim(geo_counts),
  dim(counts)
)

same_genes <- identical(
  rownames(geo_counts),
  rownames(counts)
)

same_samples <- identical(
  colnames(geo_counts),
  colnames(counts)
)

cat("Same dimensions:", same_dimensions, "\n")
cat("Same gene identifiers:", same_genes, "\n")
cat("Same sample identifiers:", same_samples, "\n")

# Compare counts only when both matrices are aligned
stopifnot(
  same_dimensions,
  same_genes,
  same_samples
)

same_counts <- all(
  as.matrix(geo_counts) == counts
)

cat("Same count values:", same_counts, "\n")

stopifnot(same_counts)

# 8. FINAL VALIDATION REPORT

cat("\n========================================\n")
cat("INPUT VALIDATION COMPLETED SUCCESSFULLY\n")
cat("========================================\n")

cat("Genes:", nrow(counts), "\n")
cat("Samples:", ncol(counts), "\n")
cat("Zero-count genes:", zero_count_genes, "\n")
cat("GEO count matrix concordance: TRUE\n")
