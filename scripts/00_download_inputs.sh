
#!/usr/bin/env bash

# ============================================================
# Project: Temporal RNA-seq Patterns in Murine Malaria
# Script 00: Download Input Data
#
# Purpose:
#   - Download author-provided bulk RNA-seq data
#   - Retrieve GEO raw count matrix
#   - Verify downloaded files
#   - Prepare data/raw for downstream analysis
# ============================================================

set -euo pipefail

# 1. SET PROJECT DIRECTORY

cd "$(dirname "$0")/.."

mkdir -p data/raw

AUTHOR_BASE="https://raw.githubusercontent.com/LuChenLab/PbNK65/main/Data"
GEO_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279789/suppl"

# 2. DOWNLOAD AND VERIFY AUTHOR FILES

download_author_file() {
  local filename="$1"
  local expected_sha="$2"
  local destination="data/raw/$filename"

  if [[ -s "$destination" ]]; then
    echo "Existing file: $destination"
  else
    echo "Downloading: $filename"

    local temporary
    temporary=$(mktemp "${destination}.part.XXXXXX")

    if ! curl -fL --retry 3 \
      --connect-timeout 20 \
      -o "$temporary" \
      "$AUTHOR_BASE/$filename"; then
      rm -f "$temporary"
      echo "ERROR: Download failed: $filename" >&2
      exit 1
    fi

    mv "$temporary" "$destination"
  fi

  # Verify file integrity against the author's GitHub blob SHA
  local actual_sha
  actual_sha=$(git hash-object "$destination")

  if [[ "$actual_sha" != "$expected_sha" ]]; then
    echo "ERROR: Integrity check failed: $filename" >&2
    exit 1
  fi

  echo "Integrity verified: $filename"
}

download_author_file \
  "01_Bulk_Matrix.Rds" \
  "beafae3aaad463595444541ba340346d37333b02"

download_author_file \
  "01_Bulk_SampleTable.Rds" \
  "c2a730e433b0202cbc0f0891ff31a13403a6e7d8"

# 3. DOWNLOAD GEO RAW COUNTS

GEO_GZ="data/raw/GSE279789_raw_counts.txt.gz"
GEO_TXT="data/raw/GSE279789_raw_counts.txt"

if [[ -s "$GEO_TXT" ]]; then
  echo "Existing GEO count matrix: $GEO_TXT"
else

  if [[ ! -s "$GEO_GZ" ]]; then
    echo "Downloading GEO raw counts..."

    temporary=$(mktemp "${GEO_GZ}.part.XXXXXX")

    if ! curl -fL --retry 3 \
      --connect-timeout 20 \
      -o "$temporary" \
      "$GEO_BASE/GSE279789_raw_counts.txt.gz"; then
      rm -f "$temporary"
      echo "ERROR: GEO download failed." >&2
      exit 1
    fi

    mv "$temporary" "$GEO_GZ"
  fi

  # Check compressed-file integrity
  gzip -t "$GEO_GZ"

  echo "Decompressing GEO count matrix..."

  temporary=$(mktemp "${GEO_TXT}.part.XXXXXX")

  if ! gzip -dc "$GEO_GZ" > "$temporary"; then
    rm -f "$temporary"
    echo "ERROR: GEO decompression failed." >&2
    exit 1
  fi

  mv "$temporary" "$GEO_TXT"
fi

# 4. FINAL VALIDATION

for file in \
  "data/raw/01_Bulk_Matrix.Rds" \
  "data/raw/01_Bulk_SampleTable.Rds" \
  "data/raw/GSE279789_raw_counts.txt"
do
  if [[ ! -s "$file" ]]; then
    echo "ERROR: Missing or empty file: $file" >&2
    exit 1
  fi
done

echo ""
echo "========================================"
echo "INPUT DATA PREPARATION COMPLETED"
echo "========================================"
echo "All required input files are available."
echo "Next step: Rscript scripts/01_validate_inputs.R"
