#!/bin/bash
#
# GTDB-Tk wrapper script for MAG workflow
# Performs taxonomic classification of genome bins
#

set -euo pipefail

echo "=== GTDB-Tk Taxonomic Classification ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Parse arguments to extract output directory
OUTPUT_DIR=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --out_dir)
            OUTPUT_DIR="$2"
            ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

# Run GTDB-Tk
gtdbtk "${ARGS[@]}"

# Copy summary to expected location
if [[ -n "$OUTPUT_DIR" ]]; then
    SAMPLE_NAME=$(basename "$OUTPUT_DIR" | sed 's/_gtdbtk$//')
    # Check for bacterial and/or archaeal summaries
    for summary in "${OUTPUT_DIR}"/classify/*.summary.tsv; do
        if [[ -f "$summary" ]]; then
            cp "$summary" "${SAMPLE_NAME}_gtdbtk.summary.tsv"
            break
        fi
    done
fi

echo "GTDB-Tk completed successfully"
