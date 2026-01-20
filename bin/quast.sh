#!/bin/bash
#
# QUAST wrapper script for MAG workflow
# Performs assembly quality assessment
#

set -euo pipefail

echo "=== QUAST Assembly QC ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Parse arguments to extract output directory
OUTPUT_DIR=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output-dir)
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

# Run QUAST
quast.py "${ARGS[@]}"

# Copy reports to expected output locations
if [[ -n "$OUTPUT_DIR" ]]; then
    SAMPLE_NAME=$(basename "$OUTPUT_DIR" | sed 's/_quast$//')
    if [[ -f "${OUTPUT_DIR}/report.tsv" ]]; then
        cp "${OUTPUT_DIR}/report.tsv" "${SAMPLE_NAME}_quast_report.tsv"
    fi
    if [[ -f "${OUTPUT_DIR}/report.html" ]]; then
        cp "${OUTPUT_DIR}/report.html" "${SAMPLE_NAME}_quast_report.html"
    fi
fi

echo "QUAST completed successfully"
