#!/bin/bash
#
# CheckM2 wrapper script for MAG workflow
# Performs quality assessment of genome bins
#

set -euo pipefail

echo "=== CheckM2 Bin Quality Assessment ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Parse arguments to extract output directory
OUTPUT_DIR=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --output-directory)
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

# Run CheckM2
checkm2 "${ARGS[@]}"

# Copy quality report to expected location
if [[ -n "$OUTPUT_DIR" && -f "${OUTPUT_DIR}/quality_report.tsv" ]]; then
    SAMPLE_NAME=$(basename "$OUTPUT_DIR" | sed 's/_checkm2$//')
    cp "${OUTPUT_DIR}/quality_report.tsv" "${SAMPLE_NAME}_checkm2_quality.tsv"
fi

echo "CheckM2 completed successfully"
