#!/bin/bash
#
# MultiQC wrapper script for MAG workflow
# Aggregates QC reports from all tools
#

set -euo pipefail

echo "=== MultiQC Report Generation ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Parse arguments to extract output directory
OUTPUT_DIR=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--outdir)
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

# Run MultiQC
multiqc "${ARGS[@]}"

# Copy reports to expected locations
if [[ -n "$OUTPUT_DIR" ]]; then
    if [[ -f "${OUTPUT_DIR}/multiqc_report.html" ]]; then
        cp "${OUTPUT_DIR}/multiqc_report.html" "multiqc_report.html"
    fi
    if [[ -f "${OUTPUT_DIR}/multiqc_data/multiqc_data.json" ]]; then
        cp "${OUTPUT_DIR}/multiqc_data/multiqc_data.json" "multiqc_data.json"
    fi
fi

echo "MultiQC completed successfully"
