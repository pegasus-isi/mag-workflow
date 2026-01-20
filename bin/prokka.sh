#!/bin/bash
#
# Prokka wrapper script for MAG workflow
# Performs genome annotation
#

set -euo pipefail

echo "=== Prokka Genome Annotation ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Parse arguments to extract output directory and prefix
OUTPUT_DIR=""
PREFIX=""
ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --outdir)
            OUTPUT_DIR="$2"
            ARGS+=("$1" "$2")
            shift 2
            ;;
        --prefix)
            PREFIX="$2"
            ARGS+=("$1" "$2")
            shift 2
            ;;
        *)
            ARGS+=("$1")
            shift
            ;;
    esac
done

# Run Prokka
prokka "${ARGS[@]}"

# Copy outputs to expected locations
if [[ -n "$OUTPUT_DIR" && -n "$PREFIX" ]]; then
    if [[ -f "${OUTPUT_DIR}/${PREFIX}.gff" ]]; then
        cp "${OUTPUT_DIR}/${PREFIX}.gff" "${PREFIX}_prokka.gff"
    fi
    if [[ -f "${OUTPUT_DIR}/${PREFIX}.gbk" ]]; then
        cp "${OUTPUT_DIR}/${PREFIX}.gbk" "${PREFIX}_prokka.gbk"
    fi
    if [[ -f "${OUTPUT_DIR}/${PREFIX}.faa" ]]; then
        cp "${OUTPUT_DIR}/${PREFIX}.faa" "${PREFIX}_prokka.faa"
    fi
fi

echo "Prokka completed successfully"
