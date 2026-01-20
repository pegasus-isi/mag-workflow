#!/bin/bash
#
# Prodigal wrapper script for MAG workflow
# Performs gene prediction on metagenomic contigs
#

set -euo pipefail

echo "=== Prodigal Gene Prediction ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Run Prodigal
prodigal "$@"

echo "Prodigal completed successfully"
