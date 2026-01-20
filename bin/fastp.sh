#!/bin/bash
#
# fastp wrapper script for MAG workflow
# Performs read trimming and quality filtering
#

set -euo pipefail

echo "=== fastp Read Trimming ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Run fastp
fastp "$@"

echo "fastp completed successfully"
