#!/bin/bash
#
# FastQC wrapper script for MAG workflow
# Performs quality control on raw sequencing reads
#

set -euo pipefail

echo "=== FastQC Quality Control ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Force headless mode - unset DISPLAY to prevent GUI initialization
unset DISPLAY
export JAVA_TOOL_OPTIONS="-Djava.awt.headless=true"

# Run FastQC with xvfb-run if available, otherwise direct
if command -v xvfb-run &> /dev/null; then
    echo "Running FastQC with xvfb-run"
    xvfb-run --auto-servernum fastqc "$@"
else
    echo "Running FastQC in headless mode"
    fastqc "$@"
fi

echo "FastQC completed successfully"
