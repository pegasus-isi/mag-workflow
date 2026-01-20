#!/bin/bash
#
# Example usage script for MAG workflow
#

set -euo pipefail

echo "=== MAG Workflow Example Usage ==="

# Check if samplesheet is provided
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <samplesheet.csv> [assembler]"
    echo ""
    echo "Arguments:"
    echo "  samplesheet.csv  - CSV file with sample information"
    echo "  assembler        - 'megahit' (default) or 'spades'"
    echo ""
    echo "Example:"
    echo "  $0 my_samples.csv megahit"
    exit 1
fi

SAMPLESHEET="$1"
ASSEMBLER="${2:-megahit}"

# Validate samplesheet exists
if [[ ! -f "$SAMPLESHEET" ]]; then
    echo "Error: Samplesheet not found: $SAMPLESHEET"
    exit 1
fi

echo "Samplesheet: $SAMPLESHEET"
echo "Assembler: $ASSEMBLER"

# Generate workflow
echo ""
echo "Step 1: Generating Pegasus workflow..."
./workflow_generator.py \
    --samplesheet "$SAMPLESHEET" \
    --assembler "$ASSEMBLER" \
    --output-dir ./output \
    --output workflow.yml

# Show generated files
echo ""
echo "Step 2: Generated files:"
ls -la output/*.yml workflow.yml 2>/dev/null || true

# Instructions for submission
echo ""
echo "Step 3: To submit the workflow, run:"
echo "  pegasus-plan --submit -s condorpool -o local workflow.yml"
echo ""
echo "Step 4: To monitor the workflow:"
echo "  pegasus-status <run_directory>"
echo "  pegasus-analyzer <run_directory>"
echo ""
echo "=== Done ==="
