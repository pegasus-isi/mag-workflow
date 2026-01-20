#!/bin/bash
#
# MetaBAT2 wrapper script for MAG workflow
# Performs genome binning of metagenomic contigs
#

set -euo pipefail

echo "=== MetaBAT2 Genome Binning ==="
echo "Arguments: $@"
echo "Date: $(date)"

# Check if this is a depth calculation or binning run
if [[ "$1" == "jgi_summarize_bam_contig_depths" ]]; then
    shift
    jgi_summarize_bam_contig_depths "$@"
else
    metabat2 "$@"
fi

echo "MetaBAT2 completed successfully"
