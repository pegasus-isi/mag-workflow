#!/usr/bin/env python3

"""
Fetch nf-core/mag test data for running the MAG Pegasus workflow.

This script downloads test FASTQ files from the nf-core test-datasets repository
and generates a compatible samplesheet for the MAG Pegasus workflow.

Usage:
    ./fetch_test_data.py [--output-dir test_data] [--samplesheet test_samplesheet.csv]
"""

import argparse
import csv
import os
import sys
import urllib.request
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple


# nf-core/mag test data URLs
TEST_DATA_BASE_URL = "https://github.com/nf-core/test-datasets/raw/mag/test_data"

# Test samples configuration
TEST_SAMPLES = [
    {
        "sample": "test_minigut",
        "group": "minigut",
        "fastq_1": f"{TEST_DATA_BASE_URL}/test_minigut_R1.fastq.gz",
        "fastq_2": f"{TEST_DATA_BASE_URL}/test_minigut_R2.fastq.gz",
    },
    {
        "sample": "test_minigut_sample2",
        "group": "minigut",
        "fastq_1": f"{TEST_DATA_BASE_URL}/test_minigut_sample2_R1.fastq.gz",
        "fastq_2": f"{TEST_DATA_BASE_URL}/test_minigut_sample2_R2.fastq.gz",
    },
]

# Optional: Mock databases for testing (small versions)
TEST_DATABASES = {
    "busco": "https://raw.githubusercontent.com/nf-core/test-datasets/mag/databases/busco/bacteria_odb10.2024-01-08.tar.gz",
    "gtdbtk": "https://raw.githubusercontent.com/nf-core/test-datasets/mag/databases/gtdbtk/gtdbtk_mockup_20250422.tar.gz",
    "cat": "https://raw.githubusercontent.com/nf-core/test-datasets/mag/databases/cat/minigut_cat.tar.gz",
}


def download_file(url: str, output_path: str, show_progress: bool = True) -> bool:
    """
    Download a file from URL to the specified path.

    Args:
        url: URL to download from
        output_path: Local path to save the file
        show_progress: Whether to show download progress

    Returns:
        True if download successful, False otherwise
    """
    try:
        # Create parent directory if needed
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        # Skip if file already exists
        if os.path.exists(output_path):
            print(f"  [SKIP] Already exists: {os.path.basename(output_path)}")
            return True

        print(f"  [DOWNLOAD] {os.path.basename(output_path)}")

        # Download with progress
        def report_progress(block_num, block_size, total_size):
            if show_progress and total_size > 0:
                downloaded = block_num * block_size
                percent = min(100, downloaded * 100 // total_size)
                bar_length = 30
                filled = int(bar_length * percent // 100)
                bar = '=' * filled + '-' * (bar_length - filled)
                print(f"\r    [{bar}] {percent}%", end='', flush=True)

        urllib.request.urlretrieve(url, output_path, reporthook=report_progress if show_progress else None)

        if show_progress:
            print()  # New line after progress bar

        return True

    except Exception as e:
        print(f"\n  [ERROR] Failed to download {url}: {e}")
        return False


def verify_file(file_path: str) -> Tuple[bool, int]:
    """
    Verify a downloaded file exists and has content.

    Returns:
        Tuple of (is_valid, file_size)
    """
    if not os.path.exists(file_path):
        return False, 0

    size = os.path.getsize(file_path)
    return size > 0, size


def fetch_test_data(output_dir: str, include_databases: bool = False) -> List[Dict]:
    """
    Download test data files to the specified directory.

    Args:
        output_dir: Directory to save downloaded files
        include_databases: Whether to also download mock databases

    Returns:
        List of sample dictionaries with local file paths
    """
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print(f"Downloading nf-core/mag test data to: {output_dir}")
    print("-" * 60)

    downloaded_samples = []

    for sample in TEST_SAMPLES:
        sample_name = sample["sample"]
        print(f"\nSample: {sample_name}")

        # Download R1
        r1_filename = os.path.basename(sample["fastq_1"])
        r1_path = os.path.join(output_dir, r1_filename)
        r1_success = download_file(sample["fastq_1"], r1_path)

        # Download R2
        r2_filename = os.path.basename(sample["fastq_2"])
        r2_path = os.path.join(output_dir, r2_filename)
        r2_success = download_file(sample["fastq_2"], r2_path)

        if r1_success and r2_success:
            downloaded_samples.append({
                "sample": sample_name,
                "group": sample["group"],
                "fastq_1": r1_path,
                "fastq_2": r2_path,
            })
        else:
            print(f"  [WARNING] Skipping sample {sample_name} due to download failures")

    # Optionally download mock databases
    if include_databases:
        print(f"\nDownloading mock databases...")
        db_dir = os.path.join(output_dir, "databases")
        os.makedirs(db_dir, exist_ok=True)

        for db_name, db_url in TEST_DATABASES.items():
            db_filename = os.path.basename(db_url)
            db_path = os.path.join(db_dir, db_filename)
            download_file(db_url, db_path)

    return downloaded_samples


def generate_samplesheet(samples: List[Dict], samplesheet_path: str) -> None:
    """
    Generate a CSV samplesheet compatible with the MAG Pegasus workflow.

    Args:
        samples: List of sample dictionaries
        samplesheet_path: Path to write the samplesheet
    """
    fieldnames = ["sample", "fastq_1", "fastq_2", "group"]

    with open(samplesheet_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for sample in samples:
            writer.writerow({
                "sample": sample["sample"],
                "fastq_1": sample["fastq_1"],
                "fastq_2": sample["fastq_2"],
                "group": sample["group"],
            })

    print(f"\nSamplesheet written to: {samplesheet_path}")


def print_summary(samples: List[Dict], samplesheet_path: str) -> None:
    """Print a summary of downloaded data and next steps."""
    print("\n" + "=" * 60)
    print("TEST DATA DOWNLOAD COMPLETE")
    print("=" * 60)

    print(f"\nDownloaded {len(samples)} samples:")
    for sample in samples:
        r1_valid, r1_size = verify_file(sample["fastq_1"])
        r2_valid, r2_size = verify_file(sample["fastq_2"])
        status = "[OK]" if (r1_valid and r2_valid) else "[ERROR]"
        print(f"  {status} {sample['sample']}: R1={r1_size//1024}KB, R2={r2_size//1024}KB")

    print(f"\nSamplesheet: {samplesheet_path}")

    print("\n" + "-" * 60)
    print("NEXT STEPS")
    print("-" * 60)
    print(f"""
1. Generate the workflow:
   ./workflow_generator.py --samplesheet {samplesheet_path} --output workflow.yml

2. Or use the --test flag directly:
   ./workflow_generator.py --test --output workflow.yml

3. Submit to HTCondor:
   pegasus-plan --submit -s condorpool -o local workflow.yml

Note: For a quick test run, consider using these flags to skip
resource-intensive steps:
   ./workflow_generator.py --samplesheet {samplesheet_path} \\
       --skip-taxonomy --skip-annotation --output workflow.yml
""")


def main():
    parser = argparse.ArgumentParser(
        description="Fetch nf-core/mag test data for Pegasus MAG workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download test data to default location
  %(prog)s

  # Specify custom output directory
  %(prog)s --output-dir /path/to/test_data

  # Also download mock databases for testing
  %(prog)s --include-databases

  # Generate samplesheet only (if data already downloaded)
  %(prog)s --output-dir existing_data --samplesheet my_samples.csv
        """
    )

    parser.add_argument(
        "--output-dir", "-d",
        type=str,
        default="./test_data",
        help="Directory to download test data (default: ./test_data)"
    )

    parser.add_argument(
        "--samplesheet", "-s",
        type=str,
        default="test_samplesheet.csv",
        help="Output samplesheet path (default: test_samplesheet.csv)"
    )

    parser.add_argument(
        "--include-databases",
        action="store_true",
        help="Also download mock databases for GTDB-Tk, CheckM2, etc."
    )

    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Suppress progress output"
    )

    args = parser.parse_args()

    # Fetch test data
    samples = fetch_test_data(
        output_dir=args.output_dir,
        include_databases=args.include_databases
    )

    if not samples:
        print("Error: No samples were downloaded successfully")
        sys.exit(1)

    # Generate samplesheet
    samplesheet_path = os.path.abspath(args.samplesheet)
    generate_samplesheet(samples, samplesheet_path)

    # Print summary
    if not args.quiet:
        print_summary(samples, samplesheet_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
