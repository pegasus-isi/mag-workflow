#!/usr/bin/env python3

"""
MAG (Metagenome-Assembled Genomes) Workflow Generator for Pegasus WMS

This script generates a Pegasus workflow for metagenomic assembly, binning,
and annotation - equivalent to the nf-core/mag Nextflow pipeline.

Pipeline steps:
1. Quality Control (FastQC, fastp)
2. Assembly (MEGAHIT or SPAdes)
3. Assembly QC (QUAST)
4. Gene Prediction (Prodigal)
5. Binning (MetaBAT2)
6. Bin Quality Assessment (CheckM2)
7. Taxonomic Classification (GTDB-Tk)
8. Genome Annotation (Prokka)
9. Report Generation (MultiQC)

Usage:
    ./workflow_generator.py --samplesheet samples.csv --output workflow.yml

    # A Slurm cluster:
    ./workflow_generator.py --samplesheet samples.csv --site-style slurm \\
        --queue cpu --project my_lab

Sites: the workflow names no scheduler. Jobs state cores, memory and a
wall-clock runtime, and plan against a site named "compute". A hosted site
catalog (--hosted-site-catalog unity.yml, ...) defines it; otherwise
custom_sites.py writes it to sites.yml as an HTCondor pool, or for Slurm via
--site-style slurm --queue --project.
"""

import argparse
import csv
import json
import os
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Pegasus imports
try:
    from Pegasus.api import (
        Arch, Job, File, OS, Properties,
        Transformation, TransformationCatalog, Container,
        ReplicaCatalog, Workflow
    )
except ImportError:
    print("Error: Pegasus Python API not found.")
    print("Install with: pip install pegasus-wms.api")
    sys.exit(1)

# Site-catalog handling. Keeps this generator free of scheduler details.
sys.path.insert(0, str(Path(__file__).parent.resolve()))
from custom_sites import (  # noqa: E402
    HOSTED_SITE, STYLES, ensure_sites_yml, hosted_catalog, is_batch_site, parse_profile,
    parse_tag_profile, worker_package_url,
)

# The container's OS, for the Pegasus worker package staged into it. Must
# run on the image's base (Apptainer/MAG_Container.def:
# mambaorg/micromamba:1.5-jammy, Ubuntu 22.04, glibc 2.35), not the submit
# host. Pegasus publishes no ubuntu_22 package, so use the oldest-glibc one,
# rhel_8 (glibc 2.28).
CONTAINER_PLATFORM = "x86_64_rhel_8"


# Default container image
# Local Apptainer image, relative to this file's directory. Build it with
# `apptainer build` (see Apptainer/MAG_Container.def); Pegasus stages the .sif
# like any other input file.
DEFAULT_CONTAINER = "Apptainer/MAG_Container.sif"

# Test data configuration
TEST_DATA_BASE_URL = "https://github.com/nf-core/test-datasets/raw/mag/test_data"
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

# Tool configurations — the only resource statements the workflow makes,
# with memory as "<n> GB" (the API converts that to MB; "16GB" is passed
# through and the planner rejects it as non-numeric),
# portable on any site. These are production-scale needs (GTDB-Tk and SPAdes
# need more than a 16 GB node); on small worker nodes cap them with
# --max-memory-gb / --max-cores. runtime is a wall-clock budget in seconds;
# batch sites (Slurm) refuse or kill jobs without one, so keep it generous.
TOOL_CONFIGS = {
    "fastqc": {"memory": "2 GB", "cores": 2, "runtime": 1800},
    "fastp": {"memory": "4 GB", "cores": 4, "runtime": 3600},
    "megahit": {"memory": "16 GB", "cores": 8, "runtime": 4 * 3600},
    "spades": {"memory": "32 GB", "cores": 16, "runtime": 8 * 3600},
    "quast": {"memory": "4 GB", "cores": 4, "runtime": 1800},
    "prodigal": {"memory": "4 GB", "cores": 1, "runtime": 3600},
    "metabat2": {"memory": "8 GB", "cores": 4, "runtime": 3600},
    "checkm2": {"memory": "16 GB", "cores": 8, "runtime": 2 * 3600},
    "gtdbtk": {"memory": "64 GB", "cores": 8, "runtime": 8 * 3600},
    "prokka": {"memory": "8 GB", "cores": 4, "runtime": 2 * 3600},
    "multiqc": {"memory": "4 GB", "cores": 2, "runtime": 1800},
}


def parse_samplesheet(samplesheet_path: str) -> List[Dict]:
    """
    Parse input samplesheet CSV.

    Expected format:
    sample,fastq_1,fastq_2,group
    sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,group1
    """
    samples = []
    with open(samplesheet_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = {
                'id': row.get('sample', row.get('id', '')),
                'fastq_1': row.get('fastq_1', row.get('R1', '')),
                'fastq_2': row.get('fastq_2', row.get('R2', '')),
                'group': row.get('group', 'default'),
                'single_end': row.get('single_end', 'false').lower() == 'true'
            }
            if sample['id'] and sample['fastq_1']:
                samples.append(sample)
    return samples


def download_test_data(output_dir: str) -> Tuple[List[Dict], str]:
    """
    Download nf-core/mag test data and generate a samplesheet.

    Args:
        output_dir: Directory to download test data to

    Returns:
        Tuple of (samples list, samplesheet path)
    """
    import urllib.request

    test_data_dir = os.path.join(output_dir, "test_data")
    os.makedirs(test_data_dir, exist_ok=True)

    print(f"Downloading nf-core/mag test data to: {test_data_dir}")
    print("-" * 60)

    downloaded_samples = []

    for sample in TEST_SAMPLES:
        sample_name = sample["sample"]
        print(f"\nSample: {sample_name}")

        # Download R1
        r1_filename = os.path.basename(sample["fastq_1"])
        r1_path = os.path.join(test_data_dir, r1_filename)
        if not os.path.exists(r1_path):
            print(f"  Downloading {r1_filename}...")
            try:
                urllib.request.urlretrieve(sample["fastq_1"], r1_path)
            except Exception as e:
                print(f"  Error downloading R1: {e}")
                continue
        else:
            print(f"  [SKIP] {r1_filename} already exists")

        # Download R2
        r2_filename = os.path.basename(sample["fastq_2"])
        r2_path = os.path.join(test_data_dir, r2_filename)
        if not os.path.exists(r2_path):
            print(f"  Downloading {r2_filename}...")
            try:
                urllib.request.urlretrieve(sample["fastq_2"], r2_path)
            except Exception as e:
                print(f"  Error downloading R2: {e}")
                continue
        else:
            print(f"  [SKIP] {r2_filename} already exists")

        downloaded_samples.append({
            'id': sample_name,
            'fastq_1': r1_path,
            'fastq_2': r2_path,
            'group': sample["group"],
            'single_end': False
        })

    # Generate samplesheet
    samplesheet_path = os.path.join(output_dir, "test_samplesheet.csv")
    with open(samplesheet_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "fastq_1", "fastq_2", "group"])
        writer.writeheader()
        for sample in downloaded_samples:
            writer.writerow({
                "sample": sample['id'],
                "fastq_1": sample['fastq_1'],
                "fastq_2": sample['fastq_2'],
                "group": sample['group'],
            })

    print(f"\nTest samplesheet written to: {samplesheet_path}")
    print("-" * 60)

    return downloaded_samples, samplesheet_path


def create_properties(sites_yml: str = "sites.yml",
                      bypass_input_staging: bool = False,
                      worker_package: Optional[str] = None,
                      hosted_site_catalog: Optional[str] = None) -> Properties:
    """Pegasus properties; which ones depend on the site in use."""
    props = Properties()
    props["pegasus.transfer.threads"] = "16"
    # Symlink rather than copy when an input already sits on the execution
    # site. A no-op otherwise, so always on.
    props["pegasus.transfer.links"] = "true"
    if hosted_site_catalog:
        # A centrally hosted site catalog (pegasushub pegasus-site-catalogs)
        # defines the execution site; pegasus-plan downloads and caches it,
        # and merges sites.yml over it.
        props["pegasus.catalog.site.repo.file"] = hosted_site_catalog
    if bypass_input_staging:
        # Jobs read inputs (notably the .sif) straight from the submit host's
        # paths. Only valid where workers share a filesystem with the submit
        # host (a Slurm cluster) — never on a condor pool staging over
        # HTCondor file transfer.
        props["pegasus.transfer.bypass.input.staging"] = "true"
    if worker_package:
        # Stage the container-compatible worker package named in the TC
        # (pegasus::worker) and never download one from inside a job: the
        # image may lack curl/wget, workers may lack internet, and the submit
        # host's kickstart may need a newer glibc than the image.
        props["pegasus.transfer.worker.package"] = "true"
        props["pegasus.transfer.worker.package.strict"] = "false"
        props["pegasus.transfer.worker.package.autodownload"] = "false"
    if os.path.isfile(sites_yml):
        # Lets pegasus-plan find sites.yml from any directory.
        props["pegasus.catalog.site"] = "YAML"
        props["pegasus.catalog.site.file"] = os.path.abspath(sites_yml)
    return props


def create_transformation_catalog(
    container_image: str,
    bind_workflow_dir: bool = False,
    worker_package: Optional[str] = None,
) -> Tuple[TransformationCatalog, Container]:
    """Create Pegasus transformation catalog with container.

    Nothing here names an execution site: the tool scripts live inside the
    container and the .sif on "local" (the submit host), shipped to wherever
    the job runs.
    """
    tc = TransformationCatalog()

    # A path ending in .sif (the default) is a locally built Apptainer image.
    # Pegasus stages the file like any other input, so image_site is the site
    # where the .sif physically lives (the submit host = "local"). A full URL
    # (docker://, https://) is passed through unchanged, and a bare name still
    # means Docker Hub. The .sif suffix is the discriminator on purpose — a bare
    # registry reference like "kthare10/mag-workflow:latest" also contains a
    # slash, so testing for a path separator would misread it as a local file.
    if "://" in container_image:
        image_url = container_image
        image_site = {"http": "web", "https": "web", "file": "local"}.get(
            container_image.split("://", 1)[0], "docker_hub")
    elif not container_image.endswith(".sif"):
        image_url = "docker://" + container_image
        image_site = "docker_hub"
    else:
        sif_path = container_image if os.path.isabs(container_image) else \
            os.path.join(os.path.dirname(os.path.abspath(__file__)), container_image)
        if not os.path.exists(sif_path):
            print(f"Warning: Apptainer image not found at {sif_path} — build it "
                  f"first with: cd mag-workflow && apptainer build {sif_path} "
                  f"Apptainer/MAG_Container.def")
        image_url = "file://" + sif_path
        image_site = "local"

    # Create container
    container = Container(
        "mag_container",
        Container.SINGULARITY,
        image=image_url,
        image_site=image_site,
    )
    if bind_workflow_dir:
        # Batch sites stage inputs as symlinks into the workflow directory,
        # and PegasusLite starts containers with --no-home, so without this
        # bind every job fails with kickstart "Unable to execute the
        # specified binary" (exit 127). Never on a condor pool: the directory
        # does not exist on its workers.
        container.add_pegasus_profile(
            container_arguments="--bind "
            + os.path.dirname(os.path.abspath(__file__)))
    tc.add_containers(container)

    # Define transformations for each tool
    tools = [
        "fastqc", "fastp", "megahit", "spades", "quast",
        "prodigal", "metabat2", "checkm2", "gtdbtk", "prokka", "multiqc"
    ]

    for tool in tools:
        config = TOOL_CONFIGS.get(
            tool, {"memory": "4 GB", "cores": 2, "runtime": 3600})
        tx = Transformation(
            tool,
            site="local",
            pfn=f"/usr/local/bin/{tool}.sh",
            is_stageable=False,  # Scripts are inside container, don't stage from submit host
            container=container
        )
        tx.add_pegasus_profile(memory=config["memory"], cores=config["cores"],
                               runtime=config["runtime"])
        # MultiQC (click) needs a UTF-8 locale; this used to sit on the
        # execution site, which no longer carries tool settings.
        tx.add_env(LANG="en_US.UTF-8")
        tc.add_transformations(tx)

    if worker_package:
        tc.add_transformations(
            Transformation(
                "worker",
                namespace="pegasus",
                site="local",
                pfn=worker_package,
                is_stageable=True,
                arch=Arch.X86_64,
                os_type=OS.LINUX,
            )
        )

    return tc, container


def create_replica_catalog(samples: List[Dict]) -> ReplicaCatalog:
    """Create replica catalog with input files."""
    rc = ReplicaCatalog()

    for sample in samples:
        # Add forward reads
        if sample['fastq_1'] and os.path.exists(sample['fastq_1']):
            rc.add_replica(
                "local",
                f"{sample['id']}_R1.fastq.gz",
                f"file://{os.path.abspath(sample['fastq_1'])}"
            )

        # Add reverse reads (if paired-end)
        if not sample['single_end'] and sample['fastq_2'] and os.path.exists(sample['fastq_2']):
            rc.add_replica(
                "local",
                f"{sample['id']}_R2.fastq.gz",
                f"file://{os.path.abspath(sample['fastq_2'])}"
            )

    return rc


def create_workflow(
    samples: List[Dict],
    assembler: str = "megahit",
    skip_binning: bool = False,
    skip_taxonomy: bool = False,
    skip_annotation: bool = False,
    skip_fastqc: bool = False,
    gtdbtk_db: Optional[str] = None,
    checkm2_db: Optional[str] = None
) -> Workflow:
    """
    Create the MAG Pegasus workflow.

    Workflow DAG:

    Input FASTQ files
         |
         v
    [FastQC] -----> QC Reports (optional)
         |
         v
    [fastp] -----> Trimmed reads
         |
         v
    [MEGAHIT/SPAdes] -----> Contigs
         |
         v
    [QUAST] -----> Assembly QC
         |
         v
    [Prodigal] -----> Gene predictions
         |
         v
    [MetaBAT2] -----> Genome bins
         |
         v
    [CheckM2] -----> Bin quality
         |
         v
    [GTDB-Tk] -----> Taxonomy
         |
         v
    [Prokka] -----> Annotations
         |
         v
    [MultiQC] -----> Final report
    """
    wf = Workflow("mag-workflow")

    # Track all QC files for MultiQC
    all_qc_files = []
    all_bin_dirs = []

    for sample in samples:
        sample_id = sample['id']
        is_paired = not sample['single_end']

        # Input files
        r1_input = File(f"{sample_id}_R1.fastq.gz")
        r2_input = File(f"{sample_id}_R2.fastq.gz") if is_paired else None

        # ============================================================
        # Step 1: Quality Control - FastQC (raw reads)
        # ============================================================
        if not skip_fastqc:
            fastqc_r1_html = File(f"{sample_id}_R1_fastqc.html")
            fastqc_r1_zip = File(f"{sample_id}_R1_fastqc.zip")

            fastqc_job = Job("fastqc")
            fastqc_job.add_args("--outdir", ".", "--threads", str(TOOL_CONFIGS["fastqc"]["cores"]))
            fastqc_job.add_inputs(r1_input)
            # The input FASTQ must also be on the command line: FastQC
            # with no file arguments starts its GUI and dies headless.
            fastqc_job.add_args(r1_input)
            fastqc_job.add_outputs(fastqc_r1_html, fastqc_r1_zip, stage_out=True)

            if is_paired:
                fastqc_r2_html = File(f"{sample_id}_R2_fastqc.html")
                fastqc_r2_zip = File(f"{sample_id}_R2_fastqc.zip")
                fastqc_job.add_inputs(r2_input)
                fastqc_job.add_args(r2_input)
                fastqc_job.add_outputs(fastqc_r2_html, fastqc_r2_zip, stage_out=True)
                all_qc_files.extend([fastqc_r1_zip, fastqc_r2_zip])
            else:
                all_qc_files.append(fastqc_r1_zip)

            wf.add_jobs(fastqc_job)

        # ============================================================
        # Step 2: Read Trimming - fastp
        # ============================================================
        trimmed_r1 = File(f"{sample_id}_trimmed_R1.fastq.gz")
        fastp_json = File(f"{sample_id}_fastp.json")
        fastp_html = File(f"{sample_id}_fastp.html")

        fastp_job = Job("fastp")
        fastp_job.add_args(
            "-i", r1_input,
            "-o", trimmed_r1,
            "--json", fastp_json,
            "--html", fastp_html,
            "--thread", str(TOOL_CONFIGS["fastp"]["cores"]),
            "--qualified_quality_phred", "20",
            "--length_required", "50"
        )
        fastp_job.add_inputs(r1_input)
        fastp_job.add_outputs(trimmed_r1, fastp_json, fastp_html, stage_out=True)

        if is_paired:
            trimmed_r2 = File(f"{sample_id}_trimmed_R2.fastq.gz")
            fastp_job.add_args("-I", r2_input, "-O", trimmed_r2)
            fastp_job.add_inputs(r2_input)
            fastp_job.add_outputs(trimmed_r2, stage_out=True)
        else:
            trimmed_r2 = None

        all_qc_files.append(fastp_json)
        wf.add_jobs(fastp_job)

        # ============================================================
        # Step 3: Assembly - MEGAHIT or SPAdes
        # ============================================================
        contigs = File(f"{sample_id}_contigs.fa")
        assembly_log = File(f"{sample_id}_assembly.log")

        if assembler == "megahit":
            assembly_job = Job("megahit")
            if is_paired:
                assembly_job.add_args(
                    "-1", trimmed_r1,
                    "-2", trimmed_r2,
                    "-o", f"{sample_id}_megahit",
                    "-t", str(TOOL_CONFIGS["megahit"]["cores"]),
                    "--min-contig-len", "1000"
                )
                assembly_job.add_inputs(trimmed_r1, trimmed_r2)
            else:
                assembly_job.add_args(
                    "-r", trimmed_r1,
                    "-o", f"{sample_id}_megahit",
                    "-t", str(TOOL_CONFIGS["megahit"]["cores"]),
                    "--min-contig-len", "1000"
                )
                assembly_job.add_inputs(trimmed_r1)
        else:  # spades
            assembly_job = Job("spades")
            if is_paired:
                assembly_job.add_args(
                    "-1", trimmed_r1,
                    "-2", trimmed_r2,
                    "-o", f"{sample_id}_spades",
                    "-t", str(TOOL_CONFIGS["spades"]["cores"]),
                    "--meta"
                )
                assembly_job.add_inputs(trimmed_r1, trimmed_r2)
            else:
                assembly_job.add_args(
                    "-s", trimmed_r1,
                    "-o", f"{sample_id}_spades",
                    "-t", str(TOOL_CONFIGS["spades"]["cores"]),
                    "--meta"
                )
                assembly_job.add_inputs(trimmed_r1)

        assembly_job.add_outputs(contigs, assembly_log, stage_out=True)
        # add_profiles passes the value through unconverted, and the planner
        # wants MB.
        assembly_job.add_profiles(
            Namespace.PEGASUS, key="memory",
            value=str(int(TOOL_CONFIGS[assembler]["memory"].split()[0]) * 1024))
        wf.add_jobs(assembly_job)

        # ============================================================
        # Step 4: Assembly QC - QUAST
        # ============================================================
        quast_report = File(f"{sample_id}_quast_report.tsv")
        quast_html = File(f"{sample_id}_quast_report.html")

        quast_job = Job("quast")
        quast_job.add_args(
            contigs,
            "-o", f"{sample_id}_quast",
            "--min-contig", "1000",
            "--threads", str(TOOL_CONFIGS["quast"]["cores"])
        )
        quast_job.add_inputs(contigs)
        quast_job.add_outputs(quast_report, quast_html, stage_out=True)
        wf.add_jobs(quast_job)

        all_qc_files.append(quast_report)

        # ============================================================
        # Step 5: Gene Prediction - Prodigal
        # ============================================================
        genes_faa = File(f"{sample_id}_genes.faa")
        genes_gff = File(f"{sample_id}_genes.gff")

        prodigal_job = Job("prodigal")
        prodigal_job.add_args(
            "-i", contigs,
            "-a", genes_faa,
            "-o", genes_gff,
            "-f", "gff",
            "-p", "meta"
        )
        prodigal_job.add_inputs(contigs)
        prodigal_job.add_outputs(genes_faa, genes_gff, stage_out=True)
        wf.add_jobs(prodigal_job)

        # ============================================================
        # Step 6: Binning - MetaBAT2
        # ============================================================
        if not skip_binning:
            bins_dir = File(f"{sample_id}_bins")
            depth_file = File(f"{sample_id}_depth.txt")

            # Calculate depth
            depth_job = Job("metabat2")
            depth_job.add_args(
                "jgi_summarize_bam_contig_depths",
                "--outputDepth", depth_file,
                contigs
            )
            depth_job.add_inputs(contigs)
            depth_job.add_outputs(depth_file, stage_out=True)
            wf.add_jobs(depth_job)

            # Run MetaBAT2
            metabat_job = Job("metabat2")
            metabat_job.add_args(
                "-i", contigs,
                "-a", depth_file,
                "-o", f"{sample_id}_bins/bin",
                "-m", "1500",
                "-t", str(TOOL_CONFIGS["metabat2"]["cores"])
            )
            metabat_job.add_inputs(contigs, depth_file)
            metabat_job.add_outputs(bins_dir, stage_out=True)
            wf.add_jobs(metabat_job)

            all_bin_dirs.append(bins_dir)

            # ============================================================
            # Step 7: Bin Quality - CheckM2
            # ============================================================
            checkm2_report = File(f"{sample_id}_checkm2_quality.tsv")

            checkm2_job = Job("checkm2")
            checkm2_job.add_args(
                "predict",
                "--input", bins_dir,
                "--output-directory", f"{sample_id}_checkm2",
                "--threads", str(TOOL_CONFIGS["checkm2"]["cores"])
            )
            if checkm2_db:
                checkm2_job.add_args("--database_path", checkm2_db)
            checkm2_job.add_inputs(bins_dir)
            checkm2_job.add_outputs(checkm2_report, stage_out=True)
            wf.add_jobs(checkm2_job)

            all_qc_files.append(checkm2_report)

            # ============================================================
            # Step 8: Taxonomy - GTDB-Tk
            # ============================================================
            if not skip_taxonomy:
                gtdbtk_summary = File(f"{sample_id}_gtdbtk.summary.tsv")

                gtdbtk_job = Job("gtdbtk")
                gtdbtk_job.add_args(
                    "classify_wf",
                    "--genome_dir", bins_dir,
                    "--out_dir", f"{sample_id}_gtdbtk",
                    "--extension", "fa",
                    "--cpus", str(TOOL_CONFIGS["gtdbtk"]["cores"])
                )
                if gtdbtk_db:
                    gtdbtk_job.add_args("--gtdbtk_data_path", gtdbtk_db)
                gtdbtk_job.add_inputs(bins_dir, checkm2_report)
                gtdbtk_job.add_outputs(gtdbtk_summary, stage_out=True)
                gtdbtk_job.add_profiles(Namespace.PEGASUS, key="memory",
                                        value=TOOL_CONFIGS["gtdbtk"]["memory"])
                wf.add_jobs(gtdbtk_job)

            # ============================================================
            # Step 9: Annotation - Prokka
            # ============================================================
            if not skip_annotation:
                prokka_gff = File(f"{sample_id}_prokka.gff")
                prokka_gbk = File(f"{sample_id}_prokka.gbk")
                prokka_faa = File(f"{sample_id}_prokka.faa")

                prokka_job = Job("prokka")
                prokka_job.add_args(
                    "--outdir", f"{sample_id}_prokka",
                    "--prefix", sample_id,
                    "--metagenome",
                    "--cpus", str(TOOL_CONFIGS["prokka"]["cores"]),
                    bins_dir
                )
                prokka_job.add_inputs(bins_dir)
                prokka_job.add_outputs(prokka_gff, prokka_gbk, prokka_faa, stage_out=True)
                wf.add_jobs(prokka_job)

    # ============================================================
    # Step 10: MultiQC Report
    # ============================================================
    multiqc_report = File("multiqc_report.html")
    multiqc_data = File("multiqc_data.json")

    multiqc_job = Job("multiqc")
    multiqc_job.add_args(
        ".",
        "-o", "multiqc_output",
        "--force"
    )
    for qc_file in all_qc_files:
        multiqc_job.add_inputs(qc_file)
    multiqc_job.add_outputs(multiqc_report, multiqc_data, stage_out=True)
    wf.add_jobs(multiqc_job)

    return wf


# Import Namespace for profiles
from Pegasus.api import Namespace


def main():
    parser = argparse.ArgumentParser(
        description="MAG Workflow Generator for Pegasus WMS",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with samplesheet
  %(prog)s --samplesheet samples.csv --output workflow.yml

  # Run with nf-core/mag test data (auto-downloads)
  %(prog)s --test --output workflow.yml

  # Quick test (skip resource-intensive steps)
  %(prog)s --test --skip-taxonomy --skip-annotation

  # Use SPAdes assembler with custom databases
  %(prog)s --samplesheet samples.csv --assembler spades \\
      --gtdbtk-db /path/to/gtdbtk_db --checkm2-db /path/to/checkm2_db

  # Skip taxonomy and annotation steps
  %(prog)s --samplesheet samples.csv --skip-taxonomy --skip-annotation

  # A Slurm cluster, or a hosted site catalog
  %(prog)s --samplesheet samples.csv --site-style slurm --queue cpu --project my_lab
  %(prog)s --samplesheet samples.csv --hosted-site-catalog unity.yml --project my_lab

Samplesheet format (CSV):
  sample,fastq_1,fastq_2,group
  sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,group1
  sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,group1
        """
    )

    # Test mode
    parser.add_argument(
        "--test", "-t",
        action="store_true",
        help="Use nf-core/mag test data (auto-downloads ~10MB)"
    )

    # Input/Output arguments
    parser.add_argument(
        "--samplesheet", "-s",
        type=str,
        help="Input samplesheet CSV with sample information (not required with --test)"
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="workflow.yml",
        help="Output workflow YAML file (default: workflow.yml)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./output",
        help="Directory for downloaded --test data (default: ./output). "
             "Results stage to output/ next to this script (the 'local' "
             "site in sites.yml)"
    )

    # Assembly options
    parser.add_argument(
        "--assembler",
        type=str,
        choices=["megahit", "spades"],
        default="megahit",
        help="Assembler to use (default: megahit)"
    )

    # Pipeline control
    parser.add_argument(
        "--skip-binning",
        action="store_true",
        help="Skip genome binning steps"
    )
    parser.add_argument(
        "--skip-taxonomy",
        action="store_true",
        help="Skip GTDB-Tk taxonomy classification"
    )
    parser.add_argument(
        "--skip-annotation",
        action="store_true",
        help="Skip Prokka annotation"
    )
    parser.add_argument(
        "--skip-fastqc",
        action="store_true",
        help="Skip FastQC QC reports"
    )

    # Database paths
    parser.add_argument(
        "--gtdbtk-db",
        type=str,
        help="Path to GTDB-Tk database"
    )
    parser.add_argument(
        "--checkm2-db",
        type=str,
        help="Path to CheckM2 database"
    )

    # Execution options. Every site option has a default, so a run with
    # only the inputs plans on "compute": the hosted catalog's site when one
    # is named (--hosted-site-catalog or ~/.pegasusrc), else an HTCondor pool
    # written to sites.yml. dest "execution_site" is what Studio matches to
    # keep its "Where it runs" choice in sync.
    parser.add_argument(
        "--execution-site", "-e",
        dest="execution_site",
        type=str,
        default=HOSTED_SITE,
        help=f"Site to plan against (default: {HOSTED_SITE!r}, the name "
             "hosted catalogs give their site)"
    )
    # -s is already --samplesheet here, so the hosted catalog has no short
    # form.
    parser.add_argument(
        "--hosted-site-catalog",
        metavar="FILE",
        help="Centrally hosted site catalog to plan against, e.g. unity.yml "
             "(github.com/pegasushub/pegasus-site-catalogs); written to "
             "pegasus.properties. Default: the one named in ~/.pegasusrc, "
             "if any."
    )
    parser.add_argument(
        "--site-style",
        choices=("auto",) + STYLES + ("none",),
        default="auto",
        help="auto (default): keep a sites.yml entry or hosted catalog if one "
             "exists, else add an HTCondor site. condor/slurm: (re)write this "
             "site's entry. none: leave sites.yml alone."
    )
    parser.add_argument("--queue", metavar="PARTITION",
                        help="Batch partition/queue jobs submit to")
    parser.add_argument("--project", metavar="ACCOUNT",
                        help="Allocation/account charged on a batch site")
    parser.add_argument("--site-scratch", metavar="DIR",
                        help="Slurm: shared scratch visible to workers and "
                             "the submit host (default: ./work)")
    parser.add_argument("--site-profile", action="append", default=[],
                        type=parse_profile, metavar="NS:KEY=VALUE",
                        help="Extra profile on the execution site; repeatable")
    parser.add_argument("--tag-profile", action="append", default=[],
                        type=parse_tag_profile, metavar="TAG:NS:KEY=VALUE",
                        help="Profile for jobs carrying a tag; repeatable")
    parser.add_argument("--shared-filesystem", choices=("auto", "yes", "no"),
                        default="auto",
                        help="Let jobs read inputs directly from the submit "
                             "host. auto: on for Slurm sites, off for HTCondor.")
    parser.add_argument("--sites-yml", metavar="FILE", default="sites.yml",
                        help="Local site catalog (default: sites.yml)")
    parser.add_argument("--skip-sites-catalog", action="store_true",
                        help="Deprecated: same as --site-style none")
    parser.add_argument(
        "--container-image",
        type=str,
        default=DEFAULT_CONTAINER,
        help=f"Container image to use (default: {DEFAULT_CONTAINER})"
    )

    parser.add_argument(
        "--max-memory-gb",
        type=int,
        default=None,
        help="Cap every job's memory profile at this many GB (for pools "
             "with small worker nodes; test-scale data needs far less than "
             "the production profiles)"
    )
    parser.add_argument(
        "--max-cores",
        type=int,
        default=None,
        help="Cap every job's cores profile (for small worker nodes)"
    )

    args = parser.parse_args()

    # Apply small-pool resource caps before any jobs are built.
    if args.max_memory_gb or args.max_cores:
        for cfg in TOOL_CONFIGS.values():
            if args.max_memory_gb:
                gb = int(cfg["memory"].split()[0])
                cfg["memory"] = f"{min(gb, args.max_memory_gb)} GB"
            if args.max_cores:
                cfg["cores"] = min(cfg["cores"], args.max_cores)

    # Validate input: either --test or --samplesheet must be provided
    if not args.test and not args.samplesheet:
        print("Error: Either --test or --samplesheet must be provided")
        parser.print_help()
        sys.exit(1)

    # Create output directory first (needed for test data download)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    # Handle test mode vs samplesheet mode
    if args.test:
        print("=" * 60)
        print("RUNNING WITH nf-core/mag TEST DATA")
        print("=" * 60)
        samples, samplesheet_path = download_test_data(output_dir)
        if not samples:
            print("Error: Failed to download test data")
            sys.exit(1)
    else:
        # Validate samplesheet exists
        if not os.path.exists(args.samplesheet):
            print(f"Error: Samplesheet not found: {args.samplesheet}")
            sys.exit(1)

        # Parse samplesheet
        print(f"Parsing samplesheet: {args.samplesheet}")
        samples = parse_samplesheet(args.samplesheet)

    if not samples:
        print("Error: No valid samples found")
        sys.exit(1)

    print(f"\nFound {len(samples)} samples:")
    for sample in samples:
        print(f"  - {sample['id']} ({'single-end' if sample['single_end'] else 'paired-end'})")

    # Site catalog and the settings that depend on the site
    print("\nCreating Pegasus catalogs...")
    if args.skip_sites_catalog:
        args.site_style = "none"
    wf_dir = str(Path(__file__).parent.resolve())
    action, style = ensure_sites_yml(
        args.sites_yml, args.execution_site, wf_dir,
        style=args.site_style, hosted=args.hosted_site_catalog,
        queue=args.queue, project=args.project,
        scratch=args.site_scratch, profiles=args.site_profile,
        tag_profiles=args.tag_profile)
    hosted = hosted_catalog(args.hosted_site_catalog)
    print(f"  Site catalog: {args.sites_yml}: {action}"
          + (f" (merged over hosted {hosted})" if hosted else ""))
    if (style is None and hosted and args.execution_site != "local"
            and args.execution_site != HOSTED_SITE):
        # Nothing was written for this site, so planning works only if the
        # hosted catalog happens to define it.
        print(f"Warning: {args.execution_site!r} is not defined in "
              f"{args.sites_yml} and hosted catalogs normally define only "
              f"{HOSTED_SITE!r}: pegasus-plan will fail unless {hosted} has "
              f"it. Use -e {HOSTED_SITE}, or --site-style condor/slurm to "
              f"describe {args.execution_site!r}.")
    if args.shared_filesystem == "auto":
        bypass = style is not None and style != "condor"
    else:
        bypass = args.shared_filesystem == "yes"
    bind = bypass or is_batch_site(style, args.hosted_site_catalog)
    worker = worker_package_url(CONTAINER_PLATFORM)
    if not worker:
        print("Warning: pegasus-version not found: Pegasus will choose the "
              "container's worker package itself (needs curl/wget in the "
              "image and internet on the workers)")

    props = create_properties(args.sites_yml, bypass, worker,
                              args.hosted_site_catalog)
    tc, container = create_transformation_catalog(
        args.container_image, bind_workflow_dir=bind, worker_package=worker)
    rc = create_replica_catalog(samples)

    # Create workflow
    print(f"\nCreating MAG workflow with {args.assembler} assembler...")
    wf = create_workflow(
        samples=samples,
        assembler=args.assembler,
        skip_binning=args.skip_binning,
        skip_taxonomy=args.skip_taxonomy,
        skip_annotation=args.skip_annotation,
        skip_fastqc=args.skip_fastqc,
        gtdbtk_db=args.gtdbtk_db,
        checkm2_db=args.checkm2_db
    )

    # Write catalogs (sites.yml was written by ensure_sites_yml above)
    props.write()
    tc.write("transformations.yml")
    rc.write("replicas.yml")

    # Write workflow
    wf.add_transformation_catalog(tc)
    wf.add_replica_catalog(rc)
    wf.write(args.output)

    print(f"\nWorkflow generated successfully!")
    print(f"  Workflow: {args.output}")
    print(f"  Site catalog: {args.sites_yml}")
    print("  Transformation catalog: transformations.yml")
    print("  Replica catalog: replicas.yml")
    print(f"\nTo submit the workflow:")
    print(f"  pegasus-plan --submit -s {args.execution_site} -o local {args.output}")


if __name__ == "__main__":
    main()
