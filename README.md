# MAG Workflow - Metagenome-Assembled Genomes Pipeline

A Pegasus workflow for metagenomic assembly, binning, and annotation - equivalent to the [nf-core/mag](https://nf-co.re/mag) Nextflow pipeline.

## Overview

This workflow processes metagenomic sequencing data to:
- Quality control and trim raw reads
- Assemble metagenomes using MEGAHIT or SPAdes
- Predict genes from assembled contigs
- Bin contigs into draft genomes (MAGs)
- Assess bin quality and completeness
- Classify bins taxonomically
- Annotate genomes functionally

### Workflow Architecture

```
Input FASTQ files
       │
       ▼
  ┌─────────┐
  │ FastQC  │────────────────────────────────────────────┐
  └────┬────┘                                            │
       ▼                                                 │
  ┌─────────┐                                            │
  │  fastp  │──► Trimmed reads                           │
  └────┬────┘                                            │
       ▼                                                 │
  ┌─────────────┐                                        │
  │MEGAHIT/SPAdes│──► Contigs                            │
  └────┬────────┘                                        │
       │                                                 │
       ├────────────────┐                                │
       ▼                ▼                                │
  ┌─────────┐     ┌──────────┐                           │
  │  QUAST  │     │ Prodigal │──► Gene predictions       │
  └────┬────┘     └──────────┘                           │
       │                                                 │
       ▼                                                 │
  ┌──────────┐                                           │
  │ MetaBAT2 │──► Genome bins                            │
  └────┬─────┘                                           │
       │                                                 │
       ├─────────────────┬─────────────────┐             │
       ▼                 ▼                 ▼             │
  ┌──────────┐     ┌──────────┐     ┌─────────┐          │
  │ CheckM2  │     │ GTDB-Tk  │     │ Prokka  │          │
  └────┬─────┘     └────┬─────┘     └────┬────┘          │
       │                │                │               │
       └────────────────┴────────────────┴───────────────┘
                                         │
                                         ▼
                                  ┌──────────┐
                                  │ MultiQC  │──► Final Report
                                  └──────────┘
```

## Data Source

This workflow processes metagenomic sequencing data from Illumina (short-read) platforms.

### Input Format

Input data is specified via a CSV samplesheet:

```csv
sample,fastq_1,fastq_2,group
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz,group1
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz,group1
sample3,/path/to/sample3_R1.fastq.gz,,group2
```

| Column | Description | Required |
|--------|-------------|----------|
| `sample` | Unique sample identifier | Yes |
| `fastq_1` | Path to forward reads (R1) | Yes |
| `fastq_2` | Path to reverse reads (R2) | No (empty for single-end) |
| `group` | Group ID for co-assembly | No |

### Supported Formats

- **FASTQ files**: `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`
- **Paired-end**: Both R1 and R2 files provided
- **Single-end**: Only R1 file provided (leave R2 empty)

## Features

- **Flexible assembly**: Choose between MEGAHIT (faster, less memory) or SPAdes (more accurate)
- **Comprehensive binning**: MetaBAT2 for robust genome binning
- **Quality assessment**: CheckM2 for bin completeness and contamination estimates
- **Taxonomic classification**: GTDB-Tk for standardized taxonomy
- **Functional annotation**: Prokka for gene annotation
- **Unified reporting**: MultiQC aggregates all QC metrics

## Prerequisites

### Pegasus/HTCondor Cluster

Before running this workflow, you need a Pegasus/HTCondor cluster. See the [FABRIC Pegasus artifact](https://artifacts.fabric-testbed.net/artifacts/53da4088-a175-4f0c-9e25-a4a371032a39) for deployment instructions.

### Required Databases

Some tools require reference databases:

| Tool | Database | Size | Download |
|------|----------|------|----------|
| CheckM2 | DIAMOND database | ~3 GB | `checkm2 database --download` |
| GTDB-Tk | GTDB reference | ~85 GB | [GTDB data](https://gtdb.ecogenomic.org/downloads) |
| Prokka | Built-in | Included | Installed with tool |

**Note**: GTDB-Tk is included in the base container. Use `--skip-taxonomy` if you do not want to run it.

### Software Requirements

- Python 3.9+
- Pegasus WMS v5.0+
- HTCondor v10.2+
- Docker or Singularity

### Python Dependencies

```bash
pip install -r requirements.txt
```

## Quick Start

### Option A: Run with Test Data (Recommended for First-Time Users)

The workflow includes built-in support for nf-core/mag test data. This is the easiest way to verify your setup:

```bash
# Run with test data (auto-downloads ~10MB from nf-core)
./workflow_generator.py --test --skip-fastqc --skip-binning --output workflow.yml

# Quick test (skip resource-intensive taxonomy/annotation steps)
./workflow_generator.py --test --skip-taxonomy --skip-annotation --output workflow.yml

# Submit to HTCondor
pegasus-plan --submit -s condorpool -o local workflow.yml
```

Alternatively, use the standalone fetch script for more control:

```bash
# Download test data separately
./fetch_test_data.py --output-dir ./test_data --samplesheet test_samples.csv

# Then generate workflow
./workflow_generator.py --samplesheet test_samples.csv --output workflow.yml
```

### Option B: Run with Your Own Data

#### 1. Prepare Samplesheet

Create a CSV file with your samples:

```csv
sample,fastq_1,fastq_2,group
gut_sample1,/data/gut1_R1.fastq.gz,/data/gut1_R2.fastq.gz,gut
gut_sample2,/data/gut2_R1.fastq.gz,/data/gut2_R2.fastq.gz,gut
soil_sample1,/data/soil1_R1.fastq.gz,/data/soil1_R2.fastq.gz,soil
```

#### 2. Build Docker Container

```bash
cd Docker
docker build -f MAG_Dockerfile -t kthare10/mag-workflow:latest .
docker push kthare10/mag-workflow:latest
```

#### 3. Generate Workflow

```bash
# Basic usage
./workflow_generator.py \
    --samplesheet samples.csv \
    --output workflow.yml

# With all options
./workflow_generator.py \
    --samplesheet samples.csv \
    --assembler megahit \
    --gtdbtk-db /path/to/gtdbtk_db \
    --checkm2-db /path/to/checkm2_db \
    --execution-site condorpool \
    --output workflow.yml
```

#### 4. Submit Workflow

```bash
pegasus-plan --submit -s condorpool -o local workflow.yml
```

#### 5. Monitor Progress

```bash
pegasus-status <run_directory>
pegasus-analyzer <run_directory>
```

## Command-Line Options

### workflow_generator.py

| Argument | Description | Default |
|----------|-------------|---------|
| `--test, -t` | Use nf-core/mag test data (auto-downloads) | False |
| `--samplesheet, -s` | Input CSV with sample information | Required (unless `--test`) |
| `--output, -o` | Output workflow YAML file | `workflow.yml` |
| `--output-dir` | Output directory for results | `./output` |
| `--assembler` | Assembler: `megahit` or `spades` | `megahit` |
| `--skip-binning` | Skip genome binning steps | False |
| `--skip-fastqc` | Skip FastQC QC reports | False |
| `--skip-taxonomy` | Skip GTDB-Tk classification | False |
| `--skip-annotation` | Skip Prokka annotation | False |
| `--gtdbtk-db` | Path to GTDB-Tk database | None |
| `--checkm2-db` | Path to CheckM2 database | None |
| `--execution-site, -e` | HTCondor execution site | `condorpool` |
| `--container-image` | Docker/Singularity image | `kthare10/mag-workflow:latest` |

### fetch_test_data.py

| Argument | Description | Default |
|----------|-------------|---------|
| `--output-dir, -d` | Directory to download test data | `./test_data` |
| `--samplesheet, -s` | Output samplesheet path | `test_samplesheet.csv` |
| `--include-databases` | Also download mock databases | False |
| `--quiet, -q` | Suppress progress output | False |

## Pipeline Steps

### 1. Quality Control (FastQC)

Generates quality metrics for raw sequencing reads:
- Per-base sequence quality
- GC content distribution
- Adapter content
- Overrepresented sequences

### 2. Read Trimming (fastp)

Filters and trims reads:
- Adapter removal
- Quality filtering (Q20 default)
- Length filtering (50bp minimum)
- Per-read quality reports

### 3. Metagenome Assembly

**MEGAHIT** (default):
- Memory-efficient de Bruijn graph assembler
- Good for large datasets
- Minimum contig length: 1000bp

**SPAdes** (alternative):
- More accurate for complex metagenomes
- Higher memory requirements
- Uses `--meta` mode for metagenomics

### 4. Assembly QC (QUAST)

Evaluates assembly quality:
- Total assembly length
- N50/N90 statistics
- Number of contigs
- GC content

### 5. Gene Prediction (Prodigal)

Identifies protein-coding genes:
- Uses `-p meta` for metagenomic mode
- Outputs GFF coordinates
- Outputs protein sequences (FAA)

### 6. Genome Binning (MetaBAT2)

Clusters contigs into genome bins:
- Uses tetranucleotide frequency
- Incorporates coverage information
- Minimum contig size: 1500bp

### 7. Bin Quality (CheckM2)

Assesses bin quality:
- Completeness estimation
- Contamination estimation
- Quality classification (high/medium/low)

### 8. Taxonomy (GTDB-Tk)

Classifies bins taxonomically:
- Uses GTDB reference database
- Provides standardized taxonomy
- Outputs classification confidence

### 9. Annotation (Prokka)

Annotates genome bins:
- Gene prediction
- Functional annotation
- Multiple output formats (GFF, GBK, FAA)

### 10. Report (MultiQC)

Aggregates all QC metrics:
- Combined HTML report
- Interactive visualizations
- Exportable data tables

## Output Files

```
output/
├── {sample}_R1_fastqc.html          # FastQC reports
├── {sample}_R2_fastqc.html
├── {sample}_fastp.html              # fastp reports
├── {sample}_fastp.json
├── {sample}_trimmed_R1.fastq.gz     # Trimmed reads
├── {sample}_trimmed_R2.fastq.gz
├── {sample}_contigs.fa              # Assembled contigs
├── {sample}_assembly.log
├── {sample}_quast_report.html       # Assembly QC
├── {sample}_quast_report.tsv
├── {sample}_genes.faa               # Predicted proteins
├── {sample}_genes.gff               # Gene coordinates
├── {sample}_bins/                   # Genome bins
│   ├── bin.1.fa
│   ├── bin.2.fa
│   └── ...
├── {sample}_checkm2_quality.tsv     # Bin quality
├── {sample}_gtdbtk.summary.tsv      # Taxonomy
├── {sample}_prokka.gff              # Annotations
├── {sample}_prokka.gbk
├── {sample}_prokka.faa
├── multiqc_report.html              # Combined report
└── multiqc_data.json
```

## Resource Requirements

| Step | Memory | CPUs | Time (per sample) |
|------|--------|------|-------------------|
| FastQC | 2 GB | 2 | ~5 min |
| fastp | 4 GB | 4 | ~10 min |
| MEGAHIT | 16 GB | 8 | ~1-4 hours |
| SPAdes | 32 GB | 16 | ~2-8 hours |
| QUAST | 4 GB | 4 | ~10 min |
| Prodigal | 4 GB | 1 | ~15 min |
| MetaBAT2 | 8 GB | 4 | ~30 min |
| CheckM2 | 16 GB | 8 | ~30 min |
| GTDB-Tk | 64 GB | 8 | ~1-2 hours |
| Prokka | 8 GB | 4 | ~30 min |
| MultiQC | 4 GB | 2 | ~5 min |

## Customization

### Adding/Modifying Tools

To add or modify tools in the workflow:

1. **Add wrapper script** in `bin/` directory:

```bash
#!/bin/bash
set -euo pipefail
echo "=== My Tool ==="
my_tool "$@"
```

2. **Update TOOL_CONFIGS** in `workflow_generator.py`:

```python
TOOL_CONFIGS = {
    # ... existing tools ...
    "mytool": {"memory": "8GB", "cores": 4},
}
```

3. **Add transformation** in `create_transformation_catalog()`:

```python
tools = [
    # ... existing tools ...
    "mytool",
]
```

4. **Add job** in `create_workflow()`:

```python
mytool_output = File(f"{sample_id}_mytool_output.txt")
mytool_job = Job("mytool")
mytool_job.add_args("--input", input_file, "--output", mytool_output)
mytool_job.add_inputs(input_file)
mytool_job.add_outputs(mytool_output, stage_out=True)
wf.add_jobs(mytool_job)
```

### Modifying Parameters

Tool parameters can be modified in the `create_workflow()` function. For example, to change fastp quality threshold:

```python
fastp_job.add_args(
    # ... other args ...
    "--qualified_quality_phred", "30",  # Change from 20 to 30
    "--length_required", "100",          # Change from 50 to 100
)
```

## Troubleshooting

### Common Issues

**1. Out of memory during assembly**
- Use MEGAHIT instead of SPAdes
- Increase memory allocation in TOOL_CONFIGS
- Consider subsampling reads

**2. GTDB-Tk fails with database errors**
- Verify database path with `--gtdbtk-db`
- Check database version compatibility
- Ensure sufficient disk space (~85 GB)

**3. No bins produced by MetaBAT2**
- Check assembly quality (N50 > 10kb recommended)
- Verify sufficient sequencing depth
- Lower minimum contig length

**4. Container pull failures**
- Verify container image exists
- Check Singularity cache permissions
- Try pre-pulling: `singularity pull docker://kthare10/mag-workflow:latest`

### Debugging

```bash
# View workflow status
pegasus-status <run_directory>

# Analyze failures
pegasus-analyzer <run_directory>

# Check job logs
cat <run_directory>/*/*.out
cat <run_directory>/*/*.err

# View Condor job details
condor_q -analyze <job_id>
```

## Comparison with nf-core/mag

| Feature | nf-core/mag | This Workflow |
|---------|-------------|---------------|
| Workflow engine | Nextflow | Pegasus WMS |
| Execution | Cloud/HPC | HTCondor pools |
| Assemblers | MEGAHIT, SPAdes, FLYE | MEGAHIT, SPAdes |
| Binners | MetaBAT2, MaxBin2, CONCOCT | MetaBAT2 |
| Bin refinement | DAS Tool | Not included |
| Virus detection | geNomad | Not included |
| Long-read support | Yes | Planned |

## References

- nf-core/mag: https://nf-co.re/mag
- MEGAHIT: https://github.com/voutcn/megahit
- SPAdes: https://github.com/ablab/spades
- MetaBAT2: https://bitbucket.org/berkeleylab/metabat
- CheckM2: https://github.com/chklovski/CheckM2
- GTDB-Tk: https://github.com/Ecogenomics/GTDBTk
- Prokka: https://github.com/tseemann/prokka
- Pegasus WMS: https://pegasus.isi.edu/

## Authors

Komal Thareja (kthare10@renci.org)

Built with the assistance of [Claude](https://claude.ai), Anthropic's AI assistant.

## License

This workflow is released under the [MIT license](./LICENSE).
