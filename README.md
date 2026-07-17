# osmotool

Screen osmoadaptation genes in metagenomic datasets using [DIAMOND](https://github.com/bbuchfink/diamond).

## Overview

`osmotool` has two subcommands:

| Subcommand | Input | Method |
|---|---|---|
| `profile` | FASTQ reads (paired or single-end) | `diamond blastx` (6-frame translation) |
| `annotate` | Assembly FASTA (or pre-called proteins) | Prodigal → `diamond blastp` |

Both modes produce per-gene raw counts and RPKM abundances.

## Installation

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate osmotool
pip install -e .
```

### pip

```bash
pip install osmotool
# requires diamond and prodigal on PATH
```

### Container (Docker / Singularity)

```bash
# Docker
docker pull ghcr.io/barbarahelena/osmotool:latest

# Singularity (HPC)
singularity pull docker://ghcr.io/barbarahelena/osmotool:latest
```

## Usage

### Profile reads (FASTQ)

```bash
osmotool profile \
  /path/to/osmodiamond.dmnd \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  --out_prefix results/sample \
  --min_identity 0.80 \
  --min_query_cover 0.80 \
  --evalue 1e-5 \
  --threads 8
```

Single-end reads:

```bash
osmotool profile \
  /path/to/osmodiamond.dmnd \
  --singles sample.fastq.gz \
  --out_prefix results/sample \
  --threads 8
```

### Annotate assembly

```bash
# With Prodigal ORF prediction (default):
osmotool annotate \
  /path/to/osmodiamond.dmnd \
  assembly.fasta \
  --out_prefix results/assembly \
  --threads 8

# With pre-called proteins (skip Prodigal):
osmotool annotate \
  /path/to/osmodiamond.dmnd \
  assembly.fasta \
  --proteins assembly_proteins.faa \
  --out_prefix results/assembly \
  --threads 8
```

## Output

| File | Description |
|---|---|
| `<prefix>.gene_counts.tsv` | Per-gene raw counts and RPKM |
| `<prefix>.aln_stats.tsv` | Summary statistics |

### `gene_counts.tsv` format

```
# total_reads    1000000
# filtered_reads 42753
gene    raw_count   rpkm
UniRef90_A0A000   12.5    3.21
UniRef90_A0A001    4.0    1.05
```

**RPKM** = `(raw_count × 10⁹) / (gene_length_bp × total_reads)`

In `profile` mode, gene lengths are not available (reads are not assembled), so RPM is
reported instead of RPKM. In `annotate` mode, protein lengths from Prodigal are used
for true RPKM.

**Paired-end counting**: each read of a pair counts as 0.5 (one read pair = one fragment).

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--min_identity` | 0.80 | Min fractional sequence identity |
| `--min_query_cover` | 0.80 | Min query coverage |
| `--min_seqlen` | 45 | Min aligned length (bp) |
| `--evalue` | 1e-5 | Max e-value |
| `--threads` | 4 | Threads for DIAMOND |
| `--tmpdir` | `$TMPDIR` | Temp directory for merged paired-end FASTQ |
| `--total_reads` | *(counted)* | Skip read counting when count is already known |

## HPC usage

osmotool is designed for use on HPC clusters with SLURM or PBS.

### Installation

```bash
# 1. Create the conda environment (no internet needed on compute nodes)
conda env create -f environment.yml

# 2. Install the package into that environment
conda activate osmotool
pip install -e .        # development
# or: pip install osmotool   (release)
```

> **Note:** `pip install -e .` is intentionally **not** run inside
> `environment.yml` because editable installs require write access to the
> source tree, which may not be available on login vs. compute nodes.

### SLURM example

```bash
#!/bin/bash
#SBATCH --job-name=osmotool_profile
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00

module load anaconda3
conda activate osmotool

osmotool profile \
    /path/to/osmodiamond.dmnd \
    -1 ${SAMPLE}_R1.fastq.gz \
    -2 ${SAMPLE}_R2.fastq.gz \
    --out_prefix results/${SAMPLE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --tmpdir ${TMPDIR} \
    --total_reads ${NREADS}   # pass from MultiQC if available
```

`$TMPDIR` is set by SLURM to fast node-local scratch — merged paired-end
FASTQ lands there instead of on the shared filesystem.

`--total_reads` accepts a pre-computed count (e.g. from FastQC/MultiQC)
to skip the `zcat | wc -l` step on very large files.

## Use in Nextflow

```nextflow
process OSMOTOOL_PROFILE {
    container 'ghcr.io/barbarahelena/osmotool:latest'

    input:
    tuple val(sample), path(r1), path(r2)
    path diamond_db

    output:
    tuple val(sample), path("${sample}.gene_counts.tsv"), emit: counts
    tuple val(sample), path("${sample}.aln_stats.tsv"),  emit: stats

    script:
    """
    osmotool profile \\
        ${diamond_db} \\
        -1 ${r1} -2 ${r2} \\
        --out_prefix ${sample} \\
        --threads ${task.cpus}
    """
}
```

## Versioning

Versions are tied to git tags using `setuptools-scm`:

```bash
git tag v0.1.0
git push --tags
# → package version becomes 0.1.0
```

## Development

```bash
git clone https://github.com/barbarahelena/osmotool
cd osmotool
pip install -e ".[dev]"
pytest
```

## License

MIT
