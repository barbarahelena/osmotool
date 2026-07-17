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

## Reference database

`osmotool` needs a reference DIAMOND database (`--database`, `.dmnd`) and,
for HMM-based detection (`annotate`'s default, or `profile`'s optional
cascade), a pressed HMM database (`--hmm_db`, `.hmm` + `.h3*` indices).
These are built and benchmarked separately by
[`osmo_refdb`](https://github.com/barbarahelena/osmo_refdb) — download the
latest release from Zenodo:

**v5** — DOI: [10.5281/zenodo.21420253](https://doi.org/10.5281/zenodo.21420253)

Cite this DOI (alongside `osmotool` itself) in any publication using these
results.

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

With the DIAMOND+HMM cascade (a second HMM opinion on families with
demonstrated DIAMOND precision problems) and scope-excluded families
hidden from the report:

```bash
osmotool profile \
  /path/to/osmodiamond.dmnd \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --out_prefix results/sample \
  --hmm_db /path/to/osmo_refdb.hmm \
  --cascade_config /path/to/osmo_refdb.profile_cascade.tsv \
  --exclude_families /path/to/osmo_refdb.profile_excluded_families.txt \
  --threads 8
```

### Annotate assembly

```bash
# HMM (default method) with Prodigal ORF prediction:
osmotool annotate \
  /path/to/osmodiamond.dmnd \
  assembly.fasta \
  --hmm_db /path/to/osmo_refdb.hmm \
  --out_prefix results/assembly \
  --threads 8

# DIAMOND instead (avoids the HMMER dependency), with a per-family
# specificity gate and decoy references hidden from the report:
osmotool annotate \
  /path/to/osmodiamond.dmnd \
  assembly.fasta \
  --method diamond \
  --diamond_cutoffs /path/to/osmo_refdb.diamond_cutoffs.tsv \
  --exclude_families /path/to/osmo_refdb.annotate_excluded_families.txt \
  --out_prefix results/assembly \
  --threads 8

# With pre-called proteins (skip Prodigal):
osmotool annotate \
  /path/to/osmodiamond.dmnd \
  assembly.fasta \
  --proteins assembly_proteins.faa \
  --hmm_db /path/to/osmo_refdb.hmm \
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
| `--min_subject_cover` | 0.0 (off) | Min fraction of the matched *reference* protein a hit must cover. Off by default (a short read is always much shorter than most full-length references) — raise cautiously (e.g. 0.1-0.3) to filter hits that only cover a small fused-in domain of a much longer multi-domain reference protein |
| `--min_seqlen` | 45 | Min aligned length (bp) (`profile`; `--min_seqlen` there is in amino acids — see `osmotool profile --help`) |
| `--evalue` | 1e-5 | Max e-value |
| `--threads` | 4 | Threads for DIAMOND |
| `--tmpdir` | `$TMPDIR` | Temp directory for merged paired-end FASTQ |
| `--total_reads` | *(counted)* | Skip read counting when count is already known |
| `--method` | `hmm` (annotate only) | `diamond`, `hmm`, or `both`. `profile` is DIAMOND-only (HMM's GA cutoffs don't apply to short read fragments); `annotate` defaults to HMM since it's meaningfully more specific on full-length ORFs |
| `--hmm_db` | none | Path to the pressed reference HMM database (`osmo_refdb.hmm` + `.h3*` indices). Required for `--method hmm`/`both` (`annotate`), or together with `--cascade_config` (`profile`) |
| `--cascade_config` | none | `profile` only. Path to `<release>.profile_cascade.tsv`: reads whose DIAMOND call lands on a family with demonstrated precision problems get a second opinion from `hmmscan` before being kept. Requires `--hmm_db`, plus `orfm`+`hmmscan` on PATH |
| `--diamond_cutoffs` | none | `annotate` only. Path to `<release>.diamond_cutoffs.tsv`: per-family minimum DIAMOND bitscore, giving DIAMOND calls a specificity gate analogous to HMM's GA cutoff |
| `--exclude_families` | none | Path to a release's excluded-families list (`profile`: `<release>.profile_excluded_families.txt`; `annotate`: `<release>.annotate_excluded_families.txt` — different files). Those families/labels are still searched normally; this only drops them from the reported `gene_counts.tsv` |
| `--proteins` | none | `annotate` only. Pre-called protein FASTA (skip Prodigal) |
| `--keep_aln` / `--keep_proteins` | off | Retain intermediate alignment/protein files instead of deleting them after the run |

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
