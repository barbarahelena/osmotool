# osmotool

Screen osmoadaptation genes in metagenomic datasets using [DIAMOND](https://github.com/bbuchfink/diamond) and [HMMER](http://hmmer.org/).

## Overview

`osmotool` has three subcommands:

| Subcommand | Input | Method |
|---|---|---|
| `download-db` | — | Downloads and unpacks an [`osmo_refdb`](https://github.com/barbarahelena/osmo_refdb) release from Zenodo |
| `profile` | FASTQ reads (paired or single-end) | `diamond blastx` (6-frame translation), optionally with a DIAMOND+HMM cascade (`--cascade`) |
| `annotate` | Assembly FASTA (or pre-called proteins) | Prodigal → `hmmscan --cut_ga` (default), or `diamond blastp` (`--method diamond`) |

Both modes produce per-gene raw counts and a normalised value: `profile`
reports real depth-normalised RPM abundances; `annotate` reports
`copies_per_kb`, a per-genome gene-copy-number statistic (see
[Output](#output) below) — the two are not directly comparable.

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
docker pull barbarahelena/osmotool:latest

# Singularity (HPC)
singularity pull docker://barbarahelena/osmotool:latest
```

Images are built and pushed to Docker Hub automatically by
[`.github/workflows/docker-publish.yml`](.github/workflows/docker-publish.yml)
whenever a version tag is pushed (see [Versioning](#versioning)) — `:latest`
and `:X.Y.Z` are both updated on each release.

## Reference database

`osmotool` takes a single `DATABASE` argument: the path to an unpacked
[`osmo_refdb`](https://github.com/barbarahelena/osmo_refdb) release
directory (e.g. `v5/`) — download and unpack it with:

```bash
osmotool download-db --release v5 --location /path/to/refdb
# -> downloads + unpacks into /path/to/refdb/v5/
osmotool profile /path/to/refdb/v5 ...
```

`--release` defaults to `latest` (currently `v5`); `--location` defaults to
the current directory. Already-extracted releases are left alone unless
`--force` is passed, and the downloaded archive is deleted after a
successful extraction unless `--keep_archive` is passed. The download is
checksum-verified against Zenodo's reported md5 before extraction.

**v5** — DOI: [10.5281/zenodo.21420253](https://doi.org/10.5281/zenodo.21420253)

Or download manually:

```bash
curl -LO https://zenodo.org/records/21420253/files/v5.tar.gz
tar -xzf v5.tar.gz    # produces a v5/ directory
osmotool profile v5 ...
```

Cite this DOI (alongside `osmotool` itself) in any publication using these
results.

Inside that directory, osmotool finds everything it needs by fixed
filename — no need to point at each file separately:

| File | Used by |
|---|---|
| `osmo_refdb.dmnd` | both (DIAMOND database) |
| `hmms/osmo_refdb.hmm` (+ `.h3*` indices) | `annotate`'s default `--method hmm`/`both`; `profile --cascade` |
| `osmo_refdb.profile_cascade.tsv` | `profile --cascade` |
| `osmo_refdb.diamond_cutoffs.tsv` | `annotate --method diamond`/`both` (auto-applied if present) |
| `osmo_refdb.profile_excluded_families.txt` | `profile` (auto-applied if present) |
| `osmo_refdb.annotate_excluded_families.txt` | `annotate` (auto-applied if present) |

Only `osmo_refdb.dmnd` is required. Everything else is optional and used
automatically when present, except the `profile` DIAMOND+HMM cascade,
which also needs `--cascade` on the command line since it pulls in
`orfm` + `hmmscan` as extra dependencies.

## Usage

### Profile reads (FASTQ)

```bash
osmotool profile \
  /path/to/osmo_refdb/releases/v5 \
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
  /path/to/osmo_refdb/releases/v5 \
  --singles sample.fastq.gz \
  --out_prefix results/sample \
  --threads 8
```

With the DIAMOND+HMM cascade (a second HMM opinion on families with
demonstrated DIAMOND precision problems). Scope-excluded families are
hidden from the report automatically whenever
`osmo_refdb.profile_excluded_families.txt` is present in `DATABASE`:

```bash
osmotool profile \
  /path/to/osmo_refdb/releases/v5 \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --out_prefix results/sample \
  --cascade \
  --threads 8
```

### Annotate assembly

```bash
# HMM (default method) with Prodigal ORF prediction:
osmotool annotate \
  /path/to/osmo_refdb/releases/v5 \
  assembly.fasta \
  --out_prefix results/assembly \
  --threads 8

# DIAMOND instead (avoids the HMMER dependency). The per-family
# specificity gate and decoy references hidden from the report are
# applied automatically whenever osmo_refdb.diamond_cutoffs.tsv /
# osmo_refdb.annotate_excluded_families.txt are present in DATABASE:
osmotool annotate \
  /path/to/osmo_refdb/releases/v5 \
  assembly.fasta \
  --method diamond \
  --out_prefix results/assembly \
  --threads 8

# With pre-called proteins (skip Prodigal):
osmotool annotate \
  /path/to/osmo_refdb/releases/v5 \
  assembly.fasta \
  --proteins assembly_proteins.faa \
  --out_prefix results/assembly \
  --threads 8
```

## Output

| File | Description |
|---|---|
| `<prefix>.gene_counts.tsv` | Per-gene raw counts and normalised value (`rpm` or `copies_per_kb`, depending on mode) |
| `<prefix>.aln_stats.tsv` | Summary statistics |

### `gene_counts.tsv` format

`profile` mode:

```
# total_reads    1000000
# filtered_reads 42753
gene    raw_count   rpm
UniRef90_A0A000   12.5    3.21
UniRef90_A0A001    4.0    1.05
```

**RPM** = `(raw_count × 10⁶) / total_reads` — reads are not assembled, so no gene
length is available and no per-length correction is applied. This is a real
sequencing-depth-normalised abundance estimate: it corrects for how many reads a
sample happened to be sequenced to.

`annotate` mode:

```
# total_proteins 4128
# filtered_reads 187
gene    raw_count   copies_per_kb
ectA    1.0         5.87
ectB    1.0         6.42
```

**copies_per_kb** = `(raw_count × 10⁹) / (gene_length_bp × total_proteins)`, using
protein lengths from Prodigal.

This reuses RPM/RPKM's formula shape, but it is **not** an abundance or expression
estimate — `annotate` works on an assembled genome, not raw reads, so there's no
sequencing depth to normalise against. `total_proteins` is the number of
Prodigal-called ORFs in that one genome, and `raw_count` is that family's copy
number in the genome (usually 0, 1, or a small integer). `copies_per_kb` is a
gene-length- and proteome-size-corrected copy-number statistic — useful for
comparing a family's representation across genomes of different sizes, but **not**
directly comparable to `profile` mode's RPM values, which measure something
different (community-level read abundance).

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
| `--method` | `hmm` (annotate only) | `diamond`, `hmm`, or `both`. `profile` is DIAMOND-only (HMM's GA cutoffs don't apply to short read fragments); `annotate` defaults to HMM since it's meaningfully more specific on full-length ORFs. Requires `hmms/osmo_refdb.hmm` present in `DATABASE` for `hmm`/`both` |
| `--cascade` | off | `profile` only. Give reads whose DIAMOND call lands on a family with demonstrated precision problems a second opinion from `hmmscan`, using `DATABASE`'s `osmo_refdb.profile_cascade.tsv` + `hmms/osmo_refdb.hmm`. Requires `orfm`+`hmmscan` on PATH, and both files present in `DATABASE` |
| `--proteins` | none | `annotate` only. Pre-called protein FASTA (skip Prodigal) |
| `--keep_aln` / `--keep_proteins` | off | Retain intermediate alignment/protein files instead of deleting them after the run |

DIAMOND per-family cutoffs (`osmo_refdb.diamond_cutoffs.tsv`) and
scope-excluded-families reporting filters
(`osmo_refdb.profile_excluded_families.txt` /
`osmo_refdb.annotate_excluded_families.txt`) apply automatically whenever
present in `DATABASE` — no flag needed.

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
    /path/to/osmo_refdb/releases/v5 \
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
    container 'barbarahelena/osmotool:latest'

    input:
    tuple val(sample), path(r1), path(r2)
    path osmo_refdb_dir

    output:
    tuple val(sample), path("${sample}.gene_counts.tsv"), emit: counts
    tuple val(sample), path("${sample}.aln_stats.tsv"),  emit: stats

    script:
    """
    osmotool profile \\
        ${osmo_refdb_dir} \\
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
