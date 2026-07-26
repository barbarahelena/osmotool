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
| `<prefix>.systems.tsv` | `annotate` only — multi-subunit system calls (see [Multi-subunit systems](#multi-subunit-systems-annotate-mode-only) below) |

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

## Multi-subunit systems (`annotate` mode only)

Several osmoadaptation gene families are not single-ORF genes but subunits of one
multi-subunit system — e.g. the OpuA glycine-betaine ABC transporter (`opuAA` +
`opuAB` + `opuAC`) or the Firmicutes Mrp/Mnh Na⁺/H⁺ antiporter complex (`mrpA`
through `mrpG`, the Bacillota substitute for `nhaA`). Scoring each subunit gene
independently (as `gene_counts.tsv` does) can't tell you whether a genome actually
has the assembled complex versus just some of its subunits scattered around an
unrelated part of the genome. `annotate` mode adds a proximity/co-occurrence layer
on top of the existing per-gene calls to answer that, written to
`<prefix>.systems.tsv`:

| system | mandatory subunits | accessory | fusion |
|---|---|---|---|
| `opuA` | opuAA, opuAB, opuAC | — | — |
| `opuB` | opuBA, opuBB, opuBC | — | — |
| `opuC` | opuCA, opuCB, opuCC | — | — |
| `mrp` | mrpA, mrpB, mrpC, mrpE | mrpD, mrpF, mrpG | mrpA+mrpB may occur as one fused ORF (`mrpA_mrpB_fused`) |
| `ect` | ectA, ectB, ectC | — | — |
| `ots` | otsA, otsB | — | — |
| `trk` | trkA, trkH | — | — |
| `ktr` | ktrA, (ktrB *or* ktrD) | — | — |
| `galE_murB_mazG` | galE, murB, mazG | — | — |

`galE_murB_mazG` is a different KIND of system from the rest of the table — not a
transcriptional operon, but a wider genomic-neighborhood finding: [Culligan et
al. 2012](https://doi.org/10.1038/ismej.2012.38) (*ISME J* 6(10):1916-1925,
PMID 22534607) found these three genes — each individually too broad a domain
to be specific on its own (see `annotate_only`/`profile_excluded_families`
scoping elsewhere in this README) — together on a single ~39kb metagenomic
fosmid clone from a *Collinsella*-related organism via transposon mutagenesis,
each independently conferring NaCl/KCl tolerance when expressed in an
osmosensitive *E. coli* strain. Confirmed directly against *C. aerofaciens*
ATCC 25986's own reference genome (NCBI GCF_010509075.1): `galE`/`murB` sit 1
gene / 1701bp apart, but `mazG` is a real 23 genes / 20810bp away — with a
strand flip in between. That's why this system alone uses a much wider window
(30 genes / 35000bp, bounded by the fosmid's own insert size) and disables the
strand-consistency check other systems use (`require_strand_consistency=False`
in `systems.py` — a real, confirmed strand flip sits inside the one validated
positive case, and there's no reason to expect one strand across a
~20-40kb neighborhood the way there is for an actual operon). Checked against
all seven other validation genomes: fires only on *C. aerofaciens* itself, not
indiscriminately on every genome carrying all three genes somewhere.

For each system, `systems.tsv` reports:

| Column | Meaning |
|---|---|
| `present` | Every mandatory subunit was found co-located (within the proximity window) on one contig |
| `subunit_count` | How many mandatory+accessory subunits were found *anywhere* in the assembly, regardless of proximity |
| `colocalized` | Whether `present` was proximity-confirmed (mirrors `present`; kept as its own column for clarity) |
| `fragmented_assembly_flag` | Subunits were found, but on ≥2 different contigs, so colocalization genuinely can't be evaluated — distinct from subunits simply being too far apart on the *same* contig, which is a real "not colocalized" call, not a fragmentation call |
| `contig` | The contig the system was called present on (`-` otherwise) |
| `mandatory_hits` / `accessory_hits` | Which families were actually found |
| `rescued_hits` | `<family><-original_label>` entries for any subunit reattributed via paralog-confusion rescue (see below) instead of found under its own label |

**Proximity window**: a system counts as colocalized when the mandatory subunits,
sorted by genomic position, form one contiguous run — every *consecutive* pair
within `max_distance_genes` intervening genes (default 3) **and**
`max_distance_bp` (default 2000), **and** strand-consistent (see below). This is
a chain constraint (adjacent pairs), not "every subunit within the window of
every other subunit" — that clique version incorrectly rejects real operons
where a non-mandatory subunit sits between two mandatory ones (e.g. `mrp`'s
accessory `mrpD` sitting directly between mandatory `mrpC` and `mrpE`).

These defaults replace the extension spec's untested 8-gene/8000bp suggestion —
checked instead against real complete-genome operons across six species:
`mrpA`-`G` (*B. subtilis* 168, *S. aureus* NCTC 8325, *L. monocytogenes* EGD-e,
*H. elongata* DSM 2581), `opuABC` (*B. subtilis*, *L. monocytogenes*), `opuCABC`
(*B. subtilis*, *S. aureus*, *L. monocytogenes*, *E. faecalis* V583), `ectABC`
(*H. elongata* — the actual ectoine-discovery organism), and `otsBA`
(*E. coli* K-12) are **all** 0 intervening genes apart (gaps of -34 to +93bp,
essentially back-to-back CDSs) regardless of species. STRING's independently
precomputed gene-neighborhood score for *B. subtilis* `mrpA`/`mrpB` corroborates
this (0.803/1.0). MacSyFinder/TXSScan's own expert-curated models
([github.com/macsy-models](https://github.com/macsy-models)) converge on a
similar range for comparably-sized, similarly tightly-packed operons —
CasFinder's CRISPR-Cas models (7-9 genes, `mrp`'s scale) use
`inter_gene_max_space=3`; TXSScan's small transport systems like T1SS/T2SS
(3-4 genes, `opuA`/`ots`'s scale) use 5. Only genuinely large, multi-locus
machinery (T6SS, the flagellum — nothing in this panel is that complex) goes
up to 20.

The same sweep also surfaced a real false-bridging risk the spec's original
window didn't guard against: several species carry a duplicate/paralogous copy
of the same operon nearby, always on the **opposite strand** — e.g. *H.
elongata*'s second `mrp`-like cluster sits immediately adjacent to the first (0
genes / 236bp), comfortably inside even the *tightened* window. Distance alone
can't reliably reject that — which is exactly what the strand-consistency check
below is for. 3 genes / 2000bp (matching CasFinder's precedent exactly) keeps a
solid margin above every real operon gap observed across all six species —
still not exhaustive, so treat it as a solid starting point rather than final
(see `--systems_config` below to override per system without touching code).

**Strand consistency**: a real operon is transcribed as one contiguous block on
a single strand, so alongside the distance window, every consecutive pair of
chosen subunits must (a) share a strand, and (b) have no gene *of any kind* on
the opposite strand sitting between them. This is calibration-free (binary, no
threshold to pick) and directly motivated by the duplicate-locus finding above
— checked against all six species with zero false rejections of a real operon.
It's a genuine complement to, not a substitute for, the distance check: e.g.
*E. coli*'s `trkA`/`trkH` (hundreds of genes apart, not a real operon) happen to
share a strand too, so distance is what actually rejects that case.

**Fusion handling**: some lineages carry `mrpA`+`mrpB` as a single fused ORF.
osmo_refdb ships a dedicated `mrpA_mrpB_fused` DIAMOND reference for this — a hit
against it satisfies both subunits from one ORF. There's no equivalent fused HMM
profile, so this is currently a **DIAMOND-only** capability; `--method hmm` will
undercount the fused case (same limitation as any single-best-hit-per-ORF scoring
scheme applied to a genuine multi-domain fusion protein) for the lineages that
actually carry the fusion — most don't (e.g. *Bacillus* itself always keeps
`mrpA`/`mrpB` split).

**Paralog-confusion rescue**: some osmo_refdb families are close enough paralogs
of each other, sharing the same Pfam domain, that families.yaml documents them
as systematically swapped for each other by DIAMOND/HMM — not just theoretically
similar. When a mandatory subunit has zero hits under its own family label, but
a hit under one of its documented siblings is positionally confirmed adjacent to
another real, already-found subunit of the SAME system (same proximity window +
strand-consistency check as everything else — never on sequence label alone),
it's reattributed to the correct slot and reported in a new `rescued_hits`
column (e.g. `opuBC<-opuCC`). This directly targets a real finding from the
six-species sweep: in 4/5 Firmicutes tested, the gene right after a confirmed
`opuBB` is a real ORF mislabeled `opuCC` — exactly the substrate-binding-domain
(SBD, PF04069) confusion — and rescuing it is now visible in the output even
when the system overall still isn't `present` for an unrelated reason (see
below).

**Important scope limit, found the hard way**: "documented cross-hit risk
between two labels" and "these are interchangeable copies of the same
biological thing" are not the same claim. `opuAC`/`opuBC`/`opuCC` really are
three structurally-identical, co-occurring *paralogous copies* of the same
recurring operon architecture (confirmed directly — *B. subtilis* carries real,
separate `opuA` and `opuC` operons in the same genome), so rescuing among them
is legitimate. `proX` was initially included in the same rescue clique, since
families.yaml's negative-query design cross-excludes it from the same reference
sets — but `proX` belongs to a different, *deliberately kept separate* system
(ProU), not a duplicate copy of Opu (the same reasoning that got `proV`/`proW`
deferred entirely rather than folded into `opuAA`/`opuAB`). That distinction
was missed initially: rescuing `opuAC` from a real `proX` hit fired on *E.
coli* and was first reported here as a validation win. It wasn't — E. coli's
real ProVWX operon getting reattributed to `opuAC` is a false positive, not
confirmation, and the same mechanism would very plausibly fire the same way on
any Proteobacterium carrying a real ProU operon. `proX` is now explicitly
excluded from the rescue clique (`systems.PARALOG_CONFUSION_CLIQUES`), with a
regression test guarding against reintroducing it. The one other rescue clique,
pore-forming `trkH`/`ktrB`/`ktrD` (PF02386), is treated differently — documented
as lineage-specific naming for the *same* underlying Trk/Ktr system rather than
two independently-arisen analogous ones — but that's this module's own reading
of the documentation, not functionally confirmed either, and worth an
independent check by whoever has closer domain knowledge of osmo_refdb's own
family curation.

**Known refdb-level limitation (`opuB`)**: across all six validation species,
`opuB` never once resolved to a full colocalized operon, even with the rescue
above — always 1-2 of its 3 mandatory subunits. The rescue fixes the `opuBC`
side (a real labeling ambiguity), but the other blocker, `opuBA`, is a different
problem: in every species checked, there's no hit *at all*, under any family, at
the position where it should be. That's a genuine HMM/DIAMOND recall gap in
osmo_refdb's `opuBA` model, not a mislabeling — no positional reasoning can
invent evidence that isn't there. Worth flagging back to osmo_refdb as a
separate, distinct finding from the SBD confusion.

Relevant flags (see [Parameters](#parameters)): `--skip_systems` to disable this
layer entirely, `--systems_config` to override the proximity window per system.

**Not available in `profile` mode**: read-based profiling has no contig-level
coordinates at typical short-read lengths, so colocalization simply isn't
checkable there. `profile` reports each gene family's abundance independently
(unchanged) and carries an explicit `note_systems` line in `aln_stats.tsv` saying
so — a system-level "present" claim requires `annotate` mode (or a dedicated
assembly+binning step) on the same samples.

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
| `--skip_systems` | off | `annotate` only. Skip the [multi-subunit systems](#multi-subunit-systems-annotate-mode-only) layer entirely — no `<prefix>.systems.tsv`, no ORF-coordinate parsing |
| `--systems_config` | none | `annotate` only. TSV overriding `max_distance_genes`/`max_distance_bp` per system (columns: `system`, `max_distance_genes`, `max_distance_bp`) — see [Multi-subunit systems](#multi-subunit-systems-annotate-mode-only) |

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
