# osmotool

Screen osmoadaptation genes in metagenomic datasets using [DIAMOND](https://github.com/bbuchfink/diamond) and [HMMER](http://hmmer.org/).

## Overview

`osmotool` has three subcommands:

| Subcommand | Input | Method |
|---|---|---|
| `download-db` | — | Downloads and unpacks an [`osmo_refdb`](https://github.com/barbarahelena/osmo_refdb) release from Zenodo |
| `profile` | FASTQ reads (paired or single-end) | Reads assigned against the release's annotated gene catalog (bwa-mem3); the rest by translated search (`diamond blastx` + a DIAMOND+HMM cascade). Needs osmo_refdb v10+ |
| `annotate` | Assembly FASTA (or pre-called proteins) | Prodigal → `hmmsearch --cut_ga` (default), or `diamond blastp` (`--method diamond`) |

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

# for NaCl-tolerance prediction (annotate --predict_nacl, predict-nacl):
pip install 'osmotool[predict]'   # adds xgboost
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
directory (e.g. `v10/`) — download and unpack it with:

```bash
osmotool download-db --release v10 --location /path/to/refdb
# -> downloads + unpacks into /path/to/refdb/v10/
osmotool profile /path/to/refdb/v10 ...
```

**`profile` needs v10 or later.** v10 is the first release that ships a gene
catalog, and `profile` assigns reads against it first; on an older release it
stops with an error rather than quietly running translated search alone.
v5–v9 are kept downloadable as **annotate-only** releases: `annotate` results
on them stay valid, and `download-db` warns when you fetch one.

> v10's Zenodo record is not up yet, so `download-db` does not list it and
> `latest` still means v9. Until then, unpack v10 by hand and pass its path.

`--release` defaults to `latest` (currently `v9`); `--location` defaults to
the current directory. Already-extracted releases are left alone unless
`--force` is passed, and the downloaded archive is deleted after a
successful extraction unless `--keep_archive` is passed. The download is
checksum-verified against Zenodo's reported md5 before extraction.

**v9** (latest on Zenodo; annotate-only) — DOI: [10.5281/zenodo.24154957](https://doi.org/10.5281/zenodo.24154957).
Adds `teaA` (ectoine/hydroxyectoine TRAP binding protein, TeaA/UehA clade)
and ships a `markers/` set of universal single-copy genes, which `profile`
uses for copies-per-genome normalisation. Every family's benchmark negatives
are now held out from cutoff calibration, and capped UniProt fetches sample
the whole population instead of the first accessions in sort order; this
lowered the ktrB and opuAC GA cutoffs and raised their genome-level recall
(ktrB called in 0.18 → 0.43 of BacDive genomes, with the ktrA partner present
in 0.96 of the new calls). A peptide ABC ATPase HMM decoy (`pepABC_decoy`)
now competes with opuBA in `annotate`'s best-hit contest and is never
reported. opuBA remains the least reliable family: low recall (HMM 0.13 on
held-out reads) and a high random-background false-positive rate, so treat
its calls with caution (see [the opuB limitation](#multi-subunit-systems-annotate-mode-only)
below). Benchmark precision is not directly
comparable to v8's, because the held-out negative pool is about half the size.

**v8** — DOI: [10.5281/zenodo.21691253](https://doi.org/10.5281/zenodo.21691253).
Calibration-only fix on top of v7's panel: corrects negative-pool
mislabeling for 8 families (murB, mrpB, mrpE, mrpF, mrpG, gshB, cspA,
mscS) and replaces 6 families' Pfam default HMM cutoffs with locally
calibrated ones, improving short-fragment-read sensitivity (e.g. gshB F1
0.47→0.88, mrpE F1 0.20→0.72; panel-wide DIAMOND F1 0.827→0.832, HMM F1
0.693→0.722). No gene/family panel changes. Annotate-only
(`--release v8`).

**v7** — DOI: [10.5281/zenodo.21600046](https://doi.org/10.5281/zenodo.21600046).
Adds the Firmicutes compatible-solute transporters (opuA/opuB/opuC), the
Mrp/Mnh Na⁺/H⁺ antiporter complex, gut-validated subsystems (bet/gsh/cspA),
and Kdp backup/disambiguation genes (trk/ktr) on top of v5's original panel
— see [Multi-subunit systems](#multi-subunit-systems-annotate-mode-only)
above for what `annotate` mode does with the new multi-subunit ones.
Annotate-only (`--release v7`).

**v5** — DOI: [10.5281/zenodo.21420253](https://doi.org/10.5281/zenodo.21420253).
Annotate-only (`--release v5`), for reproducing results generated against
the original panel.

Or download manually:

```bash
curl -LO https://zenodo.org/records/24154957/files/osmo_refdb_v9.tar.gz
tar -xzf osmo_refdb_v9.tar.gz    # produces an osmo_refdb_v9/ directory
osmotool annotate osmo_refdb_v9 ...   # v9 is annotate-only
```

Cite the DOI of whichever release you used (alongside `osmotool` itself) in
any publication using these results.

Inside that directory, osmotool finds everything it needs by fixed
filename — no need to point at each file separately:

| File | Used by |
|---|---|
| `catalogs/<name>/` | `profile` (the gene catalog; v10+, required) |
| `osmo_refdb.dmnd` | both (DIAMOND database) |
| `hmms/osmo_refdb.hmm` (no `hmmpress` indices needed) | `annotate`'s default `--method hmm`/`both`; `profile`'s cascade (on by default) |
| `osmo_refdb.profile_cascade.tsv` | `profile`'s cascade (on by default) |
| `osmo_refdb.diamond_cutoffs.tsv` | `annotate --method diamond`/`both` (auto-applied if present) |
| `osmo_refdb.profile_excluded_families.txt` | `profile` (auto-applied if present) |
| `osmo_refdb.annotate_excluded_families.txt` | `annotate` (auto-applied if present) |

`osmo_refdb.dmnd` is required, and `profile` also needs `catalogs/`. Everything else is optional and used
automatically when present, including the `profile` DIAMOND+HMM cascade
— it's on by default whenever `hmms/osmo_refdb.hmm` and
`osmo_refdb.profile_cascade.tsv` are both present in `DATABASE`, since
several families run as low as 5-25% precision on DIAMOND alone (see
osmo_refdb's own benchmark). It also needs `orfm` + `hmmsearch` on PATH; if
either file is missing or those tools aren't installed, `profile` fails
fast with an error rather than silently reporting DIAMOND-only numbers —
pass `--no_cascade` to opt out and run DIAMOND-only intentionally.

## Usage

### Profile reads (FASTQ)

Reads are assigned against the release's gene catalog first; whatever the
catalog cannot settle goes to translated search, which runs the DIAMOND+HMM
cascade by default (a second HMM opinion on families with demonstrated
DIAMOND precision problems). Needs `bwa-mem3`, `diamond`, `orfm` and
`hmmsearch` on PATH; see [Gene catalog](#gene-catalog-profile-mode) and
[Reference database](#reference-database). Scope-excluded families are
hidden from the report automatically whenever
`osmo_refdb.profile_excluded_families.txt` is present in `DATABASE`:

```bash
osmotool profile \
  /path/to/osmo_refdb/releases/v10 \
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
  /path/to/osmo_refdb/releases/v10 \
  --singles sample.fastq.gz \
  --out_prefix results/sample \
  --threads 8
```

Translated search without the cascade, e.g. when `orfm`/`hmmsearch` aren't installed:

```bash
osmotool profile \
  /path/to/osmo_refdb/releases/v10 \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --out_prefix results/sample \
  --no_cascade \
  --threads 8
```

### Annotate assembly

```bash
# HMM (default method) with Prodigal ORF prediction:
osmotool annotate \
  /path/to/osmo_refdb/releases/v9 \
  assembly.fasta \
  --out_prefix results/assembly \
  --threads 8

# DIAMOND instead (avoids the HMMER dependency). The per-family
# specificity gate and decoy references hidden from the report are
# applied automatically whenever osmo_refdb.diamond_cutoffs.tsv /
# osmo_refdb.annotate_excluded_families.txt are present in DATABASE:
osmotool annotate \
  /path/to/osmo_refdb/releases/v9 \
  assembly.fasta \
  --method diamond \
  --out_prefix results/assembly \
  --threads 8

# With pre-called proteins (skip Prodigal):
osmotool annotate \
  /path/to/osmo_refdb/releases/v9 \
  assembly.fasta \
  --proteins assembly_proteins.faa \
  --out_prefix results/assembly \
  --threads 8

# Also predict the genome's maximum NaCl tolerance (osmo_refdb v10+, needs
# xgboost -- see "NaCl tolerance prediction" below):
osmotool annotate \
  /path/to/osmo_refdb/releases/v10 \
  assembly.fasta \
  --predict_nacl \
  --out_prefix results/assembly \
  --threads 8
```

### Predict NaCl tolerance from existing annotate output

```bash
# one or more annotate --out_prefix values (made with --method hmm and the
# same release); writes <prefix>.nacl_prediction.tsv and
# <prefix>.nacl_contributions.tsv next to each, plus one summary table
osmotool predict-nacl /path/to/osmo_refdb/releases/v10 \
  results/genomeA results/genomeB \
  -o results/nacl_predictions.tsv
```

## Output

| File | Description |
|---|---|
| `<prefix>.gene_counts.tsv` | Per-gene raw counts and normalised value (`rpm` or `copies_per_kb`, depending on mode) |
| `<prefix>.aln_stats.tsv` | Summary statistics, one `metric value` row each (see [`aln_stats.tsv` format](#aln_statstsv-format-profile-mode) below) |
| `<prefix>.catalog_hits.tsv` | `profile --keep_aln` only — every catalog-settled read: gene, call/drop/ambiguous, label, identity, aligned bp, number of valid placements |
| `<prefix>.unit_counts.tsv` | `profile` only, when DATABASE ships `osmo_refdb.profile_report_levels.tsv` (osmo_refdb v10+) — the counts per reporting unit, excluded units included (see [Report units](#report-units-profile-mode) below) |
| `<prefix>.markers.tsv` | `profile` only, when DATABASE ships `markers/` — per-marker evidence behind `copies_per_genome` (see [Copies per genome](#copies-per-genome-profile-mode) below) |
| `<prefix>.systems.tsv` | `annotate` only — multi-subunit system calls (see [Multi-subunit systems](#multi-subunit-systems-annotate-mode-only) below) |
| `<prefix>.gene_coordinates.tsv` | `annotate` only — per-gene contig/start/end/strand coordinates (see [`gene_coordinates.tsv` format](#gene_coordinatestsv-format-annotate-mode-only) below) |
| `<prefix>.nacl_prediction.tsv` | `annotate --predict_nacl` / `predict-nacl` only — predicted maximum NaCl tolerance with interval and flags (see [NaCl tolerance prediction](#nacl-tolerance-prediction-annotate-mode) below) |
| `<prefix>.nacl_contributions.tsv` | same — every feature's contribution to that prediction |

### `gene_counts.tsv` format

`profile` mode:

```
# total_fragments    1000000
# filtered_fragments 42753.5
gene    raw_count   RPM     catalog_count   catalog_ambiguous   translated_count
betA    12.5        12.5    10.0            0.50                2.5
ectA     4.0         4.0     0.0            0.00                4.0
```

With `--translated_only` the three source columns are left out. A
`copies_per_genome` column follows them when `DATABASE` ships `markers/`.

**RPM** = `(raw_count × 10⁶) / total_fragments` — reads are not assembled, so no gene
length is available and no per-length correction is applied. This is a real
sequencing-depth-normalised abundance estimate: it corrects for how many reads a
sample happened to be sequenced to.

**Counts are in fragments**, for both read layouts. A fragment is a read pair
(paired-end) or a read (single-end). Each mate of a pair counts as 0.5, and
`total_fragments` is the number of pairs (R1's read count), so a paired-end
library and the same library run single-end give the same RPM.
`filtered_fragments` (fragments with any alignment) is weighted the same way.

> **Changed in v0.9.2:** earlier versions divided paired-end counts by the
> reads in R1 + R2, so paired-end RPM was half of what single-end gave, and
> `pct_reads_assigned` could not exceed 50%. Paired-end RPM is now 2× its old
> value; single-end is unchanged. The header lines and `aln_stats.tsv` keys
> were renamed (`total_reads` → `total_fragments`, `filtered_reads` →
> `filtered_fragments`, `reads_assigned` → `fragments_assigned`,
> `pct_reads_assigned` → `pct_fragments_assigned`) so old and new outputs
> can't be mixed by accident. `copies_per_genome` is a ratio to marker genes
> counted the same way and is unaffected.

### `aln_stats.tsv` format (`profile` mode)

Two columns, `metric` and `value`, one row per metric. Counts in fragments
(see above) except where a row says reads.

| Metric | Meaning |
|---|---|
| `mode`, `db` | `profile`, and the DIAMOND database path |
| `total_fragments` | Sequenced fragments (read pairs, or single-end reads): the RPM denominator |
| `filtered_fragments` | Fragments with a catalog call or any DIAMOND alignment |
| `fragments_assigned`, `pct_fragments_assigned` | Fragments counted for a reported family, and their share of `filtered_fragments` |
| `min_identity`, `min_query_cover`, `min_subject_cover`, `evalue` | DIAMOND filters used |
| `note_systems` | Reminder that multi-subunit systems are not checked for co-location from reads |
| `genome_coverage_rpk`, `markers_used`, ... | Marker normalisation summary, when `DATABASE` ships `markers/` (see [Copies per genome](#copies-per-genome-profile-mode)) |
| `report_levels`, `units_reported`, `units_excluded` | Report-levels file and its outcome, when shipped |
| `catalog` | Catalog name (`gmgc_gut`), or `none (--translated_only)` |
| `catalog_release` | osmo_refdb release whose HMMs annotated it |
| `catalog_hmm_sha256` | Checksum of those HMMs, as recorded in `catalog.json` and checked against `DATABASE` |
| `catalog_dir`, `catalog_aligner`, `catalog_index` | Where the catalog was read from, the mapping executable, and the index kind (`all`/`hits`/`neighbours`) |
| `catalog_reads_called` | **Reads** (each mate counts 1) the catalog settled as one family |
| `catalog_reads_ambiguous` | Reads with valid placements on different families: settled, split over them in `catalog_ambiguous` |
| `catalog_reads_dropped`, `catalog_dropped_by_reason` | Reads settled on a gene that is not counted (background, homolog below the cutoff, decoy), by reason |
| `catalog_reads_fallback` | Reads the catalog could not settle; these went on to translated search |
| `catalog_settled_fraction` | Share of panel-assigned fragments the catalog settled: Σ(`catalog_count` + `catalog_ambiguous`) / Σ(those + `translated_count`) over the reported families. `NA` when nothing was assigned |
| `translated_fallback` | `yes`, or `no` with `--catalog_only` |
| `count_<family>` | `raw_count` per family |

`catalog_settled_fraction` is counted over panel assignments, not over all
reads: a release's index holds only panel genes and their close neighbours,
so most reads of any sample map nowhere in it, and `catalog_reads_fallback`
is large even when the catalog fits. When the fraction is below 0.5,
`profile` warns on stderr, naming the catalog. A low fraction usually means
the sample comes from a habitat the catalog does not cover (marine or soil
reads against `gmgc_gut`): the counts are still valid, but they are mostly
translated-search counts. The 0.5 threshold is a starting point, not a
calibrated value.

### Copies per genome (`profile` mode)

RPM still depends on what else is in the sample: host DNA, eukaryotes, viruses
and the community's average genome size all shift it. When the release
directory contains a `markers/` subdirectory (osmo_refdb's `run_pipeline.sh
markers`), `profile` adds a `copies_per_genome` column:

```
gene    raw_count   RPM       copies_per_genome
betA    310.5       31.05     0.412
ectC    12.0        1.2       0.087
```

**copies_per_genome** = family reads per kb of effective gene length ÷ mean
reads per kb of ten universal single-copy marker genes (the fetchMG/mOTUs COGs:
ychF, pheS, argS, serS, cysS, leuS, valS, tsaD, ffh, ftsY) — the average number
of copies per genome in the community (for an at-most-single-copy family, the
fraction of cells carrying it). Following MicrobeCensus (Nayfach & Pollard 2015):

- Markers are searched in a **second, separate DIAMOND pass** with the same
  `--min_identity`/`--min_query_cover`/`--evalue` thresholds, so they never
  compete with osmo families for best hits.
- Gene length is the hit-weighted mean length of the reference proteins reads
  aligned to, converted to an **effective length** `3·L_aa − (2c − 1)·R + 1`
  (R = median read length, c = `--min_query_cover`): the number of read
  positions that can produce a qualifying hit. Without it short families
  (cspA, ~70 aa) would be underestimated ~1.7× against long markers at 150 bp.
- Markers with zero reads, and markers ≥ 1 scaled MAD from the median (e.g. one
  duplicated in a dominant lineage), are excluded before averaging. Fewer than 5
  detected markers → `NA`.

`<prefix>.markers.tsv` lists each marker's count, reference length, RPK and
whether it was used, excluded as an outlier, or undetected;
`aln_stats.tsv` gains the summary (`genome_coverage_rpk`, `markers_used`, ...).
Releases without `markers/` produce exactly the old output; `--no_markers`
skips the extra pass.

Limits: per-family detection sensitivity isn't corrected for (both sides use
the same thresholds, so it partly cancels) — compare `copies_per_genome`
across samples rather than reading it as an exact copy number. Run on
host-depleted reads: host mitochondrial aminoacyl-tRNA synthetases are
bacterial-like and could add marker hits.

`annotate` mode:

```
# total_proteins 4128
# genome_bp      4215606
# gc             0.43514
# n_contigs      1
# contig_n50     4215606
# filtered_reads 187
# proteins_with_hits 187
gene    raw_count   copies_per_kb
ectA    1.0         5.87
ectB    1.0         6.42
```

The four lines after `total_proteins` describe the input assembly, so genome size
and GC can be used as covariates when comparing genomes without a separate step:

- `genome_bp`: total assembly length (every character in sequence lines, `N`s included)
- `gc`: (G+C) / (A+C+G+T), five decimals; `N` and other IUPAC ambiguity codes are
  left out of both numerator and denominator (`NA` if the assembly has no A/C/G/T)
- `n_contigs`: number of FASTA records
- `contig_n50`: N50 of record lengths

They are computed from `ASSEMBLY.FASTA` even when `--proteins` is given. `profile`
mode does not write them. Parsers should skip `#` lines they don't recognise.

`proteins_with_hits` is the number of proteins with at least one hit that passed
filtering (DIAMOND identity/coverage/per-family cutoffs, or the HMM gathering
threshold), counted before each protein's best hit is chosen and before decoy
families are dropped. It can therefore be slightly higher than the sum of
`raw_count`. `filtered_reads` holds the same number; it is the name profile mode
used before it counted fragments, and is kept so existing parsers keep working.
Prefer `proteins_with_hits` in new code.

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

**Paired-end counting** (`profile`): each read of a pair counts as 0.5 (one read pair = one
fragment), and RPM's denominator is pairs; see [`gene_counts.tsv` format](#gene_countstsv-format).

### Report units (`profile` mode)

Some panel genes are close paralogs that a short read cannot reliably tell
apart (e.g. the OpuA/OpuB/OpuC permeases), and some are not reliable on short
reads at any level. osmo_refdb decides, per family or paralog group and on its
calibration reads (never its benchmark reads), at which level a count can be
trusted, and ships the decision as `osmo_refdb.profile_report_levels.tsv`.
When the release has it, `profile` writes `<prefix>.unit_counts.tsv` next to
`gene_counts.tsv`:

```
# total_fragments    1000000
# filtered_fragments 42753.5
# report_levels      releases/v10/osmo_refdb.profile_report_levels.tsv
unit          report_level  members            status    raw_count  RPM    [copies_per_genome]
ectA          gene          ectA               reported  12.5       12.5
opu_permease  group         opuAB,opuBB,opuCB  reported  4.0        4
opuAC         none          opuAC              excluded  NA         NA
```

- `gene`: the family on its own, same count as in `gene_counts.tsv`.
- `group`: the members' summed counts. Every read counts toward one family
  only (its best hit), so this is the number of reads called as any member.
  `copies_per_genome` is summed too (`NA` if any member's is unknown).
- `none`: listed with `status excluded` and `NA`, never 0. Don't read an
  excluded unit as "gene absent": the table says nothing about it.
- A label that has reads but is in no unit (e.g. a fused-ORF reference such
  as `mrpA_mrpB_fused`) gets `status not_in_manifest`, so no read disappears
  silently.

The levels were decided on translated-search calibration reads, so they
describe step 2 of `profile`; catalog-settled reads are full-length calls.

Every unit in the manifest gets a row, zeros included. Read results from this
table; `gene_counts.tsv` is unchanged and still lists every family on its own,
including families whose gene-level count is not reliable. The manifest's
other columns hold the calibration numbers behind each choice.

### Gene catalog (`profile` mode)

Translated search on 150 bp reads has to tell close paralogs apart from a
fragment, which is where it fails (see the report levels above). `profile`
therefore works in two steps:

1. **Catalog.** Reads are assigned against the gene catalog the release
   ships, annotated once at full length with the release's HMMs. A read that
   matches a catalog gene near-exactly takes that gene's family, the approach
   of catalog-based profilers such as Cayman.
2. **Translated search.** Whatever the catalog cannot settle falls through
   to DIAMOND and the cascade, as described above.

A release ships its catalogs under `catalogs/<name>/` (index, annotation
and metadata built with that release's HMMs, versioned with it as one unit;
osmo_refdb's `RELEASE_CATALOG`). v10 ships `gmgc_gut` alone, so nothing has
to be passed. Once a release ships catalogs for several habitats, choose one
with `--habitat NAME`:

```bash
osmotool profile releases/v10 -1 R1.fq.gz -2 R2.fq.gz -o sample
osmotool profile releases/v10 -1 R1.fq.gz -2 R2.fq.gz --habitat gmgc_gut -o sample
```

For a catalog of your own, build it once per release (translate, hmmsearch
with the release's HMMs under their gathering cutoffs, index the catalog
with bwa-mem3) and pass its directory with `--catalog_dir` (left out of
`--help`; a shipped catalog is versioned with its release, a local one is
yours to keep in step):

```bash
osmotool build-catalog releases/v10 --genes my_genes.fna.gz -o my_catalog_v10 -t 32
osmotool profile releases/v10 -1 R1.fq.gz -2 R2.fq.gz --catalog_dir my_catalog_v10 -o sample
```

For benchmarking, `--translated_only` skips step 1 and assigns every read by
translated search; it also runs on releases without a catalog.

- The whole catalog is indexed (`--index all`, the default), so reads map
  competitively against every gene, as in Cayman: a read from an unrelated
  gene lands on its own gene instead of being forced onto a panel family.
  `--index hits` indexes only genes with an HMM hit (sub-cutoff homologs and
  decoys included), for a much smaller index; reads from other genes then
  fall back, or settle on a hit gene they are close to. `--index neighbours`
  adds those close genes: every gene without a hit whose 150 bp windows
  (every 50 bp) align to a hit gene at >= 90% over >= 45 bp. A read from one
  of them maps to its own gene and is dropped. On a 1 M-gene slice of the
  GMGC human-gut catalog this removed every false catalog call of the hits
  index (65 per 4 M reads where the strain is catalogued) at 1.3x its size;
  the whole-catalog index is 19x (osmo_refdb `studies/catalog_index`). This
  is the index osmo_refdb releases ship.
- A read aligned at >= 97% nucleotide identity over >= 45 bp (Cayman's
  defaults; tuning options `--catalog_min_identity`, `--catalog_min_aligned`) is settled by
  its gene: counted for the gene's family wherever it lies on the gene (on a
  fused gene where two families pass, for the family whose region holds the
  read), or dropped if the gene has no hit (background), is a homolog below
  the cutoff or an HMM decoy. A gene counts as the family when its
  full-sequence HMM score reaches the gathering cutoff, exactly as in
  `annotate --method hmm`.
- A read with several valid placements (alternative hits in bwa's `XA` tag
  that also pass the thresholds) on different families is ambiguous: it is
  settled, split equally over those families in `catalog_ambiguous`, and kept
  out of `raw_count`.
- Every other read goes on to DIAMOND, the cascade and the report levels
  (`--catalog_only` skips that step).
- `gene_counts.tsv` and `unit_counts.tsv` gain `catalog_count`,
  `catalog_ambiguous` and `translated_count`; `raw_count` is
  `catalog_count + translated_count`. A unit on `none` becomes
  `status catalog_only`: its catalog counts are reported (full-length calls,
  but only for strains in the catalog), `raw_count` stays NA.
- The catalog is tied to the release whose HMMs annotated it (checked by
  checksum) and to the aligner that indexed it (`catalog.json`). A shipped
  catalog always matches its release, and `aln_stats.tsv` records which
  catalog, release and HMM checksum a run used. bwa-mem3 is
  the default (bioconda `bwa-mem3`; needs AVX2 on x86-64, or ARM NEON);
  `--aligner bwa-mem2` also works, and `--aligner_bin` points at
  another executable with the same CLI and index format.
- Catalog-bound: above 97% identity only strains whose genes are in the
  catalog are seen. That suits well-catalogued communities such as the human
  gut; for environmental and halophile samples most reads fall back to
  translated search, and `profile` warns when the catalog settled less than
  half of the assigned fragments (`catalog_settled_fraction`).
- Accuracy is that of the full-length annotation (`annotate --method hmm`),
  including its errors: paralogs the HMMs confuse at full length (e.g. some
  opu permeases) and genes just below a gathering cutoff are not fixed by the
  catalog.

### `gene_coordinates.tsv` format (`annotate` mode only)

Unlike `gene_counts.tsv`, which sums hits per family, this lists one row per
individual called gene — its contig, start/end (from Prodigal, 1-based
inclusive), strand, and Prodigal's own per-contig gene index:

```
gene    protein_id  contig  start   end     strand  gene_index
ectA    NODE_1_7    NODE_1  1801    2400    1       7
ectB    NODE_1_8    NODE_1  2410    3200    1       8
```

Rows are sorted by gene, then contig, then start position. Requires
assembly-derived (Prodigal-style) coordinates, so it's only written in
`annotate` mode, and only when those coordinates were actually parsed — it's
skipped under the same conditions as `systems.tsv` (`--skip_systems`, or
`--proteins` supplied with non-Prodigal headers). Excluded families (per
`--exclude_families`) are left out here too, same as in `gene_counts.tsv`.

## NaCl tolerance prediction (`annotate` mode)

osmo_refdb v10 and later ship a model in `nacl_model/` that predicts the
highest NaCl concentration (% w/v) at which a **salt-out** prokaryote grows,
as BacDive records it, from a genome's annotate output. It was trained in
osmo_refdb (`pipeline/nacl_model/`) on 6,264 BacDive strains; the release's
`nacl_model/README.md` describes how it works, how well, and what it weighs.

- **Input.** The HMM counts of the 43 osmo families plus `total_proteins`, `gc`
  and `n_contigs`, all from `gene_counts.tsv` (osmotool >= 0.8.0 writes the
  header lines it needs). No taxonomy. Only `--method hmm` output against the
  same release is valid: the counts depend on that release's HMM cutoffs.
- **Model.** XGBoost regression trees. Each prediction is split exactly into a
  base value (the training mean) plus one contribution per feature (TreeSHAP),
  in % NaCl.
- **Accuracy.** Tested on genomes from genera the model never saw: Spearman
  0.69, mean absolute error 2.5% NaCl, R² 0.49. A model given only genome size
  and GC reaches R² 0.06, so the signal comes from the osmo genes.
- **How to read it.** Predictions are compressed toward the middle: strains that
  barely tolerate salt are over-predicted (about 3% for a true ≤1%), extreme
  halotolerants under-predicted (about 13% for a true >15%). Treat the value as
  a ranking plus a rough level, with the interval as its uncertainty.
- **Salt-in organisms** (Halobacteria, Salinibacter, Halanaerobiales)
  accumulate KCl instead of compatible solutes and are outside the model's
  domain; a separate detector flags genomes that look like them.

`<prefix>.nacl_prediction.tsv` (`predict-nacl -o` writes the same table with one
row per genome):

```
# model     nacl_max (osmo_refdb v10, annotate --method hmm)
# ...
genome  nacl_max_pred  nacl_max_lo  nacl_max_hi  p_salt_in  flags  base_value  top_contributions
ecoli   7.17           2.90         11.11        0.000      ok     5.917       trkH=2:+0.85;betL=0:-0.53;trkA=1:+0.27;...
```

| Column | Meaning |
|---|---|
| `nacl_max_pred` | predicted maximum NaCl, % w/v |
| `nacl_max_lo`, `nacl_max_hi` | 80% prediction interval (split-conformal, from out-of-genus errors) |
| `p_salt_in` | salt-in detector probability |
| `flags` | `salt_in_like` (prediction not reliable); `out_of_range(n)` (n features outside the training 1st–99th percentile); `many_contigs` (more contigs than 99% of training genomes); `few_proteins`; `no_systems` (no `systems.tsv`, so the salt-in detector is less reliable); `missing:<feature>`; `ok` if none |
| `base_value` | the starting point every prediction is built from (training mean) |
| `top_contributions` | the six largest contributions, `feature=value:+shap` |

`<prefix>.nacl_contributions.tsv` lists every feature's value and contribution,
largest first; `base_value` plus the `shap` column sums to the prediction
(before clipping at 0).

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

**Known refdb-level limitation (`opuB`)**: `opuBA` is the least reliable
family in the panel, so treat `opuBA` calls, and any `opuB` system call that
depends on one, with caution. In the original six-species validation (run on
releases before v9), `opuB` never resolved to a full colocalized operon, even
with the rescue above: it always had 1-2 of its 3 mandatory subunits. The
rescue fixes the `opuBC` side (a real labeling ambiguity). `opuBA` was a
different problem: in every species checked there was no hit at all, under
any family, where it should be, a recall gap in osmo_refdb's `opuBA` model
that positional reasoning can't fix. v9 does not solve this and adds a
precision problem on top. Its `pepABC_decoy` now takes peptide ABC ATPase
hits away from `opuBA`, but the model is still weak (from osmo_refdb's v9
benchmark):

- HMM recall of 0.13 on held-out reads, with a GA cutoff calibrated on only
  23 positives (F1 0.46)
- about 4100 false positives per million random-background reads
- genome-level presence doesn't track BacDive's maximum tolerated NaCl
  (rho -0.02, rank 40 of 45 families)

The remaining `opuBA` calls are probably other ABC ATPases that lack a
peptide C-terminal domain (bare ABC_tran), which a decoy can't separate. A
v9 `opuBA` hit therefore doesn't confirm the gene is there, and a missing
hit doesn't rule it out.

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
| `--total_fragments` | *(counted)* | `profile` only. Skip counting the input: read **pairs** for paired-end (R1's read count), reads for single-end |
| `--total_reads` | *(counted)* | `profile` only. Alternative to `--total_fragments`: reads summed over all input files (R1 + R2), halved for paired-end |
| `--method` | `hmm` (annotate only) | `diamond`, `hmm`, or `both`. `profile` is DIAMOND-only (HMM's GA cutoffs don't apply to short read fragments); `annotate` defaults to HMM since it's meaningfully more specific on full-length ORFs. Requires `hmms/osmo_refdb.hmm` present in `DATABASE` for `hmm`/`both` |
| `--no_cascade` | off (cascade on by default) | `profile` only. Disables the DIAMOND+HMM cascade, which by default gives reads whose DIAMOND call lands on a family with demonstrated precision problems (several as low as 5-25% DIAMOND-only precision) a second opinion from `hmmsearch`, using `DATABASE`'s `osmo_refdb.profile_cascade.tsv` + `hmms/osmo_refdb.hmm`. The cascade requires `orfm`+`hmmsearch` on PATH and both files present in `DATABASE`; pass this flag to run DIAMOND-only if you don't have them |
| `--habitat` | the one catalog `DATABASE` ships | `profile` only. Which of `DATABASE`'s `catalogs/<name>/` to use; needed only when it ships several — see [Gene catalog](#gene-catalog-profile-mode) |
| `--catalog_only` | off | `profile` only. Count catalog-settled reads only; no translated search for the rest |
| `--catalog_min_identity` / `--catalog_min_aligned` | 0.97 / 45 | `profile` only (tuning). Nucleotide identity and aligned bp for a read to be settled by the catalog |
| `--translated_only` | off | `profile` only (benchmarking). Skip the catalog; translated search for every read |
| `--no_markers` | off (on whenever `DATABASE` has `markers/`) | `profile` only. Skips the marker-gene DIAMOND pass and the `copies_per_genome` column / `markers.tsv` (see [Copies per genome](#copies-per-genome-profile-mode)) |
| `--proteins` | none | `annotate` only. Pre-called protein FASTA (skip Prodigal) |
| `--keep_aln` / `--keep_proteins` | off | Retain intermediate alignment/protein files instead of deleting them after the run |
| `--skip_systems` | off | `annotate` only. Skip the [multi-subunit systems](#multi-subunit-systems-annotate-mode-only) layer entirely — no `<prefix>.systems.tsv`, no `<prefix>.gene_coordinates.tsv`, no ORF-coordinate parsing |
| `--systems_config` | none | `annotate` only. TSV overriding `max_distance_genes`/`max_distance_bp` per system (columns: `system`, `max_distance_genes`, `max_distance_bp`) — see [Multi-subunit systems](#multi-subunit-systems-annotate-mode-only) |
| `--predict_nacl` | off | `annotate` only. Also predict the maximum NaCl tolerance with DATABASE's `nacl_model/` (osmo_refdb v10+; needs `--method hmm`/`both` and xgboost) — see [NaCl tolerance prediction](#nacl-tolerance-prediction-annotate-mode) |

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
    /path/to/osmo_refdb/releases/v10 \
    -1 ${SAMPLE}_R1.fastq.gz \
    -2 ${SAMPLE}_R2.fastq.gz \
    --out_prefix results/${SAMPLE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --tmpdir ${TMPDIR} \
    --total_fragments ${NPAIRS}   # R1's read count, from MultiQC if available
```

`$TMPDIR` is set by SLURM to fast node-local scratch — merged paired-end
FASTQ lands there instead of on the shared filesystem.

`--total_fragments` accepts a pre-computed count (e.g. from FastQC/MultiQC)
to skip counting the reads of very large files: for paired-end, the read
count of R1 alone (= pairs), not R1 + R2. `--total_reads` takes R1 + R2 instead.

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

## License and citation

osmotool is released under the MIT license (see [LICENSE](LICENSE)). The
osmo_refdb releases it downloads carry their own terms; see the release's
README. [CITATIONS.md](CITATIONS.md) lists the tools and methods osmotool
builds on; cite the ones your run used, together with the DOI of the
osmo_refdb release you used.
