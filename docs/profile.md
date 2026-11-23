# `osmotool profile`: genes in metagenomic reads

`profile` counts osmoadaptation gene families in short reads (FASTQ, paired- or single-end) and reports each as reads per million (RPM) and, when the release ships marker genes, as copies per genome.

It needs `bwa-mem3`, `diamond`, `orfm` and `hmmsearch` on PATH. The conda environment and the container provide all four.

- [Running it](#running-it)
- [How reads are assigned](#how-reads-are-assigned)
- [Output](#output)
- [Options](#options)
- [Your own catalog](#your-own-catalog)
- [Limitations](#limitations)

## Running it

```bash
# paired-end
osmotool profile osmo_refdb_v10 \
  -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --out_prefix results/sample --threads 8

# single-end
osmotool profile osmo_refdb_v10 \
  --singles sample.fastq.gz \
  --out_prefix results/sample --threads 8
```

The catalog's index (3.9 GB for v10) is loaded into memory. On a cluster:

```bash
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00

osmotool profile /path/to/osmo_refdb_v10 \
    -1 ${SAMPLE}_R1.fastq.gz -2 ${SAMPLE}_R2.fastq.gz \
    --out_prefix results/${SAMPLE} \
    --threads ${SLURM_CPUS_PER_TASK} \
    --tmpdir ${TMPDIR} \
    --total_fragments ${NPAIRS}   # optional: R1's read count, e.g. from MultiQC
```

`--tmpdir` puts temporary files (merged paired-end FASTQ) on node-local scratch instead of the shared filesystem. `--total_fragments` skips counting the input, which saves time on large files.

## How reads are assigned

Translated search on 150 bp reads has to tell close paralogs apart from a fragment, and that is where it fails. `profile` therefore works in two steps.

### Step 1: the gene catalog

The release ships a gene catalog under `catalogs/<name>/`. v10 ships `gmgc_gut`: the GMGC v1.0 human-gut genes with prevalence > 0.5%, annotated once at full length with the release's HMMs.

- A read aligned at ≥ 97% nucleotide identity over ≥ 45 bp is settled by the gene it maps to.
- If that gene belongs to a panel family, the read counts for that family. A gene belongs to a family when its full-length HMM score reaches the family's gathering cutoff, exactly as in `annotate --method hmm`.
- If the gene is background, a homolog below the cutoff, or a decoy, the read is settled but dropped. The index includes these competitors, so their reads land on their own gene instead of being forced onto a family.
- A read with valid placements on different families is ambiguous: it is split equally over them in `catalog_ambiguous` and kept out of `raw_count`.
- Every other read goes on to step 2.

### Step 2: translated search

Reads the catalog can't settle are searched with `diamond blastx` against the release's references (`--min_identity`, `--min_query_cover`, `--evalue`). For families where DIAMOND alone is imprecise (listed in `osmo_refdb.profile_cascade.tsv`), the read's ORF (from `orfm`) gets a second opinion from the family's HMM. This cascade is on by default; `--no_cascade` turns it off.

### Report units

Some families are close paralogs that a short read can't tell apart (e.g. the OpuA/OpuB/OpuC permeases), and some are not reliable on short reads at all. The release decides, on its calibration reads, at which level each can be reported, in `osmo_refdb.profile_report_levels.tsv`. v10 has:

- `gene`: 26 families reported on their own;
- `group`: 1 paralog group reported as the sum of its members (opuAB + opuBB + opuCB);
- `none`: 15 families not reported. Either they fall below F1 0.5, or they are annotate-only by scope (otsA, murB, mrpC, trkH, ktrB, ktrD).

The results are in `unit_counts.tsv`. These levels were decided on translated search. Catalog-settled reads are full-length calls, so a `none` unit still gets its catalog counts, with status `catalog_only`.

### Copies per genome

RPM still depends on what else is in the sample: host DNA, eukaryotes, viruses and the community's average genome size all shift it. With the marker genes the release ships in `markers/`, `profile` adds `copies_per_genome`:

> family reads per kb of effective gene length ÷ mean reads per kb of ten universal single-copy marker genes (ychF, pheS, argS, serS, cysS, leuS, valS, tsaD, ffh, ftsY)

This is the average number of copies per genome in the community. For a family with at most one copy per genome, it is the fraction of cells carrying it. The method follows MicrobeCensus:

- Markers are searched in a separate DIAMOND pass with the same thresholds, so they never compete with the panel for best hits.
- Gene length is an effective length, `3·L_aa − (2c − 1)·R + 1` (R = median read length, c = `--min_query_cover`): the number of read positions that can produce a qualifying hit. Without it, short families such as cspA would be underestimated against the longer markers.
- Undetected markers, and markers ≥ 1 scaled MAD from the median, are left out of the mean. With fewer than 5 markers left, the value is `NA`.

Detection sensitivity per family is not corrected for. Compare `copies_per_genome` across samples rather than reading it as an exact copy number. Run on host-depleted reads: host mitochondrial aminoacyl-tRNA synthetases look bacterial and can add marker hits. `--no_markers` skips the marker pass.

## Output

| File | Contents |
|---|---|
| `<prefix>.gene_counts.tsv` | Counts per family |
| `<prefix>.unit_counts.tsv` | Counts per report unit; read results from this table |
| `<prefix>.aln_stats.tsv` | Run summary, one `metric value` row each |
| `<prefix>.markers.tsv` | Per-marker evidence behind `copies_per_genome` |
| `<prefix>.catalog_hits.tsv` | With `--keep_aln`: every catalog-settled read, its gene, call, identity and aligned length |

### Counting in fragments

A fragment is a read pair (paired-end) or a read (single-end). Each mate of a pair counts 0.5, and `total_fragments` is the number of pairs, so a paired-end library and the same reads run single-end give the same RPM.

**RPM** = `raw_count × 10⁶ / total_fragments`. There is no gene-length correction: reads are not assembled.

### `gene_counts.tsv`

```
# total_fragments    1000000
# filtered_fragments 42753.5
gene    raw_count   RPM     catalog_count   catalog_ambiguous   translated_count   copies_per_genome
betA    12.5        12.5    10.0            0.50                2.5                0.412
ectA     4.0         4.0     0.0            0.00                4.0                0.087
```

- `raw_count` = `catalog_count + translated_count`.
- `catalog_ambiguous` is kept out of `raw_count`.
- With `--translated_only` the three source columns are left out.
- `copies_per_genome` appears when the release ships `markers/`.

This table lists every family on its own, including families whose gene-level count is not reliable. Use `unit_counts.tsv` for results.

### `unit_counts.tsv`

```
# total_fragments    1000000
# filtered_fragments 42753.5
# report_levels      osmo_refdb_v10/osmo_refdb.profile_report_levels.tsv
unit          report_level  members            status     raw_count  RPM
ectA          gene          ectA               reported   12.5       12.5
opu_permease  group         opuAB,opuBB,opuCB  reported   4.0        4.0
opuAC         none          opuAC              excluded   NA         NA
```

The catalog source columns and `copies_per_genome` follow as in `gene_counts.tsv`.

- `group` sums its members' counts. Each read counts for one family only, so this is the number of reads called as any member.
- `excluded` units get `NA`, never 0. Don't read them as "gene absent": the table says nothing about them.
- `catalog_only`: a `none` unit with catalog counts. These are reported, and `raw_count` stays `NA`.
- `not_in_manifest`: a label with reads but in no unit (e.g. the fused reference `mrpA_mrpB_fused`), listed so no read disappears silently.

Every unit gets a row, zeros included.

### `aln_stats.tsv`

Two columns, `metric` and `value`. Counts are in fragments, except the `catalog_reads_*` rows, where each mate counts 1.

| Metric | Meaning |
|---|---|
| `total_fragments` | The RPM denominator |
| `filtered_fragments` | Fragments with a catalog call or any DIAMOND alignment |
| `fragments_assigned`, `pct_fragments_assigned` | Fragments counted for a reported family, and their share of `filtered_fragments` |
| `min_identity`, `min_query_cover`, `min_subject_cover`, `evalue` | The DIAMOND filters used |
| `genome_coverage_rpk`, `markers_used`, ... | Marker normalisation summary |
| `report_levels`, `units_reported`, `units_excluded` | The report-levels file and its outcome |
| `catalog`, `catalog_release`, `catalog_hmm_sha256` | Which catalog the run used, and the release and HMM checksum it was annotated with (`none (--translated_only)` without one) |
| `catalog_dir`, `catalog_aligner`, `catalog_index` | Where it was read from, the aligner, and the index kind |
| `catalog_reads_called`, `catalog_reads_ambiguous` | Reads the catalog settled on one family, or on several |
| `catalog_reads_dropped`, `catalog_dropped_by_reason` | Reads settled on a gene that isn't counted (background, homolog, decoy) |
| `catalog_reads_fallback` | Reads that went on to translated search |
| `catalog_settled_fraction` | Share of assigned fragments the catalog settled (see below) |
| `translated_fallback` | `yes`, or `no` with `--catalog_only` |
| `note_systems` | Reminder that multi-subunit systems can't be checked from reads |
| `count_<family>` | `raw_count` per family |

`catalog_settled_fraction` is Σ(`catalog_count` + `catalog_ambiguous`) / Σ(those + `translated_count`) over the reported families. It is counted over panel assignments, not all reads: the index holds only panel genes and their close neighbours, so most reads map nowhere in it, and `catalog_reads_fallback` is large even when the catalog fits. Below 0.5, `profile` warns on stderr. The usual cause is a sample from a habitat the catalog doesn't cover. The counts are still valid, but they are then mostly translated-search counts. The 0.5 threshold is a starting point, not a calibrated value.

### `markers.tsv`

Each marker's count, reference length, reads per kb, and whether it was used, excluded as an outlier, or undetected.

## Options

The options most runs touch. `osmotool profile --help` lists them all.

| Option | Default | |
|---|---|---|
| `-1`/`-2`, `--singles` | | Paired-end or single-end reads, optionally gzipped |
| `--out_prefix` | required | Prefix for the output files, e.g. `results/sample` |
| `--threads` | 4 | Threads |
| `--tmpdir` | `$TMPDIR`, `$SCRATCH`, `/tmp` | Temporary files |
| `--total_fragments` | counted | Read pairs (R1's read count), or reads for single-end; skips counting. `--total_reads` takes R1 + R2 instead |
| `--habitat` | the release's only catalog | Which of `catalogs/<name>/` to use, once a release ships several |
| `--catalog_only` | off | Count catalog-settled reads only; skip translated search. Only for samples the catalog covers well |
| `--translated_only` | off | Skip the catalog, for benchmarking; also runs on releases without one |
| `--no_cascade` | off | DIAMOND only in step 2, e.g. without `orfm`/`hmmsearch` |
| `--no_markers` | off | Skip copies-per-genome |
| `--min_identity`, `--min_query_cover`, `--evalue` | 0.80, 0.80, 1e-5 | DIAMOND filters for step 2 |
| `--keep_aln` | off | Keep the DIAMOND alignments and write `catalog_hits.tsv` |

The catalog thresholds (`--catalog_min_identity`, `--catalog_min_aligned`) can be changed, but the defaults are what the release was validated with.

## Your own catalog

To profile against another gene set (e.g. UHGG genes, or a catalog for another habitat), annotate it once with the release's HMMs and pass its directory with `--catalog_dir`:

```bash
osmotool build-catalog osmo_refdb_v10 --genes my_genes.fna.gz -o my_catalog_v10 -t 32
osmotool profile osmo_refdb_v10 -1 R1.fq.gz -2 R2.fq.gz \
  --catalog_dir my_catalog_v10 --out_prefix results/sample
```

A shipped catalog is versioned with its release. A catalog of your own is yours to keep in step: `profile` refuses it if the release's HMMs have changed since it was built (checked by checksum in `catalog.json`).

`build-catalog --index` sets which genes are indexed:

- `all` (default): every gene, so reads map competitively.
- `hits`: only genes with an HMM hit (sub-cutoff homologs and decoys included). Much smaller, but reads from other genes can only fall back or land on a hit gene they resemble.
- `neighbours`: the hit genes plus every gene that 150 bp windows of a hit gene align to at ≥ 90% over ≥ 45 bp. Reads from those genes map to their own gene and are dropped. On a 1 M-gene slice of GMGC human gut, this removed every false catalog call of `hits` at 1.3× its size, against 19× for `all`. This is the index osmo_refdb releases ship.

A full index built once with `osmotool index-catalog` can be reused across releases with `profile --catalog_index`.

## Limitations

- **The catalog is human gut.** Above 97% identity it only sees strains whose genes are in GMGC. For environmental or halophile samples most reads fall back to translated search, and `profile` warns when the catalog settled less than half of the assigned fragments.
- **Catalog calls carry the full-length annotation's errors.** Paralogs the HMMs confuse at full length (e.g. some opu permeases) and genes just below a gathering cutoff are not fixed by the catalog.
- **The report levels describe translated search.** A unit at `none` may still be measurable through the catalog for catalogued strains.
- **No systems from reads.** Reads carry no contig positions, so multi-subunit systems can't be checked for co-location. That needs `annotate` on assemblies of the same samples.
