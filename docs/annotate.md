# `osmotool annotate`: genes in a genome

`annotate` calls osmoadaptation gene families in an assembled genome or MAG. It reports their copy numbers, checks whether multi-subunit systems (e.g. the OpuA transporter or the Mrp antiporter) sit together on a contig, and can predict the genome's maximum NaCl tolerance.

It needs `prodigal` and `hmmsearch` on PATH, or `diamond` with `--method diamond`. NaCl prediction also needs xgboost.

- [Running it](#running-it)
- [How genes are called](#how-genes-are-called)
- [Report units](#report-units)
- [Output](#output)
- [Multi-subunit systems](#multi-subunit-systems)
- [NaCl tolerance prediction](#nacl-tolerance-prediction)
- [Options](#options)

## Running it

```bash
# HMM (default), with Prodigal gene calling
osmotool annotate osmo_refdb_v10 assembly.fasta --out_prefix results/assembly

# with proteins you already called (skips Prodigal)
osmotool annotate osmo_refdb_v10 assembly.fasta \
  --proteins assembly_proteins.faa --out_prefix results/assembly

# also predict the maximum NaCl tolerance
osmotool annotate osmo_refdb_v10 assembly.fasta --predict_nacl \
  --out_prefix results/assembly

# predict NaCl tolerance later, from existing annotate output
osmotool predict-nacl osmo_refdb_v10 results/genomeA results/genomeB \
  -o results/nacl_predictions.tsv
```

## How genes are called

1. Prodigal (`-p meta`) calls the proteins, unless `--proteins` is given.
2. Each protein is searched against the release's HMMs with `hmmsearch --cut_ga`, using each family's calibrated gathering cutoff. Full-length proteins are the scale those cutoffs were calibrated on.
3. Each protein counts for its best-scoring family. Decoy families compete for proteins but are never reported (`osmo_refdb.annotate_excluded_families.txt`).

`--method diamond` uses `diamond blastp` instead, with the release's per-family bitscore cutoffs (`osmo_refdb.diamond_cutoffs.tsv`). It avoids the HMMER dependency, but HMM is more specific on full-length proteins. `--method both` runs the two side by side and writes `<prefix>.hmm.*` and `<prefix>.diamond.*`.

## Report units

Some families are close paralogs that the HMMs confuse even at full length (e.g. the OpuA/OpuB/OpuC permeases), and some are not reliable at all. The release decides, on its calibration proteins, how each can be reported, in `osmo_refdb.annotate_report_levels.tsv`:

- `gene`: the family is reported on its own;
- `group`: a paralog group is reported as the sum of its members;
- `none`: not reported. The `reason` column says why: below the F1 threshold, too few calibration sequences, or not stable under the bootstrap.

The levels were decided for `--method hmm`, so only HMM output gets `<prefix>.unit_counts.tsv` (with `--method both`, `<prefix>.hmm.unit_counts.tsv`). Read results from this table. `gene_counts.tsv` still lists every family on its own, and the NaCl model uses it as it is.

## Output

| File | Contents |
|---|---|
| `<prefix>.gene_counts.tsv` | Copies per family, with assembly statistics in the header |
| `<prefix>.unit_counts.tsv` | With `--method hmm`: copies per report unit; read results from this table |
| `<prefix>.systems.tsv` | Multi-subunit system calls |
| `<prefix>.gene_coordinates.tsv` | One row per called gene: contig, start, end, strand |
| `<prefix>.aln_stats.tsv` | Run summary |
| `<prefix>.nacl_prediction.tsv` | With `--predict_nacl`: predicted NaCl tolerance |
| `<prefix>.nacl_contributions.tsv` | With `--predict_nacl`: every feature's contribution |

### `gene_counts.tsv`

```
# total_proteins 4128
# genome_bp      4215606
# gc             0.43514
# n_contigs      1
# contig_n50     4215606
# filtered_reads 187
# proteins_with_hits 187
gene    raw_count   COPIES_PER_KB
ectA    1.0         5.87
ectB    1.0         6.42
```

- `raw_count` is the family's copy number in the genome.
- `COPIES_PER_KB` = `raw_count × 10⁹ / (gene_length_bp × total_proteins)`, a copy number corrected for gene length and proteome size, for comparing genomes of different sizes. It is not an abundance and can't be compared with `profile`'s RPM.
- `genome_bp` counts every sequence character, `N`s included. `gc` is (G+C) / (A+C+G+T), leaving ambiguity codes out. `n_contigs` and `contig_n50` describe the FASTA records. These four are computed from the assembly even when `--proteins` is given, so genome size and GC can be used as covariates without a separate step.
- `proteins_with_hits` counts proteins with a hit that passed filtering, before best-hit selection and decoy removal, so it can be slightly higher than the sum of `raw_count`. `filtered_reads` holds the same number under its older name; prefer `proteins_with_hits` in new code.

Parsers should skip `#` lines they don't recognise.

### `unit_counts.tsv`

```
# total_proteins      4128
# proteins_with_hits  187
# report_levels       osmo_refdb_v11/osmo_refdb.annotate_report_levels.tsv
unit          report_level  members            status    raw_count  COPIES_PER_KB
ectA          gene          ectA               reported  1.0        5.87
opu_permease  group         opuAB,opuBB,opuCB  reported  1.0        4.12
opuAC         none          opuAC              excluded  NA         NA
```

- `group` sums its members' copies. Each protein counts for one family only, so this is the number of proteins called as any member.
- `excluded` units get `NA`, never 0. Don't read them as "gene absent": the table says nothing about them.
- `not_in_manifest`: a label with a call but in no unit, listed so no protein disappears silently.

Every unit gets a row, zeros included. `aln_stats.tsv` gains `report_levels`, `units_reported` and `units_excluded`.

### `gene_coordinates.tsv`

```
gene    protein_id  contig  start   end     strand  gene_index
ectA    NODE_1_7    NODE_1  1801    2400    1       7
ectB    NODE_1_8    NODE_1  2410    3200    1       8
```

Coordinates come from Prodigal (1-based, inclusive). Rows are sorted by gene, contig and start. The file is skipped, like `systems.tsv`, with `--skip_systems` or when `--proteins` has no Prodigal-style headers.

## Multi-subunit systems

Several families are subunits of one complex. Counting each subunit, as `gene_counts.tsv` does, can't tell a complete operon from subunits scattered over the genome. `systems.tsv` adds that check:

| System | Mandatory subunits | Accessory | Fusion |
|---|---|---|---|
| `opuA` | opuAA, opuAB, opuAC | | |
| `opuB_opuC` | opuBA *or* opuCA, opuBB *or* opuCB, opuBC *or* opuCC | | |
| `mrp` | mrpA, mrpB, mrpC, mrpE | mrpD, mrpF, mrpG | mrpA+mrpB as one ORF (`mrpA_mrpB_fused`) |
| `ect` | ectA, ectB, ectC | | |
| `ots` | otsA, otsB | | |
| `trk` | trkA, trkH | | |
| `ktr` | ktrA, ktrB *or* ktrD | | |
| `galE_murB_mazG` | galE, murB, mazG | | |

| Column | Meaning |
|---|---|
| `present` | Every mandatory subunit was found together on one contig |
| `subunit_count` | Subunits found anywhere in the assembly, one per mandatory slot (ktrB and ktrD together count once) plus accessory subunits |
| `colocalized` | Same as `present`, kept as its own column |
| `fragmented_assembly_flag` | Subunits were found, but on different contigs, so co-location can't be judged. Subunits too far apart on the same contig are a real "not colocalized" |
| `contig` | The contig the system is on (`-` otherwise) |
| `mandatory_hits`, `accessory_hits` | The subunits found |
| `rescued_hits` | Subunits recovered from a paralog's label, as `family<-original_label` |

### When subunits count as together

The mandatory subunits, sorted by position, must form one chain:

- **Distance.** Every consecutive pair is within 3 intervening genes and 2000 bp. In complete genomes of six species (*B. subtilis*, *S. aureus*, *L. monocytogenes*, *E. faecalis*, *H. elongata*, *E. coli*), the `mrp`, `opuA`, OpuB/OpuC-type, `ect` and `ots` operons all had 0 genes between subunits. MacSyFinder's CRISPR-Cas models use the same 3-gene limit. `--systems_config` overrides the window per system.
- **Strand.** Consecutive subunits share a strand, with no gene on the opposite strand between them. Several species carry a paralogous copy of an operon right next to the first, on the opposite strand (e.g. *H. elongata*'s second `mrp`-like cluster, 236 bp away), and this check keeps the two apart.
- **`galE_murB_mazG` is a neighbourhood, not an operon.** Culligan et al. (2012, *ISME J* 6:1916–1925, [doi:10.1038/ismej.2012.38](https://doi.org/10.1038/ismej.2012.38)) found the three genes on one ~39 kb fosmid from a *Collinsella*-related organism, each conferring salt tolerance in *E. coli*. In *C. aerofaciens* ATCC 25986, mazG sits 23 genes (20.8 kb) from galE/murB, across a strand flip. This system uses a 30-gene / 35 kb window and no strand check. Of the validation genomes, it fires only on *C. aerofaciens*.

### Fusions and paralog rescue

- **mrpA+mrpB fusion.** A hit on the `mrpA_mrpB_fused` DIAMOND reference satisfies both subunits. There is no fused HMM, so `--method hmm` undercounts the lineages that carry the fusion; most don't.
- **Paralog rescue.** When a mandatory subunit has no hit under its own label, a hit under a known paralog that sits in the right place, next to another subunit of the same system, is moved to the missing slot and listed in `rescued_hits` (e.g. `trkH<-ktrB`). It applies to the Trk/Ktr pore proteins (trkH, ktrB, ktrD). The opu substrate-binding proteins are not rescued: swaps between OpuB and OpuC labels are covered by `opuB_opuC`'s either/or slots, and a rescue between opuA and `opuB_opuC` would rest on ATPase and permease labels that are just as ambiguous. proX is never rescued either: it belongs to ProU, a different system, and rescuing from it called *E. coli*'s ProVWX operon an `opuA`.

### opuB_opuC

OpuB and OpuC are paralogous operons from a gene duplication (ATPase, two permeases, substrate-binding protein), and the HMMs can't tell their subunits apart. osmo_refdb v11's report levels don't report opuBA, opuBC, opuCA or opuCC at any level (calibration F1 0.05–0.35), and report opuBB and opuCB only within the `opu_permease` group. Called as two systems, the operon was found but its identity followed single label flips. *E. faecalis* V583's one operon, labelled opuBB, opuCC, opuCB, opuBA, came out as `opuB` with v11 and as `opuC` with v9, where its ATPase was labelled opuCA. *B. subtilis* 168 carries both operons side by side and was called `opuC` only.

`opuB_opuC` takes either label in each slot, the way `ktr` takes ktrB or ktrD. It is present in *B. subtilis*, *L. monocytogenes*, *S. aureus* and *E. faecalis*, and absent in *H. elongata*, *E. coli* K-12, *C. aerofaciens*, *H. volcanii* and *R. baltica*. In *E. coli* the yehWXYZ (osmF) transporter gives an opuBA and an opuBB call but no binding protein, so the call stays at 2 subunits.

It tells you an OpuB/OpuC-type transporter is there, not which one, and not how many copies: the substrate range sits in the binding protein, the worst-separated label (opuCC F1 0.07). `mandatory_hits` lists the labels of the operon that satisfied the system. Don't read them as OpuB versus OpuC. `opuA` is still called on its own.

## NaCl tolerance prediction

The release ships a model in `nacl_model/` that predicts the highest NaCl concentration (% w/v) at which a **salt-out** prokaryote grows, as BacDive records it. The release's `nacl_model/README.md` describes it in full.

- **Input.** The copy numbers of the 44 families from `--method hmm`, plus `total_proteins`, `gc` and `n_contigs`, all from `gene_counts.tsv`. No taxonomy. Only `--method hmm` output against the same release is valid: the counts depend on that release's cutoffs.
- **Model.** XGBoost regression trees, trained on 6,262 BacDive strains from 1,570 genera. Each prediction splits exactly into the training mean plus one contribution per feature (TreeSHAP), in % NaCl.
- **Accuracy.** On genera the model never saw: Spearman 0.68, mean absolute error 2.5% NaCl, R² 0.49. Genome size and GC alone reach R² 0.05, so the signal comes from the osmo genes.
- **How to read it.** Predictions are pulled toward the middle: strains with a true maximum below 1% are predicted at about 3%, extreme halotolerants (true mean 22%) at about 13%. Treat the value as a ranking and a rough level, with the interval as its uncertainty.
- **Salt-in organisms** (Halobacteria, Salinibacter, Halanaerobiales) use KCl instead of compatible solutes and are outside the model's domain. A separate detector flags genomes that look like them.

`<prefix>.nacl_prediction.tsv` (`predict-nacl -o` writes one row per genome):

```
# model     nacl_max (osmo_refdb v10, annotate --method hmm)
genome  nacl_max_pred  nacl_max_lo  nacl_max_hi  p_salt_in  flags  base_value  top_contributions
ecoli   7.17           2.90         11.11        0.000      ok     5.917       trkH=2:+0.85;betL=0:-0.53;trkA=1:+0.27;...
```

| Column | Meaning |
|---|---|
| `nacl_max_pred` | Predicted maximum NaCl, % w/v |
| `nacl_max_lo`, `nacl_max_hi` | 80% prediction interval |
| `p_salt_in` | Salt-in detector probability |
| `flags` | `ok`, or: `salt_in_like` (prediction not reliable); `out_of_range(n)` (n features outside the training 1st–99th percentile); `many_contigs`; `few_proteins`; `no_systems` (salt-in detector less reliable); `missing:<feature>` |
| `base_value` | The training mean every prediction starts from |
| `top_contributions` | The six largest contributions, `feature=value:+shap` |

`<prefix>.nacl_contributions.tsv` lists every feature's value and contribution, largest first. `base_value` plus the contributions gives the prediction (before clipping at 0).

`predict-nacl` takes the `--out_prefix` of earlier `annotate --method hmm` runs. For `--method both` output, pass `PREFIX.hmm`.

## Options

The options most runs touch. `osmotool annotate --help` lists them all.

| Option | Default | |
|---|---|---|
| `--out_prefix` | required | Prefix for the output files |
| `--method` | `hmm` | `hmm`, `diamond`, or `both` |
| `--proteins` | | Proteins you called yourself; skips Prodigal |
| `--predict_nacl` | off | Also predict the maximum NaCl tolerance (needs xgboost) |
| `--skip_systems` | off | No `systems.tsv` or `gene_coordinates.tsv` |
| `--systems_config` | | TSV (`system`, `max_distance_genes`, `max_distance_bp`) overriding the window per system |
| `--threads` | 4 | Threads |
| `--min_identity`, `--min_query_cover`, `--evalue` | 0.80, 0.80, 1e-5 | DIAMOND filters, for `--method diamond`/`both` |
| `--keep_proteins`, `--keep_aln` | off | Keep the Prodigal proteins and the alignments |
