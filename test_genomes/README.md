# test_genomes/

Real, complete reference genomes used to validate `osmotool annotate`'s
complex/system-aware detection layer (`src/osmotool/systems.py`) against
actual biology, not just synthetic unit-test fixtures. All are NCBI
RefSeq assemblies, downloaded via `https://ftp.ncbi.nlm.nih.gov/genomes/all/...`
(accessions verified through NCBI E-utilities, not guessed from memory).

## Genomes

| File | Organism | Assembly | Sequence(s) | Size |
|---|---|---|---|---|
| `ecoli_k12.fna` | *Escherichia coli* str. K-12 substr. MG1655 | GCF_000005845.2 | NC_000913.3 | 4.5 MB |
| `bsubtilis_168.fna` | *Bacillus subtilis* subsp. subtilis str. 168 | GCF_000009045.1 (ASM904v1) | NC_000964.3 | 4.1 MB |
| `saureus.fna` | *Staphylococcus aureus* subsp. aureus NCTC 8325 | GCF_000013425.1 (ASM1342v1) | NC_007795.1 | 2.7 MB |
| `lmonocytogenes.fna` | *Listeria monocytogenes* EGD-e | GCF_000196035.1 (ASM19603v1) | NC_003210.1 | 2.8 MB |
| `efaecalis.fna` | *Enterococcus faecalis* V583 | GCF_000007785.1 (ASM778v1) | NC_004668.1 + 3 more (plasmids) | 3.2 MB |
| `helongata.fna` | *Halomonas elongata* DSM 2581 | GCF_000196875.2 (ASM19687v2) | NC_014532.2 | 3.9 MB |
| `collinsella.fna` | *Collinsella aerofaciens* ATCC 25986 (strain JCM 10188) | GCF_010509075.1 (ASM1050907v1) | NZ_CP048433.1 | 2.4 MB |
| `hvolcanii.fna` | *Haloferax volcanii* DS2 (archaeon) | GCF_000025685.1 | NC_013967.1 + 4 plasmids | 3.9 MB |
| `rbaltica.fna` | *Rhodopirellula baltica* SH 1 (Planctomycetota) | — (NC_005027, no versioned GCF recorded) | NC_005027.1 | 6.9 MB |

**Total: ~34 MB.** `ecoli_k12`/`hvolcanii`/`rbaltica` predate this validation
work — they're osmo_refdb's own real-genome validation set, documented in
`pipeline/real_genome_validation/README.md` in that repo. H. volcanii was
specifically chosen there for "zero training-set overlap" (an archaeon
sharing no lineage with anything osmo_refdb's HMMs/DIAMOND references were
built from) — the same role it plays here: a genuinely independent negative
control, not just "some other domain of life." The other six were pulled
specifically for this validation, each chosen as the standing NCBI RefSeq
"reference genome" (or "Complete Genome" where no "reference genome" flag
exists) for its organism.

## Why these particular species

- **B. subtilis, S. aureus, L. monocytogenes, E. faecalis** — Firmicutes
  reference organisms, chosen to check `opuA`/`opuB`/`opuC`/`mrp`/`trk`/`ktr`
  colocalization across multiple independent lineages, not just one. B.
  subtilis specifically is the organism most of these systems were originally
  characterized in.
- **H. elongata** — the actual ectoine-biosynthesis discovery organism;
  the only genome in this set that gave a real `ect` (ectA/B/C) positive.
- **C. aerofaciens** — the organism (or closest relative) the
  `galE_murB_mazG` locus finding (Culligan et al. 2012, *ISME J*,
  doi:10.1038/ismej.2012.38) is from.
- **E. coli** — cross-domain negative control for the Firmicutes-specific
  systems (`mrp`/`opuA`/`opuB`/`opuC`/`ect`), and the real positive for
  `ots` (otsBA) and the real negative-but-present-elsewhere case for `trk`
  (trkA/trkH genuinely both present, genuinely not colocalized).
- **H. volcanii** — an archaeon with, per osmo_refdb's own real-genome
  validation notes, zero training-set overlap with anything its HMMs/DIAMOND
  references were built from: a genuinely independent negative control, not
  just "a different domain of life." Still produced a plausible `mrp` true
  positive here (Mrp/Mnh-type antiporters aren't Firmicutes-exclusive).
- **R. baltica** — Planctomycetota cross-phylum control, also from
  osmo_refdb's own validation set.

## What was found

Full reasoning and evidence live in `src/osmotool/systems.py`'s module
docstring and `README.md`'s "Multi-subunit systems" section (the
calibration numbers, the paralog-confusion rescue mechanism, and the
`galE_murB_mazG` system are all cited against specific measurements from
these genomes) — this is a pointer, not a duplicate:

- Every real operon found (`mrp`, `opuA`, `opuC` across multiple species;
  `ect` and `galE_murB_mazG` in their respective source organisms; `otsBA`
  in E. coli) sits within a few genes / under 100bp of itself, informing the
  3-gene/2000bp default window.
- Several species carry a duplicate/paralogous operon copy nearby on the
  *opposite* strand — motivated the strand-consistency check.
- `opuB` never fully resolves in any species tested; traced to two distinct
  causes (a real refdb HMM recall gap in `opuBA`, and a genuine sequence
  mislabeling of `opuBC` as `opuCC`) rather than one.
- `galE_murB_mazG` fires only on *C. aerofaciens*, confirming the wide
  (30-gene/35000bp), strand-agnostic window is specific, not just permissive.
- A paralog-confusion rescue mechanism initially included `proX` as a
  rescue-eligible sibling for `opuAC`/`opuBC`/`opuCC`; this produced a false
  positive on E. coli (`opuA` "present" via E. coli's own real, but
  non-orthologous, ProU operon) and was corrected — `proX` is now excluded.

## Reproducing

Each genome can be re-fetched directly from its accession, e.g.:

```bash
curl -sL "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1_genomic.fna.gz" \
  | gunzip > bsubtilis_168.fna
```

Substitute the accession path from the table above (NCBI's FTP layout is
`GCF/<3>/<3>/<3>/<accession>_<assembly_name>/<accession>_<assembly_name>_genomic.fna.gz`,
split from the accession's own digits).
