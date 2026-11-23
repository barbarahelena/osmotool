# osmotool

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21787565.svg)](https://doi.org/10.5281/zenodo.21787565)

Detect bacterial osmoadaptation genes (compatible-solute synthesis and uptake, ion transporters, mechanosensitive channels) in metagenomes and in assembled genomes. The gene panel comes from an [osmo_refdb](https://github.com/barbarahelena/osmo_refdb) release, which osmotool downloads for you.

## What it does

**`profile`: metagenomic reads.** Counts each gene family in FASTQ reads, paired- or single-end. Reads are first matched against a gene catalog that the release has annotated at full length (human gut for v10). Reads the catalog can't place go to translated search with DIAMOND, and an HMM double-checks the families DIAMOND confuses. Counts come as reads per million and as copies per genome, the share of cells carrying the gene. Close paralogs that short reads can't separate are reported together, or not at all. → [docs/profile.md](docs/profile.md)

**`annotate`: genomes and MAGs.** Calls genes with Prodigal and assigns them with the release's HMMs, or with DIAMOND (`--method diamond`). It reports copies per family and checks whether multi-subunit systems, such as the OpuA transporter or the Mrp antiporter, sit together on a contig. With `--predict_nacl` (or `predict-nacl` afterwards) it predicts the genome's maximum NaCl tolerance. → [docs/annotate.md](docs/annotate.md)

| Subcommand | Input | Needs on PATH |
|---|---|---|
| `download-db` | | |
| `profile` | FASTQ reads | `bwa-mem3`, `diamond`, `orfm`, `hmmsearch` |
| `annotate` | Assembly FASTA (or proteins) | `prodigal`, `hmmsearch` (or `diamond`) |
| `predict-nacl` | Earlier `annotate` output | |

`build-catalog` and `index-catalog` prepare a gene catalog of your own for `profile`; see [docs/profile.md](docs/profile.md#your-own-catalog).

## Installation

### Container (recommended)

The [Docker image](https://hub.docker.com/r/barbarahelena/osmotool) has osmotool and every tool above installed, so nothing else needs installing. Its entry point is `osmotool`, so the subcommand follows the image name:

```bash
# mount the working directory so the container sees your files
docker run --rm -v "$PWD":/data -w /data barbarahelena/osmotool:0.10.1 \
  profile refdb/v10 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --out_prefix results/sample
```

Every osmotool release is published as `:X.Y.Z` and `:latest`. Pin the version, so you know which one you ran. The image is built for x86-64 (`linux/amd64`).

### Conda

```bash
conda env create -f environment.yml
conda activate osmotool
pip install osmotool          # or, from a clone: pip install -e .
```

With pip alone, install the tools in the table yourself, and add `pip install 'osmotool[predict]'` for NaCl prediction (xgboost).

## Reference database

```bash
osmotool download-db --release v10 --location /path/to/refdb
```

The release unpacks into `<location>/v10/`. Pass that directory to any subcommand as `DATABASE`; osmotool finds every file in it by name. `download-db` checks the archive's md5 against Zenodo before unpacking.

The current release is **osmo_refdb v10** ([10.5281/zenodo.24442164](https://doi.org/10.5281/zenodo.24442164)), with the human-gut gene catalog `gmgc_gut`, the report levels for `profile` and the NaCl tolerance model. It needs osmotool 0.10.0 or later. The release's README describes its contents and limitations.

To download by hand:

```bash
curl -LO https://zenodo.org/records/24442164/files/osmo_refdb_v10.tar.gz
curl -LO https://zenodo.org/records/24442164/files/osmo_refdb_v10.tar.gz.sha256
shasum -a 256 -c osmo_refdb_v10.tar.gz.sha256
tar -xzf osmo_refdb_v10.tar.gz      # -> osmo_refdb_v10/
```

## Quick start

```bash
osmotool download-db --release v10 --location refdb

# metagenome
osmotool profile refdb/v10 -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
  --out_prefix results/sample --threads 8

# genome, with NaCl tolerance prediction
osmotool annotate refdb/v10 genome.fasta --predict_nacl \
  --out_prefix results/genome --threads 8
```

With the container, put `docker run --rm -v "$PWD":/data -w /data barbarahelena/osmotool:0.10.1` where these commands say `osmotool`.

`profile` writes `<prefix>.unit_counts.tsv` (read results from this one), `gene_counts.tsv`, `aln_stats.tsv` and `markers.tsv`. `annotate` writes `gene_counts.tsv`, `systems.tsv` and `gene_coordinates.tsv`, plus `nacl_prediction.tsv` with `--predict_nacl`. The mode docs describe every column.

## License and citation

osmotool is released under the MIT license (see [LICENSE](LICENSE)). The osmo_refdb releases carry their own terms (v10 is CC BY 4.0); see the release's README.

If you use osmotool, cite the DOI of the version you ran. Every release is archived on Zenodo with its own DOI, for example v0.10.1: [10.5281/zenodo.24442151](https://doi.org/10.5281/zenodo.24442151). The DOI [10.5281/zenodo.21787565](https://doi.org/10.5281/zenodo.21787565) covers all versions and always resolves to the latest one; `osmotool --version` tells you which version you ran. Also cite the DOI of the osmo_refdb release you used. [CITATIONS.md](CITATIONS.md) lists the tools and methods osmotool builds on: cite the ones your run used.
