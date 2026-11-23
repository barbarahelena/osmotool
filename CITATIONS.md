# Citations

osmotool runs the tools below and follows the methods listed after them. If you publish results produced with osmotool, please cite:

- osmotool itself, by the Zenodo DOI of the version you ran. For v0.10.1 that is [10.5281/zenodo.24442151](https://doi.org/10.5281/zenodo.24442151); [10.5281/zenodo.21787565](https://doi.org/10.5281/zenodo.21787565) lists every version.
- the osmo_refdb release you ran it with, by its DOI (see the README's [Reference database](README.md#reference-database) section).
- the tools your run used, below.

## Tools

### Translated search and HMMs

- **DIAMOND** (`profile` translated search, `annotate --method diamond`)

  > Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods*. 2021;18:366–368. doi:[10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)

- **HMMER** (`annotate --method hmm`, the `profile` cascade, `build-catalog`)

  > Eddy SR. Accelerated profile HMM searches. *PLoS Comput Biol*. 2011;7(10):e1002195. doi:[10.1371/journal.pcbi.1002195](https://doi.org/10.1371/journal.pcbi.1002195)

- **OrfM** (open reading frames of reads, for the `profile` cascade)

  > Woodcroft BJ, Boyd JA, Tyson GW. OrfM: a fast open reading frame predictor for metagenomic data. *Bioinformatics*. 2016;32(17):2702–2703. doi:[10.1093/bioinformatics/btw241](https://doi.org/10.1093/bioinformatics/btw241)

- **Prodigal** (gene calling in `annotate`)

  > Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics*. 2010;11:119. doi:[10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119)

### Read mapping (`profile` catalog step, `build-catalog`, `index-catalog`)

- **bwa-mem3** (default aligner). A derivative of bwa-mem2 maintained by Fulcrum Genomics; no paper, cite the repository together with bwa-mem2.

  > <https://github.com/fg-labs/bwa-mem3>

- **BWA-MEM2** (bwa-mem3's basis, and the `--aligner bwa-mem2` option)

  > Vasimuddin M, Misra S, Li H, Aluru S. Efficient architecture-aware acceleration of BWA-MEM for multicore systems. In: *2019 IEEE International Parallel and Distributed Processing Symposium (IPDPS)*. 2019:314–324. doi:[10.1109/IPDPS.2019.00041](https://doi.org/10.1109/IPDPS.2019.00041)

- **BWA-MEM** (the algorithm both implement)

  > Li H. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*. 2013;1303.3997. [arXiv:1303.3997](https://arxiv.org/abs/1303.3997)

### NaCl tolerance prediction (`annotate --predict_nacl`, `predict-nacl`)

- **XGBoost** (the model)

  > Chen T, Guestrin C. XGBoost: a scalable tree boosting system. In: *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*. 2016:785–794. doi:[10.1145/2939672.2939785](https://doi.org/10.1145/2939672.2939785)

- **TreeSHAP** (per-feature contributions, computed by XGBoost's `pred_contribs`)

  > Lundberg SM, Erion G, Chen H, et al. From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell*. 2020;2:56–67. doi:[10.1038/s42256-019-0138-9](https://doi.org/10.1038/s42256-019-0138-9)

### Python libraries

- **pandas**

  > McKinney W. Data structures for statistical computing in Python. In: *Proceedings of the 9th Python in Science Conference*. 2010:56–61. doi:[10.25080/Majora-92bf1922-00a](https://doi.org/10.25080/Majora-92bf1922-00a)

- **NumPy**

  > Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. *Nature*. 2020;585:357–362. doi:[10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)

## Methods osmotool follows

- **MicrobeCensus**: copies-per-genome normalisation against universal single-copy genes, with MAD outlier removal

  > Nayfach S, Pollard KS. Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. *Genome Biol*. 2015;16:51. doi:[10.1186/s13059-015-0611-7](https://doi.org/10.1186/s13059-015-0611-7)

- **fetchMG / mOTUs marker genes**: the universal single-copy marker COGs behind `copies_per_genome`

  > Sunagawa S, Mende DR, Zeller G, et al. Metagenomic species profiling using universal phylogenetic marker genes. *Nat Methods*. 2013;10:1196–1199. doi:[10.1038/nmeth.2693](https://doi.org/10.1038/nmeth.2693)

## Reference data

osmotool's reference data come from an osmo_refdb release (CC BY 4.0): cite that release's Zenodo DOI. Its README lists the sources to credit: UniProt for the reference sequences, the Global Microbial Gene Catalog for the genes in `catalogs/`, and the prevalence-filtered GMGC catalog the `gmgc_gut` gene set is taken from ([Zenodo 10473258](https://zenodo.org/records/10473258)):

> Coelho LP, Alves R, del Río ÁR, et al. Towards the biogeography of prokaryotic genes. *Nature*. 2022;601:252–256. doi:[10.1038/s41586-021-04233-4](https://doi.org/10.1038/s41586-021-04233-4)

> Ducarmon QR, Karcher N, Giri S, et al. Cayman enables large-scale analysis of gut microbiome carbohydrate-active enzyme repertoires. *Nat Microbiol*. 2026;11:1739–1753. doi:[10.1038/s41564-026-02318-2](https://doi.org/10.1038/s41564-026-02318-2)

> The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2025. *Nucleic Acids Res*. 2025;53(D1):D609–D617. doi:[10.1093/nar/gkae1010](https://doi.org/10.1093/nar/gkae1010)
