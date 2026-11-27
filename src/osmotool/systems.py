"""
systems.py — complex/system-aware detection layer for `annotate` mode

Per-gene HMM/DIAMOND calls (functions.py, quantifier.py) score each gene
family independently. This module adds a proximity/quorum layer on top: for
multi-subunit systems (e.g. the OpuA transporter, the Mrp/Mnh Na+/H+
antiporter complex), a set of individually-called subunit genes only means
something biologically if they're actually co-located on the same contig,
close enough to plausibly be one operon.

This is deliberately bespoke rather than a MacSyFinder wrapper: MacSyFinder
expects to run its own hmmsearch stage against its own model format, so
adopting it would mean re-deriving its inputs from osmotool's own
already-filtered, already-calibrated hits (or running detection twice) --
not worth the new dependency for a proximity/quorum check this small.

`osmo_refdb` doesn't declare system membership itself (families.yaml
documents each family individually) -- DEFAULT_SYSTEMS below is osmotool's
own layer on top, grounded in what v7 actually built full subunit coverage
for (checked against families.yaml/CHANGELOG.md/hmms/, not assumed from the
extension spec's suggested table, which named a five-system list including a
`kdp` system that isn't buildable -- kdpB/kdpC were never added to the panel,
only kdpA exists -- and guessed the wrong Mrp fusion pair, mrpA/mrpD instead
of the actually-confirmed mrpA/mrpB).

`profile` (read) mode deliberately gets none of this -- see run_profile's
docstring in runners.py and docs/profile.md for why (no contig-level
coordinates at typical short-read lengths, so proximity/colocalization simply
isn't checkable there).

max_distance_genes/max_distance_bp defaults (see SystemDef) were checked
against real complete-genome operons rather than shipped as the extension
spec's untested guess, across six species: mrpA-G (B. subtilis, S. aureus,
L. monocytogenes, H. elongata), opuABC (B. subtilis, L. monocytogenes),
opuCABC (B. subtilis, S. aureus, L. monocytogenes, E. faecalis), ectABC
(H. elongata -- the actual ectoine-discovery organism), and otsBA (E. coli)
are ALL 0 intervening genes apart (gaps of -34 to +93bp -- essentially
back-to-back CDSs, typical of a real co-transcribed operon), independent of
species. STRING's precomputed "neighborhood" evidence score for B. subtilis
mrpA/mrpB corroborates this from an entirely separate source (0.803 out of
1.0, aggregated across STRING's own much larger reference-genome set).
MacSyFinder/TXSScan's own expert-curated models (github.com/macsy-models)
independently converge on a similar range for comparably-sized, similarly
tightly-packed operons -- CasFinder's CRISPR-Cas models (7-9 genes, like
`mrp`) use inter_gene_max_space=3; TXSScan's small transport systems like
T1SS/T2SS (3-4 genes, like `opuA`/`ots`) use 5. Only genuinely large,
multi-locus machinery (T6SS, the flagellum -- nothing in this panel is that
complex) goes up to 20.

The same multi-species sweep also surfaced a real false-bridging risk the
spec's original 8-gene/8000bp suggestion didn't guard against: several
species carry a duplicate/paralogous copy of the same operon nearby, always
on the OPPOSITE strand -- e.g. H. elongata's second mrp-like cluster sits
immediately adjacent to the first (0 genes / 236bp), well inside even the
original window. Distance alone can't reliably reject that; strand
information can, for free, with no threshold to calibrate (see
_strand_consistent) -- a real operon is one uninterrupted same-strand
transcriptional unit, and every genuine operon in this sweep was.

Tightened to 3 genes / 2000bp (matching CasFinder's precedent exactly): a
solid margin above every real operon gap observed across all six species,
and a solid margin below the false-bridging distances seen. Still not
exhaustive -- loosen via --systems_config if your own genomes show real
operons this misses.

OpuB and OpuC are one system here, `opuB_opuC`, not two. They are
paralogous operons from a gene duplication (ATPase, two permeases, a
substrate-binding protein that sets the substrate range), and the HMMs
can't tell their subunits apart: osmo_refdb v11's
annotate_report_levels.tsv reports none of opuBA/opuBC/opuCA/opuCC at any
level (calibration F1_all 0.05-0.35; opuCC 0.07) and opuBB/opuCB only as
the opu_permease group. As separate systems, the operon was found but its
identity followed single label flips. E. faecalis V583's one operon
(labelled opuBB, opuCC, opuCB, opuBA) was called opuB under v11, after
rescuing opuBC<-opuCC, and opuC under v9, where the ATPase was labelled
opuCA. B. subtilis 168 carries both operons side by side and got opuC only.
`opuB_opuC` asks only what the labels can answer -- is there an
OpuB/OpuC-type operon -- with an either/or slot per subunit, the way `ktr`
takes ktrB or ktrD. It doesn't say which paralogue it is, or how many
copies the genome has. `opuA` stays its own system: it resolved correctly
in every test genome.
"""

from __future__ import annotations

import csv
import itertools
import logging
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from osmotool.functions import OrfCoord
from osmotool.utils import gene_family_from_header

log = logging.getLogger("osmotool")


# ---------------------------------------------------------------------------
# System definitions
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SystemDef:
    """
    One multi-subunit system's definition.

    mandatory:
        Tuple of "slots"; each slot is a frozenset of family labels, ANY one
        of which satisfies that slot. Usually a single-member frozenset
        (e.g. {"opuAA"}) -- a multi-member slot models a real either/or
        paralog relationship without inventing separate machinery for it
        (ktr's pore-forming subunit is satisfied by either ktrB or ktrD,
        confirmed as paralogs of each other in families.yaml; every
        opuB_opuC slot takes either the OpuB or the OpuC label).
        A system is only `present` when EVERY slot is satisfied by a hit
        co-located with the others (see detect_systems).
    accessory:
        Families that raise confidence / subunit_count but are never
        required for `present`.
    max_distance_genes / max_distance_bp:
        Proximity window for colocalization -- BOTH must hold (genes as the
        primary, length-robust unit; bp as a sanity bound so a handful of
        unusually short/long intervening ORFs can't silently blow the
        window open). Set PER SYSTEM, not a single global default -- most
        systems here are real single-operon claims calibrated tight (3
        genes / 2000bp, see this module's docstring for the evidence), but
        `galE_murB_mazG` is a different KIND of claim (a wider genomic
        neighborhood, not one transcriptional unit -- see its own def
        below) and is calibrated on its own, much looser scale. See
        docs/annotate.md, and --systems_config for overriding a system's
        distances without touching code.
    require_strand_consistency:
        Whether _strand_consistent (see detect_systems) is enforced for
        this system. True (default) for real single-operon claims, where
        "one contiguous same-strand transcriptional unit" is exactly what's
        being tested. False for wider claims like `galE_murB_mazG`, where
        there's no reason to expect one strand across a ~20-40kb genomic
        region and enforcing it would incorrectly reject the one confirmed
        real instance (a strand flip sits between mazG and murB/galE in
        both the source clone and *C. aerofaciens*'s own reference genome).
    fusion_label / fusion_slots:
        Some lineages carry two adjacent mandatory subunits as a single
        fused ORF. When set, a hit whose family label is exactly
        `fusion_label` (a real, distinct reference family in osmo_refdb,
        not inferred) satisfies every slot named in `fusion_slots`
        simultaneously, from that one ORF's coordinates.
    """
    name: str
    mandatory: tuple[frozenset[str], ...]
    accessory: frozenset[str] = frozenset()
    max_distance_genes: int = 3
    max_distance_bp: int = 2000
    require_strand_consistency: bool = True
    fusion_label: str | None = None
    fusion_slots: frozenset[str] = frozenset()


DEFAULT_SYSTEMS: dict[str, SystemDef] = {
    "opuA": SystemDef(
        name="opuA",
        mandatory=(frozenset({"opuAA"}), frozenset({"opuAB"}), frozenset({"opuAC"})),
    ),
    # OpuB and OpuC together (see the module docstring): each slot takes
    # either paralogue's label. Named so it can't be mistaken for the opuBC
    # family.
    "opuB_opuC": SystemDef(
        name="opuB_opuC",
        mandatory=(
            frozenset({"opuBA", "opuCA"}),  # ATPase
            frozenset({"opuBB", "opuCB"}),  # permease
            frozenset({"opuBC", "opuCC"}),  # substrate-binding protein
        ),
    ),
    "mrp": SystemDef(
        name="mrp",
        mandatory=(
            frozenset({"mrpA"}), frozenset({"mrpB"}),
            frozenset({"mrpC"}), frozenset({"mrpE"}),
        ),
        accessory=frozenset({"mrpD", "mrpF", "mrpG"}),
        fusion_label="mrpA_mrpB_fused",
        fusion_slots=frozenset({"mrpA", "mrpB"}),
    ),
    "ect": SystemDef(
        name="ect",
        mandatory=(frozenset({"ectA"}), frozenset({"ectB"}), frozenset({"ectC"})),
    ),
    "ots": SystemDef(
        name="ots",
        mandatory=(frozenset({"otsA"}), frozenset({"otsB"})),
    ),
    "trk": SystemDef(
        name="trk",
        mandatory=(frozenset({"trkA"}), frozenset({"trkH"})),
    ),
    "ktr": SystemDef(
        name="ktr",
        mandatory=(frozenset({"ktrA"}), frozenset({"ktrB", "ktrD"})),
    ),
    # Not a transcriptional operon -- a wider genomic-neighborhood claim.
    # Culligan EP, Sleator RD, Marchesi JR, Hill C. "Functional metagenomics
    # reveals novel salt tolerance loci from the human gut microbiome."
    # ISME J. 2012;6(10):1916-1925. doi:10.1038/ismej.2012.38, PMID
    # 22534607. A transposon-mutagenesis screen of a human gut metagenomic
    # fosmid library found galE/murB/mazG (individually flagged elsewhere in
    # this panel as too promiscuous a domain to be specific on their own --
    # see quantifier.py / families.yaml) together on one ~39kb fosmid clone
    # (clone "SMG3", closest relative Collinsella aerofaciens; locus tags
    # COLAER_01955/01957/01981) -- each individually conferred NaCl/KCl
    # tolerance when expressed in an osmosensitive E. coli strain. galE and
    # murB sit close (226bp apart in the source clone; confirmed again at 1
    # gene / 1701bp in C. aerofaciens ATCC 25986's own reference genome,
    # NCBI GCF_010509075.1), but mazG is real genomic distance away -- 24
    # locus tags in the source clone, confirmed at 23 genes / 20810bp in
    # ATCC 25986, WITH a strand flip in between (hence
    # require_strand_consistency=False here specifically). Distances set
    # with margin above that confirmed real gap, bounded by the fosmid
    # insert size itself (~39kb is close to a hard ceiling for "found
    # together on one clone").
    "galE_murB_mazG": SystemDef(
        name="galE_murB_mazG",
        mandatory=(frozenset({"galE"}), frozenset({"murB"}), frozenset({"mazG"})),
        max_distance_genes=30,
        max_distance_bp=35000,
        require_strand_consistency=False,
    ),
}


# ---------------------------------------------------------------------------
# Paralog-confusion rescue
# ---------------------------------------------------------------------------
#
# Some osmo_refdb families are documented as genuinely hard to separate by
# sequence alone -- close paralogs sharing the same Pfam domain, explicitly
# flagged in families.yaml as cross-confused rather than just theoretically
# similar. When one of these families comes up empty for a system, but a
# hit under one of its documented siblings sits positionally right where
# that subunit should be (adjacent to the system's OTHER already-confirmed
# real subunits -- see _rescue_missing_mandatory_families), that's much
# stronger evidence of what the gene actually is than the ambiguous
# sequence-based label alone -- but it is still a POSITIONAL/sequence
# inference, not a functional or evolutionary-identity confirmation
# (the kind Culligan et al. actually did, by heterologous expression +
# phenotype assay -- osmotool has no equivalent and can't produce one).
# Deliberately narrow: only cliques families.yaml explicitly calls out as
# systematically swapped for each other, not every family sharing a
# broad/promiscuous domain (e.g. NOT opuAA/opuBA/opuCA or opuAB/opuBB/opuCB,
# which families.yaml gives only the generic ABC-transporter
# cross-exclusion).
#
# The opu substrate-binding clique (opuAC, opuBC, opuCC; PF04069) was
# dropped when opuB and opuC became one system:
#   - opuBC<->opuCC, the only rescue that ever fired on the test genomes,
#     is now an either/or slot of opuB_opuC, and a same-slot alternative is
#     never a rescue (see _rescue_missing_mandatory_families).
#   - What was left were A<->B/C rescues: moving a binding protein into
#     opuA from an opuBC/opuCC hit, or into opuB_opuC from an opuAC hit, on
#     the strength of ATPase/permease anchors that are themselves
#     A/B/C-ambiguous (osmo_refdb v11 reports opuAA/opuBA/opuCA at no level
#     and opuAB/opuBB/opuCB only as the opu_permease group). That could
#     pull a real OpuB/OpuC operon's binding protein into opuA. opuA never
#     needed it: on the nine test genomes, B. subtilis and L. monocytogenes
#     have opuA present from their own opuAC hits, and no opuA call used a
#     rescue.
# proX (ProU) was never rescue-eligible either. ProU is a separate,
# non-orthologous system that only shares the binding-protein domain, and
# rescuing opuAC from proX once called E. coli's ProVWX operon an opuA (see
# test_proX_never_rescues_opuAC).
_CONFUSION_CLIQUES = (
    # Pore-forming subunit, PF02386. trkH (Proteobacteria naming) and
    # ktrB/ktrD (Firmicutes naming) are documented in families.yaml as the
    # SAME underlying Trk/Ktr K+-uptake system under lineage-specific naming
    # conventions, not two independently-arisen analogous systems the way
    # ProU/Opu are -- but this is this module's own reading of that
    # documentation, not something functionally confirmed. Worth an
    # independent check by whoever has closer domain knowledge of
    # osmo_refdb's own family curation.
    frozenset({"trkH", "ktrB", "ktrD"}),
)

PARALOG_CONFUSION_CLIQUES: dict[str, frozenset[str]] = {
    fam: frozenset(clique - {fam})
    for clique in _CONFUSION_CLIQUES
    for fam in clique
}


def load_systems_config(
    path: str | Path, base: dict[str, SystemDef] = DEFAULT_SYSTEMS
) -> dict[str, SystemDef]:
    """
    Override max_distance_genes/max_distance_bp for one or more systems from
    a TSV (columns: system, max_distance_genes, max_distance_bp), without
    touching mandatory/accessory/fusion membership -- that's structural, not
    a calibration knob, so it stays code-defined in DEFAULT_SYSTEMS.

    A system named in the file but not in *base* is ignored (with a
    warning); a system in *base* but absent from the file keeps its
    default distances unchanged.
    """
    overrides: dict[str, tuple[int, int]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            overrides[row["system"]] = (
                int(row["max_distance_genes"]), int(row["max_distance_bp"])
            )

    systems: dict[str, SystemDef] = {}
    for name, sys_def in base.items():
        if name in overrides:
            genes, bp = overrides[name]
            systems[name] = SystemDef(
                name=sys_def.name,
                mandatory=sys_def.mandatory,
                accessory=sys_def.accessory,
                max_distance_genes=genes,
                max_distance_bp=bp,
                require_strand_consistency=sys_def.require_strand_consistency,
                fusion_label=sys_def.fusion_label,
                fusion_slots=sys_def.fusion_slots,
            )
        else:
            systems[name] = sys_def

    unknown = set(overrides) - set(base)
    if unknown:
        log.warning("systems_config names unknown system(s), ignored: %s", ", ".join(sorted(unknown)))

    return systems


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SystemCall:
    system: str
    present: bool
    subunit_count: int
    colocalized: bool
    fragmented_assembly_flag: bool
    contig: str | None
    mandatory_hits: tuple[str, ...]
    accessory_hits: tuple[str, ...]
    rescued_hits: tuple[str, ...] = ()


@dataclass
class _Hit:
    family: str
    orf_id: str
    coord: OrfCoord


def _collect_hits(
    best_hits: dict[str, dict], orf_coords: dict[str, OrfCoord], sys_def: SystemDef
) -> dict[str, list[_Hit]]:
    """
    {family: [_Hit, ...]} for every mandatory+accessory family of *sys_def*
    that has >=1 best-hit with known ORF coordinates. A fusion_label hit
    (e.g. "mrpA_mrpB_fused") is expanded into a _Hit under every family in
    fusion_slots, all sharing that one ORF's coordinates.
    """
    wanted = frozenset.union(*sys_def.mandatory) | sys_def.accessory
    hits: dict[str, list[_Hit]] = defaultdict(list)

    for orf_id, row in best_hits.items():
        coord = orf_coords.get(orf_id)
        if coord is None:
            continue
        family = gene_family_from_header(row["sseqid"])

        if sys_def.fusion_label is not None and family == sys_def.fusion_label:
            for fused_family in sys_def.fusion_slots:
                if fused_family in wanted:
                    hits[fused_family].append(_Hit(fused_family, orf_id, coord))
            continue

        if family in wanted:
            hits[family].append(_Hit(family, orf_id, coord))

    return dict(hits)


def _collect_rescue_hits(
    best_hits: dict[str, dict], orf_coords: dict[str, OrfCoord], sys_def: SystemDef
) -> dict[str, list[_Hit]]:
    """
    {mandatory_family: [_Hit, ...]} -- rescue candidates for each of
    *sys_def*'s mandatory families that has documented paralog-confusion
    siblings (PARALOG_CONFUSION_CLIQUES). Each _Hit keeps its own TRUE
    family label (e.g. "ktrB"), not the family it might rescue -- these
    are only accepted once independently confirmed positionally adjacent to
    a real subunit of the SAME system (see
    _rescue_missing_mandatory_families), never on label alone.
    """
    mandatory_families = frozenset(fam for slot in sys_def.mandatory for fam in slot)
    rescuable = {
        fam: PARALOG_CONFUSION_CLIQUES[fam]
        for fam in mandatory_families if fam in PARALOG_CONFUSION_CLIQUES
    }
    if not rescuable:
        return {}

    candidates: dict[str, list[_Hit]] = defaultdict(list)
    for orf_id, row in best_hits.items():
        coord = orf_coords.get(orf_id)
        if coord is None:
            continue
        family = gene_family_from_header(row["sseqid"])
        for target_fam, siblings in rescuable.items():
            if family in siblings:
                candidates[target_fam].append(_Hit(family, orf_id, coord))

    return dict(candidates)


def _rescue_missing_mandatory_families(
    hits_by_family: dict[str, list[_Hit]],
    rescue_candidates: dict[str, list[_Hit]],
    sys_def: SystemDef,
    strand_index: dict[str, dict[int, int]],
) -> tuple[dict[str, list[_Hit]], dict[str, str]]:
    """
    For each mandatory family with zero real hits, accept the first rescue
    candidate (see _collect_rescue_hits) confirmed positionally adjacent
    (same contig, within the proximity window AND strand-consistent -- see
    _within_window / _strand_consistent) to a REAL hit of some OTHER
    mandatory family of this same system. This is the only thing that makes
    a rescue legitimate: an isolated, far-away hit under a confusable
    sibling label has no more evidence of being this family than of being
    the label it already carries -- only adjacency to an already-confirmed
    subunit of the SAME operon does.

    A confirmed rescue is merged into the returned hits_by_family under its
    TARGET family (e.g. a mislabeled "ktrB" hit satisfying trk's trkH slot
    is stored as a trkH-family _Hit, not ktrB) so every downstream step
    (subunit_count, the colocalization search, mandatory_hits reporting)
    treats it exactly like a real hit. Returns (augmented_hits_by_family,
    {target_family: original_sibling_label}) -- the second dict is only
    for SystemCall.rescued_hits' transparency, not used for detection logic.

    Operates per SLOT, not per individual family: a multi-member OR-slot
    (e.g. ktr's {ktrB, ktrD}) is skipped entirely once ANY of its own
    members has a real hit, exactly like _slot_satisfied_on_contig's OR
    semantics -- this matters because ktrB/ktrD are simultaneously (a)
    interchangeable alternatives within ktr's own slot AND (b) each other's
    documented confusion siblings for the SEPARATE trk/ktr cross-system
    clique. Rescuing "ktrB" from a real "ktrD" hit whenever ktrD already
    satisfies the slot directly would be redundant at best and, worse,
    creates two candidate labels for the same physical ORF competing in the
    colocalization search below -- the OR-slot mechanism already covers
    that case with no rescue needed.
    """
    augmented = {fam: list(hs) for fam, hs in hits_by_family.items()}
    rescued_from: dict[str, str] = {}

    all_mandatory_members = [fam for slot in sys_def.mandatory for fam in slot]
    anchors = [h for fam in all_mandatory_members if fam in hits_by_family for h in hits_by_family[fam]]

    for slot in sys_def.mandatory:
        if any(fam in hits_by_family for fam in slot):
            continue  # slot already satisfied directly -- no rescue needed

        target = sorted(slot)[0]  # canonical attribution when this slot needs rescuing
        for fam in sorted(slot):
            for candidate in rescue_candidates.get(fam, ()):
                if candidate.family in slot:
                    continue  # a same-slot alternative isn't a real rescue -- see docstring
                confirmed = any(
                    candidate.coord.contig == anchor.coord.contig
                    and _within_window(candidate, anchor, sys_def)
                    and _strand_consistent(candidate, anchor, strand_index)
                    for anchor in anchors
                )
                if confirmed:
                    augmented.setdefault(target, []).append(_Hit(target, candidate.orf_id, candidate.coord))
                    rescued_from[target] = candidate.family
                    break
            if target in rescued_from:
                break

    return augmented, rescued_from


def _slot_satisfied_on_contig(
    slot: frozenset[str], hits_by_family: dict[str, list[_Hit]], contig: str
) -> list[_Hit]:
    """Hits on *contig* for any family in *slot*, in sorted family order so
    the combination _contig_is_colocalized picks (and so mandatory_hits)
    doesn't depend on string-hash order when a slot has several labels."""
    out = []
    for family in sorted(slot):
        out.extend(h for h in hits_by_family.get(family, ()) if h.coord.contig == contig)
    return out


def _build_strand_index(orf_coords: dict[str, OrfCoord]) -> dict[str, dict[int, int]]:
    """{contig: {gene_index: strand}} for EVERY predicted ORF (not just
    system-family hits) -- lets _strand_consistent check whether a gene of
    ANY kind sitting between two hits breaks the single-strand run a real
    operon should be, not just the strand of the two hits themselves."""
    index: dict[str, dict[int, int]] = defaultdict(dict)
    for coord in orf_coords.values():
        index[coord.contig][coord.gene_index] = coord.strand
    return dict(index)


def _within_window(a: _Hit, b: _Hit, sys_def: SystemDef) -> bool:
    intervening_genes = abs(a.coord.gene_index - b.coord.gene_index) - 1
    bp_distance = max(a.coord.start, b.coord.start) - min(a.coord.end, b.coord.end)
    return (
        intervening_genes <= sys_def.max_distance_genes
        and bp_distance <= sys_def.max_distance_bp
    )


def _strand_consistent(a: _Hit, b: _Hit, strand_index: dict[str, dict[int, int]]) -> bool:
    """
    A real operon is transcribed as one contiguous block on a single strand.
    Calibration-free complement to the distance window (see detect_systems'
    docstring): requires *a* and *b* to be on the same strand, AND no ORF of
    ANY kind strictly between them (by gene_index, same contig) to sit on
    the opposite strand -- a stronger check than "the two hits share a
    strand", which alone wouldn't catch e.g. two same-strand hits that are
    actually in different transcriptional units separated by an
    opposite-strand gene.

    Checked against real B. subtilis/E. coli data: every real operon tested
    (mrp, opuA, opuB/opuC, otsBA) is a single unbroken same-strand run. This is
    a complementary signal, not a substitute for the distance window --
    e.g. E. coli's trkA/trkH (hundreds of genes apart, not a real operon)
    also happen to share a strand, so distance is still what rejects it.
    """
    if a.coord.strand != b.coord.strand:
        return False
    lo, hi = sorted((a.coord.gene_index, b.coord.gene_index))
    contig_strands = strand_index.get(a.coord.contig, {})
    return all(
        contig_strands.get(idx, a.coord.strand) == a.coord.strand
        for idx in range(lo + 1, hi)
    )


def _contig_is_colocalized(
    contig: str,
    sys_def: SystemDef,
    hits_by_family: dict[str, list[_Hit]],
    strand_index: dict[str, dict[int, int]],
) -> tuple[bool, list[_Hit]]:
    """
    Return (colocalized, chosen_hits) for *contig*: whether every mandatory
    slot has >=1 hit on this contig with SOME combination (one hit per slot)
    forming one contiguous run -- i.e. every CONSECUTIVE pair, after sorting
    the combination by gene_index, is within the proximity window AND
    strand-consistent (see _within_window / _strand_consistent). This is
    deliberately a chain constraint, not "every pair of mandatory hits
    within the window" (a clique constraint): real calibration against
    Bacillus subtilis 168's mrp operon (mrpA-G) showed mrpA and mrpE --
    both mandatory, but with the non-mandatory mrpD sitting directly
    between them -- span more bp than any flat per-system budget should
    have to cover, even though every adjacent gene pair in the real operon
    is essentially back-to-back. A clique constraint would incorrectly
    reject that real operon; the chain constraint doesn't, while still
    catching genuinely scattered hits (each adjacent gap still has to be
    tight). Systems here have at most a handful of mandatory slots and
    typically 1-2 hits per family on a given contig, so exhaustively trying
    every combination is cheap and avoids the false negatives a single-
    candidate-per-slot greedy pick would produce.
    """
    per_slot_candidates = [
        _slot_satisfied_on_contig(slot, hits_by_family, contig) for slot in sys_def.mandatory
    ]
    if any(not candidates for candidates in per_slot_candidates):
        return False, []

    for combo in itertools.product(*per_slot_candidates):
        ordered = sorted(combo, key=lambda h: h.coord.gene_index)
        if all(
            _within_window(ordered[i], ordered[i + 1], sys_def)
            and (
                not sys_def.require_strand_consistency
                or _strand_consistent(ordered[i], ordered[i + 1], strand_index)
            )
            for i in range(len(ordered) - 1)
        ):
            return True, list(combo)
    return False, []


def detect_systems(
    best_hits: dict[str, dict],
    orf_coords: dict[str, OrfCoord],
    systems: dict[str, SystemDef] = DEFAULT_SYSTEMS,
) -> dict[str, SystemCall]:
    """
    Run every system in *systems* against one assembly's best-hit calls.

    best_hits:
        {orf_id: alignment_row} from quantifier.select_best_hits() -- the
        same per-gene calls already used for gene_counts.tsv, not a
        separate detection pass.
    orf_coords:
        {orf_id: OrfCoord} for EVERY predicted ORF (not just osmoadaptation
        family hits), from functions.parse_orf_coordinates() -- the strand-
        consistency check (_strand_consistent) needs the full set to tell
        whether some *other* gene sitting between two hits breaks a
        same-strand run. An empty dict (e.g. non-Prodigal --proteins
        headers) means no ORF has known coordinates, so every system comes
        back with subunit_count=0-only detection is impossible -- callers
        should skip writing systems output entirely in that case (see
        runners.run_annotate).
    """
    strand_index = _build_strand_index(orf_coords)
    calls: dict[str, SystemCall] = {}
    for name, sys_def in systems.items():
        hits_by_family = _collect_hits(best_hits, orf_coords, sys_def)

        rescue_candidates = _collect_rescue_hits(best_hits, orf_coords, sys_def)
        rescued_from: dict[str, str] = {}
        if rescue_candidates:
            hits_by_family, rescued_from = _rescue_missing_mandatory_families(
                hits_by_family, rescue_candidates, sys_def, strand_index
            )
        rescued_hits = tuple(f"{fam}<-{orig}" for fam, orig in sorted(rescued_from.items()))

        mandatory_families_found = sorted(
            fam for slot in sys_def.mandatory for fam in slot if fam in hits_by_family
        )
        accessory_families_found = sorted(fam for fam in sys_def.accessory if fam in hits_by_family)
        # one per satisfied slot, not per label: opuB_opuC with all six
        # labels found is still 3 subunits
        mandatory_slots_found = sum(
            1 for slot in sys_def.mandatory if any(fam in hits_by_family for fam in slot)
        )
        subunit_count = mandatory_slots_found + len(accessory_families_found)

        if subunit_count == 0:
            calls[name] = SystemCall(
                system=name, present=False, subunit_count=0, colocalized=False,
                fragmented_assembly_flag=False, contig=None,
                mandatory_hits=(), accessory_hits=(), rescued_hits=(),
            )
            continue

        candidate_contigs = {
            h.coord.contig for hits in hits_by_family.values() for h in hits
        }

        colocalized_contig = None
        chosen_hits: list[_Hit] = []
        for contig in sorted(candidate_contigs):
            ok, chosen = _contig_is_colocalized(contig, sys_def, hits_by_family, strand_index)
            if ok:
                colocalized_contig = contig
                chosen_hits = chosen
                break

        if colocalized_contig is not None:
            calls[name] = SystemCall(
                system=name, present=True, subunit_count=subunit_count, colocalized=True,
                fragmented_assembly_flag=False, contig=colocalized_contig,
                mandatory_hits=tuple(sorted({h.family for h in chosen_hits})),
                accessory_hits=tuple(accessory_families_found),
                rescued_hits=rescued_hits,
            )
            continue

        # Not colocalized on any single contig. Distinguish "mandatory hits
        # genuinely span multiple contigs" (fragmented assembly -- can't be
        # jointly evaluated at all) from "all on one contig but farther
        # apart than the window" (a real not-colocalized call).
        mandatory_contigs = {
            h.coord.contig
            for slot in sys_def.mandatory
            for fam in slot
            for h in hits_by_family.get(fam, ())
        }
        fragmented = len(mandatory_contigs) > 1

        calls[name] = SystemCall(
            system=name, present=False, subunit_count=subunit_count, colocalized=False,
            fragmented_assembly_flag=fragmented, contig=None,
            mandatory_hits=tuple(mandatory_families_found),
            accessory_hits=tuple(accessory_families_found),
            rescued_hits=rescued_hits,
        )

    return calls
