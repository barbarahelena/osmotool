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
docstring in runners.py and this project's README for why (no contig-level
coordinates at typical short-read lengths, so proximity/colocalization simply
isn't checkable there).
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
        of which satisfies that slot. Almost always a single-member
        frozenset (e.g. {"opuAA"}) -- a multi-member slot models a real
        either/or paralog relationship (e.g. ktr's pore-forming subunit is
        satisfied by either ktrB or ktrD, confirmed as paralogs of each
        other in families.yaml) without inventing separate machinery for it.
        A system is only `present` when EVERY slot is satisfied by a hit
        co-located with the others (see detect_systems).
    accessory:
        Families that raise confidence / subunit_count but are never
        required for `present`.
    max_distance_genes / max_distance_bp:
        Proximity window for colocalization -- BOTH must hold (genes as the
        primary, length-robust unit; bp as a sanity bound so a handful of
        unusually short/long intervening ORFs can't silently blow the
        window open). These are the extension spec's suggested defaults
        (8 genes / 8000bp), shipped as configurable-but-UNCALIBRATED
        starting points -- see this project's README for the calibration
        caveat and --systems_config for overriding them per system without
        touching code.
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
    max_distance_genes: int = 8
    max_distance_bp: int = 8000
    fusion_label: str | None = None
    fusion_slots: frozenset[str] = frozenset()


DEFAULT_SYSTEMS: dict[str, SystemDef] = {
    "opuA": SystemDef(
        name="opuA",
        mandatory=(frozenset({"opuAA"}), frozenset({"opuAB"}), frozenset({"opuAC"})),
    ),
    "opuB": SystemDef(
        name="opuB",
        mandatory=(frozenset({"opuBA"}), frozenset({"opuBB"}), frozenset({"opuBC"})),
    ),
    "opuC": SystemDef(
        name="opuC",
        mandatory=(frozenset({"opuCA"}), frozenset({"opuCB"}), frozenset({"opuCC"})),
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


def _slot_satisfied_on_contig(
    slot: frozenset[str], hits_by_family: dict[str, list[_Hit]], contig: str
) -> list[_Hit]:
    """Hits on *contig* for any family in *slot*."""
    out = []
    for family in slot:
        out.extend(h for h in hits_by_family.get(family, ()) if h.coord.contig == contig)
    return out


def _within_window(a: _Hit, b: _Hit, sys_def: SystemDef) -> bool:
    intervening_genes = abs(a.coord.gene_index - b.coord.gene_index) - 1
    bp_distance = max(a.coord.start, b.coord.start) - min(a.coord.end, b.coord.end)
    return (
        intervening_genes <= sys_def.max_distance_genes
        and bp_distance <= sys_def.max_distance_bp
    )


def _contig_is_colocalized(
    contig: str, sys_def: SystemDef, hits_by_family: dict[str, list[_Hit]]
) -> tuple[bool, list[_Hit]]:
    """
    Return (colocalized, chosen_hits) for *contig*: whether every mandatory
    slot has >=1 hit on this contig with SOME combination (one hit per slot)
    all pairwise within the proximity window. Systems here have at most a
    handful of mandatory slots and typically 1-2 hits per family on a given
    contig, so exhaustively trying every combination is cheap and avoids the
    false negatives a single-candidate-per-slot greedy pick would produce.
    """
    per_slot_candidates = [
        _slot_satisfied_on_contig(slot, hits_by_family, contig) for slot in sys_def.mandatory
    ]
    if any(not candidates for candidates in per_slot_candidates):
        return False, []

    for combo in itertools.product(*per_slot_candidates):
        if all(
            _within_window(combo[i], combo[j], sys_def)
            for i in range(len(combo)) for j in range(i + 1, len(combo))
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
        {orf_id: OrfCoord} from functions.parse_orf_coordinates(). An empty
        dict (e.g. non-Prodigal --proteins headers) means no ORF has known
        coordinates, so every system comes back with subunit_count=0-only
        detection is impossible -- callers should skip writing systems
        output entirely in that case (see runners.run_annotate).
    """
    calls: dict[str, SystemCall] = {}
    for name, sys_def in systems.items():
        hits_by_family = _collect_hits(best_hits, orf_coords, sys_def)

        mandatory_families_found = sorted(
            fam for slot in sys_def.mandatory for fam in slot if fam in hits_by_family
        )
        accessory_families_found = sorted(fam for fam in sys_def.accessory if fam in hits_by_family)
        subunit_count = len(mandatory_families_found) + len(accessory_families_found)

        if subunit_count == 0:
            calls[name] = SystemCall(
                system=name, present=False, subunit_count=0, colocalized=False,
                fragmented_assembly_flag=False, contig=None,
                mandatory_hits=(), accessory_hits=(),
            )
            continue

        candidate_contigs = {
            h.coord.contig for hits in hits_by_family.values() for h in hits
        }

        colocalized_contig = None
        chosen_hits: list[_Hit] = []
        for contig in sorted(candidate_contigs):
            ok, chosen = _contig_is_colocalized(contig, sys_def, hits_by_family)
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
        )

    return calls
