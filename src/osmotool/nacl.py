"""
nacl.py — predict a genome's maximum NaCl tolerance from annotate output

The model ships with an osmo_refdb release as ``<release>/nacl_model/``
(osmo_refdb pipeline/nacl_model/, packaged by 13_package_nacl_model.py):

  feature_schema.json   feature order and definitions, target transform,
                        prediction-interval table, applicability ranges,
                        salt-in detector coefficients, provenance
  nacl_max.xgb.json     XGBoost regression model
  README.md             what the model predicts, how well, what it weighs

It predicts the highest NaCl concentration (% w/v) at which a salt-out
prokaryote grows, as BacDive records it, from the files `annotate` writes
for one genome:

  <prefix>.gene_counts.tsv   per-family raw_count and the header lines
                             '# total_proteins', '# gc', '# n_contigs'
  <prefix>.systems.tsv       per-system calls (used by the salt-in detector)

The model was trained on `annotate --method hmm` output made with the same
release, so it is only applied to HMM counts from that release.

Feature construction and the outputs here must match osmo_refdb's
pipeline/nacl_model/features.py and predict.py exactly; the tests in
tests/test_nacl.py and a parity run on real genomes check that.

xgboost is an optional dependency (``pip install osmotool[predict]``),
imported only when a model is loaded.
"""

from __future__ import annotations

import csv
import json
import logging
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

log = logging.getLogger("osmotool")

SCHEMA_FILE = "feature_schema.json"
SELF_CHECK_FILE = "self_check.json"
# xgboost < 3.1 loads models saved by 3.x without error but misreads their
# base_score (osmo_refdb v10's model: 0.48 instead of 5.92), shifting every
# prediction; see SELF_CHECK_FILE for the guard that catches this generally
MIN_XGBOOST = (3, 1)
SUPPORTED_SCHEMA_VERSIONS = {1}
DECOY_SUFFIX = "_decoy"
# features the model cannot do without; a gene_counts.tsv from osmotool
# < 0.8.0 lacks the gc / n_contigs header lines
REQUIRED_HEADER = ("total_proteins", "gc", "n_contigs")
N_TOP = 6


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _version_tuple(v: str) -> tuple[int, ...]:
    out = []
    for part in v.split(".")[:2]:
        digits = "".join(ch for ch in part if ch.isdigit())
        out.append(int(digits) if digits else 0)
    return tuple(out)


def require_xgboost():
    """Import xgboost (>= MIN_XGBOOST) or raise a RuntimeError that says how to
    install it."""
    want = ".".join(map(str, MIN_XGBOOST))
    try:
        import xgboost
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError(
            "NaCl prediction needs the xgboost package, which is not installed. "
            f"Install it with `pip install 'osmotool[predict]'` (or `pip install "
            f"'xgboost>={want}'`, or `conda install -c conda-forge 'xgboost>={want}'`)."
        ) from exc
    if _version_tuple(xgboost.__version__) < MIN_XGBOOST:
        raise RuntimeError(
            f"NaCl prediction needs xgboost >= {want}, found {xgboost.__version__}: "
            "older versions load the model but misread it, giving wrong "
            f"predictions. Upgrade with `pip install -U 'xgboost>={want}'`."
        )
    return xgboost


@dataclass
class NaclModel:
    schema: dict[str, Any]
    booster: Any            # xgboost.Booster
    model_dir: Path

    @property
    def release(self) -> str:
        return self.schema["refdb_release"]

    @property
    def features(self) -> list[str]:
        return self.schema["features"]


def load_schema(model_dir: str | Path) -> dict[str, Any]:
    model_dir = Path(model_dir)
    schema = json.loads((model_dir / SCHEMA_FILE).read_text())
    version = schema.get("schema_version")
    if version not in SUPPORTED_SCHEMA_VERSIONS:
        raise ValueError(
            f"{model_dir / SCHEMA_FILE}: schema_version {version!r} is not supported by "
            f"this osmotool (supports {sorted(SUPPORTED_SCHEMA_VERSIONS)}); update osmotool"
        )
    method = schema.get("annotate", {}).get("method")
    if method != "hmm":
        raise ValueError(f"{model_dir / SCHEMA_FILE}: model expects annotate method "
                         f"{method!r}; this osmotool only supports 'hmm'")
    return schema


def load_nacl_model(model_dir: str | Path) -> NaclModel:
    """Load feature_schema.json and the XGBoost model it names."""
    model_dir = Path(model_dir)
    schema = load_schema(model_dir)
    xgb = require_xgboost()
    booster = xgb.Booster()
    booster.load_model(str(model_dir / schema["model_file"]))
    names = booster.feature_names
    if names is not None and list(names) != list(schema["features"]):
        raise ValueError(f"{model_dir}: the model's feature names do not match "
                         f"{SCHEMA_FILE}; repackage the release")
    self_check(model_dir, booster, xgb)
    return NaclModel(schema=schema, booster=booster, model_dir=model_dir)


def self_check(model_dir: Path, booster, xgb) -> None:
    """Recompute the packaged reference outputs (self_check.json, written by
    osmo_refdb's 13_package_nacl_model.py) and refuse a model this xgboost
    does not reproduce. Skipped when the release has no self_check.json."""
    path = model_dir / SELF_CHECK_FILE
    if not path.exists():
        log.debug("no %s in %s; model self-check skipped", SELF_CHECK_FILE, model_dir)
        return
    import numpy as np
    chk = json.loads(path.read_text())
    got = booster.predict(xgb.DMatrix(np.array(chk["rows"], dtype=float),
                                      feature_names=chk["features"]), output_margin=True)
    diff = float(np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(chk["output_margin"]))))
    if diff > chk.get("tolerance", 1e-4):
        raise RuntimeError(
            f"{model_dir}: this xgboost ({xgb.__version__}) does not reproduce the "
            f"model's reference outputs (packaged with xgboost {chk.get('xgboost_version')}; "
            f"max difference {diff:.3g}), so its predictions would be wrong. Install "
            f"xgboost {chk.get('xgboost_version')} or a compatible version.")


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------

def read_gene_counts(path: str | Path) -> tuple[dict[str, str], dict[str, float]]:
    """-> ({header_key: value}, {family: raw_count}); decoys are skipped."""
    header: dict[str, str] = {}
    counts: dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                parts = line[1:].strip().split("\t")
                if len(parts) >= 2:
                    header[parts[0].strip()] = parts[1].strip()
            elif not line or line.startswith("gene\t"):
                continue
            else:
                f = line.split("\t")
                if len(f) >= 2 and not f[0].endswith(DECOY_SUFFIX):
                    counts[f[0]] = float(f[1])
    return header, counts


def read_systems(path: str | Path | None) -> dict[str, tuple[int, int, int]]:
    """-> {system: (present 0/1, subunit_count, fragmented 0/1)}; {} if no file."""
    out: dict[str, tuple[int, int, int]] = {}
    if path is None or not Path(path).exists():
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["system"]] = (int(r["present"] == "True"), int(r["subunit_count"]),
                                int(r.get("fragmented_assembly_flag") == "True"))
    return out


def _num(x: str | None) -> float | None:
    if x is None or x in ("", "NA"):
        return None
    return float(x)


def features_from_annotate(
    gene_counts: str | Path,
    systems: str | Path | None,
    schema: dict[str, Any],
) -> dict[str, float | None]:
    """Raw feature dict for one genome, as osmo_refdb's features.py builds it.

    A family absent from gene_counts.tsv counts as 0. Every name the schema
    or the salt-in detector can ask for is present; base features missing
    from the header are None. 'fragmented_systems' (systems with
    fragmented_assembly_flag set) is informational, not a model feature or
    a flag: in the BacDive training set it is set in 75% of genomes (mostly
    galE_murB_mazG, whose genes are rarely adjacent) and is not linked to
    larger prediction errors.
    """
    header, counts = read_gene_counts(gene_counts)
    sysd = read_systems(systems)
    families, system_names = schema["families"], schema["systems"]
    bp = _num(header.get("genome_bp"))
    feat: dict[str, float | None] = {
        "total_proteins": _num(header.get("total_proteins")),
        "gc": _num(header.get("gc")),
        "n_contigs": _num(header.get("n_contigs")),
        "genome_mb": bp / 1e6 if bp is not None else None,
    }
    for f in families:
        feat[f] = counts.get(f, 0.0)
    feat["n_families_present"] = float(sum(1 for f in families if counts.get(f, 0) > 0))
    feat["total_osmo_copies"] = float(sum(counts.get(f, 0.0) for f in families))
    for s in system_names:
        p, n, _ = sysd.get(s, (0, 0, 0))
        feat[f"sys_{s}_present"] = float(p)
        feat[f"sys_{s}_subunits"] = float(n)
    feat["fragmented_systems"] = float(sum(v[2] for v in sysd.values()))
    # no systems.tsv (annotate --skip_systems, or proteins without ORF
    # coordinates): the sys_* features read as 0, which the salt-in detector
    # was not trained on -- flagged, not a model feature
    feat["systems_available"] = 1.0 if sysd else 0.0
    return feat


def check_required(feat: dict[str, float | None], gene_counts: str | Path) -> None:
    missing = [k for k in REQUIRED_HEADER if feat.get(k) is None]
    if missing:
        raise ValueError(
            f"{gene_counts} has no {', '.join('# ' + m for m in missing)} header "
            "line(s); NaCl prediction needs them. Re-run `osmotool annotate` with "
            "osmotool >= 0.8.0."
        )


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def salt_in_probability(detector: dict[str, Any], feat: dict[str, float | None]) -> float:
    z = detector["intercept"]
    for n, m, s, c in zip(detector["features"], detector["mean"], detector["scale"],
                          detector["coef"]):
        v = feat.get(n)
        x = math.log1p(max(v, 0.0)) if v is not None else m
        z += c * (x - m) / s
    return 1 / (1 + math.exp(-z))


def prediction_interval(schema: dict[str, Any], pred: float) -> tuple[float, float]:
    bins = schema["intervals"]["bins"]
    for b in bins:
        lo = -math.inf if b["pred_gt"] is None else b["pred_gt"]
        hi = math.inf if b["pred_le"] is None else b["pred_le"]
        if lo < pred <= hi:
            return max(pred + b["resid_lo"], 0.0), pred + b["resid_hi"]
    b = bins[0]
    return max(pred + b["resid_lo"], 0.0), pred + b["resid_hi"]


def flags_for(schema: dict[str, Any], feat: dict[str, float | None], p_salt_in: float) -> list[str]:
    app = schema["applicability"]
    names = schema["features"]
    flags = []
    if p_salt_in >= schema["salt_in_detector"]["threshold"]:
        flags.append("salt_in_like")
    n_out = sum(1 for n in names if feat.get(n) is not None
                and not (app["feature_q01"][n] <= feat[n] <= app["feature_q99"][n]))
    if n_out > app["max_features_outside"]:
        flags.append(f"out_of_range({n_out})")
    if (feat.get("n_contigs") or 0) > app["n_contigs_warn"]:
        flags.append("many_contigs")
    if (feat.get("total_proteins") or 0) < app["min_total_proteins"]:
        flags.append("few_proteins")
    if feat.get("systems_available") == 0.0:
        flags.append("no_systems")
    flags += [f"missing:{n}" for n in names if feat.get(n) is None]
    return flags


@dataclass
class NaclPrediction:
    genome: str
    pred: float
    lo: float
    hi: float
    p_salt_in: float
    flags: list[str]
    base_value: float
    contributions: list[tuple[str, float, float]]   # (feature, value, shap), by |shap| desc

    def top(self, n: int = N_TOP) -> str:
        return ";".join(f"{f}={v:g}:{s:+.2f}" for f, v, s in self.contributions[:n])


def predict(model: NaclModel, feats: list[dict[str, float | None]],
            genomes: list[str]) -> list[NaclPrediction]:
    """Predictions with exact TreeSHAP contributions (XGBoost pred_contribs)."""
    import numpy as np
    xgb = require_xgboost()
    names = model.features
    X = np.array([[np.nan if f.get(n) is None else f[n] for n in names] for f in feats],
                 dtype=float)
    contrib = model.booster.predict(xgb.DMatrix(X, feature_names=names), pred_contribs=True)
    target = model.schema["target"]
    out = []
    for i, (f, g) in enumerate(zip(feats, genomes)):
        raw = float(contrib[i].sum())
        if target["transform"] == "sqrt":
            pred = max(raw, 0.0) ** 2
        else:
            pred = max(raw, target.get("clip_min", 0.0))
        lo, hi = prediction_interval(model.schema, pred)
        p_si = salt_in_probability(model.schema["salt_in_detector"], f)
        sv = contrib[i, :-1]
        order = sorted(range(len(names)), key=lambda j: -abs(sv[j]))
        out.append(NaclPrediction(
            genome=g, pred=pred, lo=lo, hi=hi, p_salt_in=p_si,
            flags=flags_for(model.schema, f, p_si),
            base_value=float(contrib[i, -1]),
            contributions=[(names[j], float(X[i, j]), float(sv[j])) for j in order],
        ))
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

PREDICTION_COLUMNS = ["genome", "nacl_max_pred", "nacl_max_lo", "nacl_max_hi", "p_salt_in",
                      "flags", "base_value", "top_contributions"]


def _header_lines(model: NaclModel) -> list[str]:
    s = model.schema
    cv = s["cv"]["genus_blocked"]
    iv = s["intervals"]
    return [
        f"# model\tnacl_max (osmo_refdb {s['refdb_release']}, annotate --method "
        f"{s['annotate']['method']})",
        "# predicts\thighest NaCl (% w/v) a salt-out prokaryote grows at, as BacDive records it",
        f"# interval\t{int(round((iv['alpha_hi'] - iv['alpha_lo']) * 100))}% split-conformal "
        "(nacl_max_lo, nacl_max_hi)",
        f"# cv_genus_blocked\tSpearman {cv['spearman']:.3f}, MAE {cv['mae']:.2f} % NaCl "
        "(genomes from genera the model never saw)",
        "# note\tpredictions are compressed toward the middle (low-tolerance strains "
        "over-, extreme strains under-predicted); salt-in organisms are out of domain "
        "(flag salt_in_like). See the model's README.md.",
    ]


def _row(p: NaclPrediction) -> list[str]:
    return [p.genome, f"{p.pred:.2f}", f"{p.lo:.2f}", f"{p.hi:.2f}", f"{p.p_salt_in:.3f}",
            ",".join(p.flags) or "ok", f"{p.base_value:.3f}", p.top()]


def write_prediction_table(path: str | Path, model: NaclModel,
                           preds: list[NaclPrediction]) -> Path:
    path = Path(path)
    with path.open("w", newline="") as fh:
        for line in _header_lines(model):
            fh.write(line + "\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(PREDICTION_COLUMNS)
        for p in preds:
            w.writerow(_row(p))
    return path


def write_nacl_outputs(out_prefix: str | Path, model: NaclModel,
                       pred: NaclPrediction) -> tuple[Path, Path]:
    """<prefix>.nacl_prediction.tsv (one row) and <prefix>.nacl_contributions.tsv
    (every feature's value and SHAP contribution, largest first; base_value plus
    the shap column sums to the unclipped prediction)."""
    pred_path = write_prediction_table(f"{out_prefix}.nacl_prediction.tsv", model, [pred])
    contrib_path = Path(f"{out_prefix}.nacl_contributions.tsv")
    with contrib_path.open("w", newline="") as fh:
        fh.write(f"# base_value\t{pred.base_value:.4f}\n")
        fh.write(f"# nacl_max_pred\t{pred.pred:.4f}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["feature", "value", "shap"])
        for f, v, s in pred.contributions:
            w.writerow([f, f"{v:g}", f"{s:.4f}"])
    return pred_path, contrib_path


def log_prediction(pred: NaclPrediction) -> None:
    log.info("NaCl max prediction (%s): %.1f%% [%.1f-%.1f], flags: %s; largest "
             "contributions: %s", pred.genome, pred.pred, pred.lo, pred.hi,
             ",".join(pred.flags) or "ok", pred.top(3))
    if "salt_in_like" in pred.flags:
        log.warning("%s looks like a salt-in organism (p=%.2f); the NaCl prediction is "
                    "not reliable for salt-in lineages.", pred.genome, pred.p_salt_in)


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def predict_from_prefix(model: NaclModel, prefix: str | Path,
                        gene_counts: str | Path | None = None,
                        systems: str | Path | None = None) -> NaclPrediction:
    """Predict for one annotate output prefix and write its two output files."""
    prefix = Path(prefix)
    gc_path = Path(gene_counts) if gene_counts else Path(f"{prefix}.gene_counts.tsv")
    sys_path = Path(systems) if systems else Path(f"{prefix}.systems.tsv")
    feat = features_from_annotate(gc_path, sys_path, model.schema)
    check_required(feat, gc_path)
    pred = predict(model, [feat], [prefix.name])[0]
    write_nacl_outputs(prefix, model, pred)
    log_prediction(pred)
    return pred


def annotate_method_and_release(prefix: str | Path) -> tuple[str | None, str | None]:
    """(method, release dir name) from <prefix>.aln_stats.tsv, if present.
    The release name is the directory above hmms/ or holding osmo_refdb.dmnd."""
    path = Path(f"{prefix}.aln_stats.tsv")
    if not path.exists():
        return None, None
    stats = {}
    with path.open() as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if len(r) >= 2:
                stats[r[0]] = r[1]
    mode = stats.get("mode", "")
    method = mode.split(":", 1)[1] if mode.startswith("annotate:") else None
    db = stats.get("db")
    release = None
    if db:
        p = Path(db)
        release = p.parent.parent.name if p.parent.name == "hmms" else p.parent.name
    return method, release
