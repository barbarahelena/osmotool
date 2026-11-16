"""
test_nacl.py — NaCl-max prediction (nacl.py, annotate --predict_nacl, predict-nacl)

Uses a tiny synthetic model and schema built here, so the tests do not need
a real osmo_refdb release. Parity with osmo_refdb's reference implementation
(pipeline/nacl_model/predict.py) on real genomes is checked outside CI.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from osmotool import nacl
from osmotool.__main__ import main
from osmotool.utils import resolve_refdb

FAMILIES = ["betL", "ectA", "teaA"]
SYSTEMS = ["ect"]
FEATURES = ["total_proteins", "gc", "n_contigs", "betL", "ectA", "teaA"]
DETECTOR_FEATURES = ["total_proteins", "gc", "betL", "ectA", "teaA",
                     "sys_ect_present", "sys_ect_subunits"]


# ---------------------------------------------------------------------------
# fixtures: annotate output, schema, model, release directory
# ---------------------------------------------------------------------------

def write_annotate_output(prefix: Path, *, betl=2.0, ecta=1.0, header=True,
                          systems=True, method="hmm", db="/refdb/vTest/hmms/osmo_refdb.hmm"):
    lines = ["# total_proteins\t4000"]
    if header:
        lines += ["# genome_bp\t4200000", "# gc\t0.52000", "# n_contigs\t12",
                  "# contig_n50\t900000"]
    lines += ["# filtered_reads\t40", "gene\traw_count\tCOPIES_PER_KB",
              f"betL\t{betl}\t1.0", f"ectA\t{ecta}\t1.0", "betL_decoy\t5.0\t1.0"]
    # teaA deliberately absent (osmotool <= 0.8 omits it when zero)
    Path(f"{prefix}.gene_counts.tsv").write_text("\n".join(lines) + "\n")
    if systems:
        Path(f"{prefix}.systems.tsv").write_text(
            "system\tpresent\tsubunit_count\tcolocalized\tfragmented_assembly_flag\t"
            "contig\tmandatory_hits\taccessory_hits\trescued_hits\n"
            "ect\tTrue\t3\tTrue\tFalse\tc1\tectA,ectB,ectC\t-\t-\n")
    Path(f"{prefix}.aln_stats.tsv").write_text(
        f"metric\tvalue\nmode\tannotate:{method}\ndb\t{db}\n")


def make_schema() -> dict:
    return {
        "model": "nacl_max", "schema_version": 1, "model_file": "nacl_max.xgb.json",
        "refdb_release": "vTest",
        "annotate": {"method": "hmm"},
        "features": FEATURES, "families": FAMILIES, "systems": SYSTEMS,
        "target": {"name": "nacl_max_pct", "transform": "raw", "clip_min": 0.0},
        "cv": {"genus_blocked": {"spearman": 0.5, "mae": 2.0}},
        "intervals": {"alpha_lo": 0.1, "alpha_hi": 0.9, "bins": [
            {"pred_gt": None, "pred_le": 3.0, "resid_lo": -2.0, "resid_hi": 2.0},
            {"pred_gt": 3.0, "pred_le": None, "resid_lo": -4.0, "resid_hi": 5.0},
        ]},
        "applicability": {
            "feature_q01": {f: 0.0 for f in FEATURES},
            "feature_q99": {**{f: 10.0 for f in FEATURES}, "total_proteins": 9000.0,
                            "n_contigs": 500.0},
            "max_features_outside": 0, "min_total_proteins": 500.0, "n_contigs_warn": 400.0,
        },
        "salt_in_detector": {
            "features": DETECTOR_FEATURES,
            "mean": [8.0, 0.4, 0.5, 0.3, 0.1, 0.1, 0.5],
            "scale": [0.5, 0.1, 0.5, 0.5, 0.3, 0.3, 0.8],
            "coef": [0.2, 1.0, -1.0, -0.5, 0.0, -0.3, 0.1],
            "intercept": -2.0, "threshold": 0.5,
        },
    }


@pytest.fixture
def model_dir(tmp_path: Path) -> Path:
    xgb = pytest.importorskip("xgboost")
    np = pytest.importorskip("numpy")
    rng = np.random.default_rng(0)
    n = 300
    X = np.column_stack([rng.integers(1000, 8000, n), rng.uniform(0.3, 0.7, n),
                         rng.integers(1, 300, n), rng.integers(0, 4, n),
                         rng.integers(0, 2, n), rng.integers(0, 2, n)]).astype(float)
    y = 2 + 2.5 * X[:, 3] + 3 * X[:, 4] + rng.normal(0, 0.5, n)
    # xgboost's native API: XGBRegressor would need scikit-learn, which
    # neither osmotool nor its [dev] extra installs (CI has no sklearn)
    b = xgb.train({"max_depth": 3, "seed": 0, "objective": "reg:squarederror"},
                  xgb.DMatrix(X, label=y, feature_names=FEATURES), num_boost_round=20)
    d = tmp_path / "vTest" / "nacl_model"
    d.mkdir(parents=True)
    b.save_model(str(d / "nacl_max.xgb.json"))
    (d / nacl.SCHEMA_FILE).write_text(json.dumps(make_schema()))
    rows = X[:3].tolist()
    margin = b.predict(xgb.DMatrix(X[:3], feature_names=FEATURES), output_margin=True)
    (d / nacl.SELF_CHECK_FILE).write_text(json.dumps({
        "xgboost_version": xgb.__version__, "tolerance": 1e-4, "features": FEATURES,
        "rows": rows, "output_margin": [float(v) for v in margin]}))
    (tmp_path / "vTest" / "osmo_refdb.dmnd").write_text("")
    return d


# ---------------------------------------------------------------------------
# features
# ---------------------------------------------------------------------------

def test_read_gene_counts_header_and_decoys(tmp_path):
    write_annotate_output(tmp_path / "g")
    header, counts = nacl.read_gene_counts(tmp_path / "g.gene_counts.tsv")
    assert header["total_proteins"] == "4000" and header["gc"] == "0.52000"
    assert counts == {"betL": 2.0, "ectA": 1.0}          # decoy skipped


def test_features_absent_family_is_zero_and_systems(tmp_path):
    write_annotate_output(tmp_path / "g")
    f = nacl.features_from_annotate(tmp_path / "g.gene_counts.tsv",
                                    tmp_path / "g.systems.tsv", make_schema())
    assert f["teaA"] == 0.0
    assert f["total_proteins"] == 4000 and f["gc"] == 0.52 and f["n_contigs"] == 12
    assert f["genome_mb"] == pytest.approx(4.2)
    assert f["n_families_present"] == 2 and f["total_osmo_copies"] == 3.0
    assert f["sys_ect_present"] == 1.0 and f["sys_ect_subunits"] == 3.0
    assert f["systems_available"] == 1.0


def test_features_without_systems_file(tmp_path):
    write_annotate_output(tmp_path / "g", systems=False)
    f = nacl.features_from_annotate(tmp_path / "g.gene_counts.tsv",
                                    tmp_path / "g.systems.tsv", make_schema())
    assert f["sys_ect_present"] == 0.0 and f["systems_available"] == 0.0


def test_missing_header_is_an_error(tmp_path):
    write_annotate_output(tmp_path / "g", header=False)
    path = tmp_path / "g.gene_counts.tsv"
    f = nacl.features_from_annotate(path, None, make_schema())
    with pytest.raises(ValueError, match="osmotool >= 0.8.0"):
        nacl.check_required(f, path)


# ---------------------------------------------------------------------------
# interval, salt-in detector, flags, schema checks
# ---------------------------------------------------------------------------

def test_prediction_interval_bins_and_clip():
    s = make_schema()
    assert nacl.prediction_interval(s, 1.0) == (0.0, 3.0)      # lo clipped at 0
    assert nacl.prediction_interval(s, 3.0) == (1.0, 5.0)      # bin edge: (.., 3]
    assert nacl.prediction_interval(s, 8.0) == (4.0, 13.0)


def test_salt_in_probability_matches_formula():
    s = make_schema()
    det = s["salt_in_detector"]
    feat = {"total_proteins": 2000.0, "gc": 0.65, "betL": 0.0, "ectA": 0.0, "teaA": 0.0,
            "sys_ect_present": 0.0, "sys_ect_subunits": 0.0}
    z = det["intercept"] + sum(
        c * (math.log1p(feat[n]) - m) / sd
        for n, m, sd, c in zip(det["features"], det["mean"], det["scale"], det["coef"]))
    assert nacl.salt_in_probability(det, feat) == pytest.approx(1 / (1 + math.exp(-z)))


def test_flags():
    s = make_schema()
    feat = {"total_proteins": 400.0, "gc": 0.5, "n_contigs": 450.0, "betL": 12.0,
            "ectA": 0.0, "teaA": 0.0, "systems_available": 0.0}
    flags = nacl.flags_for(s, feat, p_salt_in=0.9)
    assert flags == ["salt_in_like", "out_of_range(1)", "many_contigs", "few_proteins",
                     "no_systems"]


def test_schema_version_and_method_checked(tmp_path):
    s = make_schema()
    s["schema_version"] = 99
    (tmp_path / nacl.SCHEMA_FILE).write_text(json.dumps(s))
    with pytest.raises(ValueError, match="schema_version"):
        nacl.load_schema(tmp_path)
    s["schema_version"], s["annotate"]["method"] = 1, "diamond"
    (tmp_path / nacl.SCHEMA_FILE).write_text(json.dumps(s))
    with pytest.raises(ValueError, match="annotate method"):
        nacl.load_schema(tmp_path)


# ---------------------------------------------------------------------------
# prediction and output files
# ---------------------------------------------------------------------------

def test_predict_shap_sums_to_prediction_and_outputs(tmp_path, model_dir):
    model = nacl.load_nacl_model(model_dir)
    write_annotate_output(tmp_path / "g")
    pred = nacl.predict_from_prefix(model, tmp_path / "g")

    raw = pred.base_value + sum(s for _, _, s in pred.contributions)
    assert pred.pred == pytest.approx(max(raw, 0.0), abs=1e-4)
    # matches xgboost's own prediction
    import numpy as np
    import xgboost as xgb
    x = np.array([[4000, 0.52, 12, 2, 1, 0]], dtype=float)
    direct = model.booster.predict(xgb.DMatrix(x, feature_names=FEATURES))[0]
    assert pred.pred == pytest.approx(max(float(direct), 0.0), abs=1e-4)
    # contributions sorted by |shap|, values carried through
    shaps = [abs(s) for _, _, s in pred.contributions]
    assert shaps == sorted(shaps, reverse=True)
    assert dict((f, v) for f, v, _ in pred.contributions)["teaA"] == 0.0
    assert pred.lo <= pred.pred <= pred.hi

    table = (tmp_path / "g.nacl_prediction.tsv").read_text().splitlines()
    assert table[0].startswith("# model\tnacl_max (osmo_refdb vTest")
    rows = [line for line in table if not line.startswith("#")]
    assert rows[0].split("\t") == nacl.PREDICTION_COLUMNS
    assert rows[1].split("\t")[0] == "g"
    contrib = (tmp_path / "g.nacl_contributions.tsv").read_text().splitlines()
    assert contrib[2] == "feature\tvalue\tshap"
    assert len(contrib) == 3 + len(FEATURES)


def test_self_check_rejects_a_model_it_cannot_reproduce(model_dir):
    chk = json.loads((model_dir / nacl.SELF_CHECK_FILE).read_text())
    chk["output_margin"] = [v + 5.0 for v in chk["output_margin"]]   # a misread base score
    (model_dir / nacl.SELF_CHECK_FILE).write_text(json.dumps(chk))
    with pytest.raises(RuntimeError, match="does not reproduce"):
        nacl.load_nacl_model(model_dir)


def test_self_check_optional(model_dir):
    (model_dir / nacl.SELF_CHECK_FILE).unlink()
    assert nacl.load_nacl_model(model_dir).release == "vTest"


def test_old_xgboost_refused(monkeypatch):
    xgb = pytest.importorskip("xgboost")
    monkeypatch.setattr(xgb, "__version__", "3.0.5")
    with pytest.raises(RuntimeError, match="xgboost >= 3.1"):
        nacl.require_xgboost()


def test_version_tuple():
    assert nacl._version_tuple("3.1.0") == (3, 1)
    assert nacl._version_tuple("3.2.0rc1") == (3, 2)
    assert nacl._version_tuple("2.0.3") < nacl.MIN_XGBOOST


# ---------------------------------------------------------------------------
# release resolution and CLI
# ---------------------------------------------------------------------------

def test_resolve_refdb_finds_nacl_model(model_dir):
    root = model_dir.parent
    assert resolve_refdb(str(root), mode="annotate").nacl_model == model_dir
    assert resolve_refdb(str(root), mode="profile").nacl_model is None


def test_resolve_refdb_without_nacl_model(tmp_path):
    (tmp_path / "osmo_refdb.dmnd").write_text("")
    assert resolve_refdb(str(tmp_path), mode="annotate").nacl_model is None


def test_resolve_refdb_schema_without_model_file(tmp_path):
    (tmp_path / "osmo_refdb.dmnd").write_text("")
    (tmp_path / "nacl_model").mkdir()
    (tmp_path / "nacl_model" / nacl.SCHEMA_FILE).write_text(json.dumps(make_schema()))
    with pytest.raises(FileNotFoundError, match="nacl_max.xgb.json"):
        resolve_refdb(str(tmp_path), mode="annotate")


def test_cli_predict_nacl(tmp_path, model_dir, capsys):
    root = model_dir.parent
    write_annotate_output(tmp_path / "a", betl=0.0, ecta=0.0)
    write_annotate_output(tmp_path / "b", betl=3.0, ecta=1.0)
    out = tmp_path / "summary.tsv"
    assert main(["predict-nacl", str(root), str(tmp_path / "a"), str(tmp_path / "b"),
                 "-o", str(out)]) == 0
    rows = [r.split("\t") for r in out.read_text().splitlines() if not r.startswith("#")]
    assert [r[0] for r in rows[1:]] == ["a", "b"]
    pred = {r[0]: float(r[1]) for r in rows[1:]}
    assert pred["b"] > pred["a"]           # the synthetic model rewards betL/ectA
    assert (tmp_path / "a.nacl_contributions.tsv").exists()


def test_cli_predict_nacl_rejects_diamond_output(tmp_path, model_dir, capsys):
    write_annotate_output(tmp_path / "a", method="diamond")
    assert main(["predict-nacl", str(model_dir.parent), str(tmp_path / "a")]) == 1
    assert "needs HMM counts" in capsys.readouterr().err


def test_cli_predict_nacl_warns_on_other_release(tmp_path, model_dir, caplog):
    write_annotate_output(tmp_path / "a", db="/refdb/v9/hmms/osmo_refdb.hmm")
    assert main(["predict-nacl", str(model_dir.parent), str(tmp_path / "a"),
                 "-o", str(tmp_path / "s.tsv")]) == 0
    assert "annotated with a database in 'v9/'" in caplog.text


def test_cli_annotate_predict_nacl_needs_hmm(tmp_path, model_dir):
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c1\nACGT\n")
    with pytest.raises(SystemExit):
        main(["annotate", str(model_dir.parent), str(fasta), "--method", "diamond",
              "--predict_nacl", "--out_prefix", str(tmp_path / "g")])


def test_cli_annotate_predict_nacl_needs_model(tmp_path):
    (tmp_path / "osmo_refdb.dmnd").write_text("")
    (tmp_path / "hmms").mkdir()
    (tmp_path / "hmms" / "osmo_refdb.hmm").write_text("")
    fasta = tmp_path / "g.fna"
    fasta.write_text(">c1\nACGT\n")
    with pytest.raises(SystemExit):
        main(["annotate", str(tmp_path), str(fasta), "--predict_nacl",
              "--out_prefix", str(tmp_path / "g")])
