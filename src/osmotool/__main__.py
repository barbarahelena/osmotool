"""
__main__.py — osmotool CLI entry point
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

from osmotool.handle_args import build_parser
from osmotool.runners import run_profile, run_annotate
from osmotool.download import download_refdb
from osmotool.utils import setup_logging, resolve_refdb, check_tool, choose_catalog
from osmotool.catalog import build_catalog, index_catalog, load_catalog
from osmotool.functions import fragments_from_reads, load_excluded_families
from osmotool.nacl import (annotate_method_and_release, load_nacl_model,
                           predict_from_prefix, write_prediction_table)

log = logging.getLogger("osmotool")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    setup_logging(verbose=args.verbose)

    try:
        if args.subcommand == "profile":
            if args.reads2 is not None and args.reads1 is None:
                parser.error("-2/--reads2 requires -1/--reads1")
            total_fragments = args.total_fragments
            if args.total_reads is not None:
                total_fragments = fragments_from_reads(
                    args.total_reads, paired=args.reads2 is not None)
            refdb = resolve_refdb(args.database, mode="profile")
            cascade_enabled = not args.no_cascade
            if cascade_enabled and refdb.hmm is None:
                parser.error(
                    f"the DIAMOND+HMM cascade (on by default) requires "
                    f"hmms/osmo_refdb.hmm in {args.database}, but it wasn't "
                    "found. Pass --no_cascade to run DIAMOND-only instead."
                )
            if cascade_enabled and refdb.cascade_config is None:
                parser.error(
                    f"the DIAMOND+HMM cascade (on by default) requires "
                    f"osmo_refdb.profile_cascade.tsv in {args.database}, but "
                    "it wasn't found. Pass --no_cascade to run DIAMOND-only "
                    "instead."
                )
            markers = None if args.no_markers else refdb.markers
            if markers is None and not args.no_markers:
                log.info("No markers/ in %s: copies_per_genome not reported "
                         "(release predates marker normalisation).", args.database)
            if markers is not None and refdb.train_refs is None:
                parser.error(
                    f"copies-per-genome normalisation needs "
                    f"osmo_refdb.train_refs.faa in {args.database} (for "
                    "reference lengths), but it wasn't found. Pass "
                    "--no_markers to skip normalisation."
                )
            catalog = None
            try:
                catalog_dir = choose_catalog(refdb, habitat=args.habitat,
                                             catalog_dir=args.catalog_dir,
                                             translated_only=args.translated_only)
            except ValueError as e:
                parser.error(str(e))
            if catalog_dir is not None:
                if refdb.hmm is None:
                    parser.error(f"the catalog needs hmms/osmo_refdb.hmm in {args.database} "
                                 "to check it was annotated with this release")
                catalog = load_catalog(catalog_dir, hmm_db=refdb.hmm,
                                       aligner_bin=args.aligner_bin,
                                       index_dir=args.catalog_index)
                if args.catalog_dir is None:
                    # shipped with the release, so versioned with it
                    catalog.release = refdb.release
                log.info("Assigning reads against the %s gene catalog (osmo_refdb %s) first; "
                         "the rest go to translated search.", catalog.name, catalog.release)
            elif args.catalog_index is not None:
                parser.error("--catalog_index needs the catalog, which --translated_only skips")
            elif args.catalog_only:
                parser.error("--catalog_only and --translated_only exclude each other")
            if refdb.report_levels is None:
                log.info("No osmo_refdb.profile_report_levels.tsv in %s: no unit_counts.tsv "
                         "(release predates report levels).", args.database)
            run_profile(
                str(refdb.dmnd),
                args.out_prefix,
                reads1=args.reads1,
                reads2=args.reads2,
                singles=args.singles,
                min_identity=args.min_identity,
                min_query_cover=args.min_query_cover,
                min_subject_cover=args.min_subject_cover,
                min_seqlen=args.min_seqlen,
                evalue=args.evalue,
                threads=args.threads,
                keep_aln=args.keep_aln,
                tmpdir=args.tmpdir,
                total_fragments=total_fragments,
                hmm_db=str(refdb.hmm) if cascade_enabled else None,
                cascade_config=str(refdb.cascade_config) if cascade_enabled else None,
                exclude_families=str(refdb.exclude_families) if refdb.exclude_families else None,
                markers=markers,
                train_refs=str(refdb.train_refs) if markers is not None else None,
                report_levels=str(refdb.report_levels) if refdb.report_levels else None,
                catalog=catalog,
                catalog_min_identity=args.catalog_min_identity,
                catalog_min_aligned=args.catalog_min_aligned,
                translated_fallback=not args.catalog_only,
            )

        elif args.subcommand == "annotate":
            refdb = resolve_refdb(args.database, mode="annotate")
            if args.method in ("hmm", "both") and refdb.hmm is None:
                parser.error(
                    f"--method {args.method} requires hmms/osmo_refdb.hmm in "
                    f"{args.database}, but it wasn't found (pass --method diamond "
                    "instead if you don't have it)"
                )
            nacl_model = None
            if args.predict_nacl:
                if args.method == "diamond":
                    parser.error("--predict_nacl uses the HMM counts the model was "
                                 "trained on: pass --method hmm (default) or both")
                if refdb.nacl_model is None:
                    parser.error(f"--predict_nacl needs nacl_model/ in {args.database}, "
                                 "which this release doesn't ship (osmo_refdb v10 and "
                                 "later do)")
                # load before the (long) annotation run, so a missing xgboost
                # or an unsupported model fails fast
                nacl_model = load_nacl_model(refdb.nacl_model)
            results = run_annotate(
                str(refdb.dmnd),
                args.assembly,
                args.out_prefix,
                proteins=args.proteins,
                method=args.method,
                hmm_db=str(refdb.hmm) if refdb.hmm else None,
                diamond_cutoffs=str(refdb.diamond_cutoffs) if refdb.diamond_cutoffs else None,
                min_identity=args.min_identity,
                min_query_cover=args.min_query_cover,
                min_subject_cover=args.min_subject_cover,
                evalue=args.evalue,
                threads=args.threads,
                keep_aln=args.keep_aln,
                keep_proteins=args.keep_proteins,
                tmpdir=args.tmpdir,
                exclude_families=str(refdb.exclude_families) if refdb.exclude_families else None,
                skip_systems=args.skip_systems,
                systems_config=args.systems_config,
            )
            if nacl_model is not None:
                counts_path, _, systems_path, _ = results["hmm"]
                prefix = str(counts_path)[: -len(".gene_counts.tsv")]
                predict_from_prefix(nacl_model, prefix, counts_path, systems_path)

        elif args.subcommand == "predict-nacl":
            refdb = resolve_refdb(args.database, mode="annotate")
            if refdb.nacl_model is None:
                parser.error(f"no nacl_model/ in {args.database}: this release doesn't "
                             "ship the NaCl model (osmo_refdb v10 and later do)")
            nacl_model = load_nacl_model(refdb.nacl_model)
            preds = []
            for prefix in args.prefixes:
                if not Path(f"{prefix}.gene_counts.tsv").exists():
                    raise FileNotFoundError(f"{prefix}.gene_counts.tsv not found -- PREFIX "
                                            "is an annotate --out_prefix")
                method, release = annotate_method_and_release(prefix)
                if method is not None and method != "hmm":
                    raise ValueError(f"{prefix} is annotate --method {method} output; the "
                                     "model needs HMM counts (for --method both output, "
                                     "pass PREFIX.hmm)")
                if release is not None and release != nacl_model.release:
                    log.warning("%s was annotated with a database in '%s/'; the model was "
                                "trained on osmo_refdb %s counts. Predictions are only valid "
                                "if that is the same release.", prefix, release,
                                nacl_model.release)
                preds.append(predict_from_prefix(nacl_model, prefix))
            if args.output:
                write_prediction_table(args.output, nacl_model, preds)
                log.info("NaCl predictions for %d genome(s): %s", len(preds), args.output)
            else:
                write_prediction_table("/dev/stdout", nacl_model, preds)

        elif args.subcommand == "build-catalog":
            refdb = resolve_refdb(args.database, mode="annotate")
            if refdb.hmm is None:
                parser.error(f"build-catalog needs hmms/osmo_refdb.hmm in {args.database}")
            check_tool("hmmsearch")
            check_tool(args.aligner_bin or args.aligner)
            decoys = (load_excluded_families(refdb.exclude_families)
                      if refdb.exclude_families else set())
            print(build_catalog(
                refdb.hmm, args.genes, args.out_dir, proteins=args.proteins,
                decoy_models=decoys, release=args.database, threads=args.threads,
                min_protein_length=args.min_protein_length, evalue=args.evalue,
                aligner=args.aligner, aligner_bin=args.aligner_bin, index=args.index,
            ))

        elif args.subcommand == "index-catalog":
            check_tool(args.aligner_bin or args.aligner)
            print(index_catalog(args.genes, args.out_dir, aligner=args.aligner,
                                aligner_bin=args.aligner_bin))

        elif args.subcommand == "download-db":
            out_dir = download_refdb(
                args.release,
                args.location,
                force=args.force,
                keep_archive=args.keep_archive,
            )
            print(out_dir)

    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
