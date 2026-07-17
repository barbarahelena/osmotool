"""
__main__.py — osmotool CLI entry point
"""

from __future__ import annotations

import sys

from osmotool.handle_args import build_parser
from osmotool.runners import run_profile, run_annotate
from osmotool.utils import setup_logging, resolve_refdb


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    setup_logging(verbose=args.verbose)

    try:
        if args.subcommand == "profile":
            if args.reads2 is not None and args.reads1 is None:
                parser.error("-2/--reads2 requires -1/--reads1")
            refdb = resolve_refdb(args.database, mode="profile")
            if args.cascade and refdb.hmm is None:
                parser.error(
                    f"--cascade requires hmms/osmo_refdb.hmm in {args.database}, "
                    "but it wasn't found"
                )
            if args.cascade and refdb.cascade_config is None:
                parser.error(
                    f"--cascade requires osmo_refdb.profile_cascade.tsv in "
                    f"{args.database}, but it wasn't found"
                )
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
                total_reads=args.total_reads,
                hmm_db=str(refdb.hmm) if args.cascade else None,
                cascade_config=str(refdb.cascade_config) if args.cascade else None,
                exclude_families=str(refdb.exclude_families) if refdb.exclude_families else None,
            )

        elif args.subcommand == "annotate":
            refdb = resolve_refdb(args.database, mode="annotate")
            if args.method in ("hmm", "both") and refdb.hmm is None:
                parser.error(
                    f"--method {args.method} requires hmms/osmo_refdb.hmm in "
                    f"{args.database}, but it wasn't found (pass --method diamond "
                    "instead if you don't have it)"
                )
            run_annotate(
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
            )

    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
