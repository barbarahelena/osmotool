"""
__main__.py — osmotool CLI entry point
"""

from __future__ import annotations

import sys

from osmotool.handle_args import build_parser
from osmotool.runners import run_profile, run_annotate
from osmotool.utils import setup_logging


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    setup_logging(verbose=args.verbose)

    try:
        if args.subcommand == "profile":
            if args.reads2 is not None and args.reads1 is None:
                parser.error("-2/--reads2 requires -1/--reads1")
            run_profile(
                args.db,
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
                hmm_db=args.hmm_db,
                cascade_config=args.cascade_config,
                exclude_families=args.exclude_families,
            )

        elif args.subcommand == "annotate":
            run_annotate(
                args.db,
                args.assembly,
                args.out_prefix,
                proteins=args.proteins,
                method=args.method,
                hmm_db=args.hmm_db,
                diamond_cutoffs=args.diamond_cutoffs,
                min_identity=args.min_identity,
                min_query_cover=args.min_query_cover,
                min_subject_cover=args.min_subject_cover,
                evalue=args.evalue,
                threads=args.threads,
                keep_aln=args.keep_aln,
                keep_proteins=args.keep_proteins,
                tmpdir=args.tmpdir,
            )

    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
