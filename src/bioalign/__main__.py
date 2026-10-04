"""Command line entry point.

    python3 demo.py
    PYTHONPATH=src python3 -m bioalign
    PYTHONPATH=src python3 -m bioalign align GATTACA GCATGCT --matrix
    PYTHONPATH=src python3 -m bioalign patterns ATGCGTACGT
    PYTHONPATH=src python3 -m bioalign fasta data/sample_sequences.fasta
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from bioalign.alignment import format_alignment, format_matrix
from bioalign.demo import run_demo
from bioalign.needleman_wunsch import global_score_matrix, needleman_wunsch
from bioalign.scoring import Scoring
from bioalign.sequence import clean_sequence, read_fasta
from bioalign.smith_waterman import local_score_matrix, smith_waterman
from bioalign.structure import format_report, predict_patterns


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "align":
            print(_align(args))
        elif args.command == "patterns":
            print(format_report(predict_patterns(args.sequence, args.window, args.step)))
        elif args.command == "fasta":
            print(_fasta(args.path))
        else:
            print(run_demo())
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioalign",
        description="Align two DNA sequences, or scan one sequence for structural patterns.",
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("demo", help="run the four short lessons (this is the default)")

    align = sub.add_parser("align", help="globally and locally align two DNA strings")
    align.add_argument("sequence_a")
    align.add_argument("sequence_b")
    align.add_argument("--matrix", action="store_true", help="print both score grids")
    align.add_argument("--match", type=int, default=2)
    align.add_argument("--mismatch", type=int, default=-1)
    align.add_argument("--gap", type=int, default=-2)

    patterns = sub.add_parser("patterns", help="scan one DNA string for structural patterns")
    patterns.add_argument("sequence")
    patterns.add_argument("--window", type=int, default=40)
    patterns.add_argument("--step", type=int, default=10)

    fasta = sub.add_parser("fasta", help="align the first two FASTA records and scan every record")
    fasta.add_argument("path")
    return parser


def _align(args: argparse.Namespace) -> str:
    scoring = Scoring(match=args.match, mismatch=args.mismatch, gap=args.gap)
    sequence_a = clean_sequence(args.sequence_a, "sequence A")
    sequence_b = clean_sequence(args.sequence_b, "sequence B")
    blocks = [
        format_alignment(needleman_wunsch(sequence_a, sequence_b, scoring)),
        "",
        format_alignment(smith_waterman(sequence_a, sequence_b, scoring)),
    ]
    if args.matrix:
        blocks.extend(
            [
                "",
                "Global score grid",
                format_matrix(
                    sequence_a,
                    sequence_b,
                    global_score_matrix(sequence_a, sequence_b, scoring),
                ),
                "",
                "Local score grid",
                format_matrix(
                    sequence_a,
                    sequence_b,
                    local_score_matrix(sequence_a, sequence_b, scoring),
                ),
            ]
        )
    return "\n".join(blocks)


def _fasta(path: str) -> str:
    text = Path(path).read_text()
    records = read_fasta(text)
    blocks: list[str] = []
    if len(records) >= 2:
        name_a, sequence_a = records[0]
        name_b, sequence_b = records[1]
        blocks.append(f"Aligning {name_a} vs {name_b}")
        blocks.append(format_alignment(needleman_wunsch(sequence_a, sequence_b)))
        blocks.append("")
        blocks.append(format_alignment(smith_waterman(sequence_a, sequence_b)))
        blocks.append("")
    for name, sequence in records:
        blocks.append(f"Patterns in {name}")
        blocks.append(format_report(predict_patterns(sequence)))
        blocks.append("")
    return "\n".join(blocks).rstrip()


if __name__ == "__main__":
    raise SystemExit(main())
