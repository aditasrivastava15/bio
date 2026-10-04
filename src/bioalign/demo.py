"""Printed lessons. Each one runs the real functions on a tiny example."""

from __future__ import annotations

from bioalign.alignment import format_alignment, format_matrix
from bioalign.needleman_wunsch import global_score_matrix, needleman_wunsch
from bioalign.samples import (
    FLANKED_A,
    FLANKED_B,
    INDEL_A,
    INDEL_B,
    TINY_A,
    TINY_B,
    pattern_sequence,
)
from bioalign.scoring import Scoring
from bioalign.smith_waterman import local_score_matrix, smith_waterman
from bioalign.structure import format_report, predict_patterns


def run_demo() -> str:
    scoring = Scoring()
    parts = [
        _lesson_grid(TINY_A, TINY_B, scoring),
        _lesson_local_vs_global(FLANKED_A, FLANKED_B, scoring),
        _lesson_indel(INDEL_A, INDEL_B, scoring),
        _lesson_patterns(pattern_sequence()),
    ]
    return "\n\n".join(parts)


def _lesson_grid(sequence_a: str, sequence_b: str, scoring: Scoring) -> str:
    global_matrix = global_score_matrix(sequence_a, sequence_b, scoring)
    local_matrix = local_score_matrix(sequence_a, sequence_b, scoring)
    global_alignment = needleman_wunsch(sequence_a, sequence_b, scoring)
    local_alignment = smith_waterman(sequence_a, sequence_b, scoring)
    return "\n".join(
        [
            "1. Score grids for GATTACA vs GCATGCT",
            "   Match +2, mismatch -1, gap -2.",
            "   Rows are sequence A. Columns are sequence B.",
            "   The extra blank row and column are the empty prefixes.",
            "",
            "Global grid (Needleman-Wunsch). Bottom-right is the full-sequence score.",
            format_matrix(sequence_a, sequence_b, global_matrix),
            "",
            format_alignment(global_alignment),
            "",
            "Local grid (Smith-Waterman). The best cell can sit anywhere.",
            "Negative choices are replaced by 0, which means start over.",
            format_matrix(sequence_a, sequence_b, local_matrix),
            "",
            format_alignment(local_alignment),
            "",
            _interpret(global_alignment, local_alignment),
        ]
    )


def _lesson_local_vs_global(sequence_a: str, sequence_b: str, scoring: Scoring) -> str:
    global_alignment = needleman_wunsch(sequence_a, sequence_b, scoring)
    local_alignment = smith_waterman(sequence_a, sequence_b, scoring)
    return "\n".join(
        [
            "2. Same conserved motif, unrelated flanks",
            f"   A: {sequence_a}",
            f"   B: {sequence_b}",
            "   Both contain ATGCGTAAGCTT. The ends were written so they do not match.",
            "",
            format_alignment(global_alignment),
            "",
            format_alignment(local_alignment),
            "",
            _interpret(global_alignment, local_alignment),
        ]
    )


def _lesson_indel(sequence_a: str, sequence_b: str, scoring: Scoring) -> str:
    global_alignment = needleman_wunsch(sequence_a, sequence_b, scoring)
    local_alignment = smith_waterman(sequence_a, sequence_b, scoring)
    return "\n".join(
        [
            "3. One extra letter (an insertion)",
            f"   A: {sequence_a}",
            f"   B: {sequence_b}",
            "   B is A with a T inserted after the fourth letter.",
            "",
            format_alignment(global_alignment),
            "",
            format_alignment(local_alignment),
            "",
            _interpret(global_alignment, local_alignment),
        ]
    )


def _lesson_patterns(sequence: str) -> str:
    report = predict_patterns(sequence)
    return "\n".join(
        [
            "4. Statistical structural patterns on one synthetic sequence",
            "   Layout: 60 bp AT repeat, 60 bp CG repeat, an A/T run, then mixed AGCT.",
            "   This does not fold the DNA into 3D. It measures composition, a Markov chain,",
            "   stacking energy, and a few classic cutoffs, window by window.",
            "",
            format_report(report),
        ]
    )


def _interpret(global_alignment, local_alignment) -> str:
    if local_alignment.score == 0:
        return (
            "Read: there is no stretch with a positive score, so the local alignment is empty. "
            "The global alignment still lines the sequences up, because it is required to."
        )
    local_is_tighter = (
        local_alignment.identity > global_alignment.identity + 0.05
        and len(local_alignment.aligned_a) < len(global_alignment.aligned_a)
    )
    if local_is_tighter:
        return (
            "Read: the local alignment is the shared stretch only, so its identity is higher. "
            "The global alignment must also place the unmatched ends, so mismatches (dots) appear there."
        )
    if global_alignment.gaps and not local_alignment.gaps:
        return (
            "Read: the dash is the extra letter. Global alignment keeps every letter, so it "
            "pays for a gap. Local alignment can stop at the best matching piece and skip that gap."
        )
    if global_alignment.gaps:
        return (
            "Read: a dash is a gap, a bar is a match, and a dot is a mismatch. "
            "The gap is how the aligner accounts for an insertion or deletion."
        )
    return (
        "Read: a bar is a match and a dot is a mismatch. "
        "Global uses every letter. Local keeps the best-scoring stretch."
    )
