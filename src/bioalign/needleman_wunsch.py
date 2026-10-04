"""Needleman-Wunsch: line up two whole DNA sequences from end to end.

Picture two words written on the top and left of a spreadsheet. Every cell
asks the same question: what is the best score for the prefixes that end here?

The answer is the best of three choices:

- diagonal: pair the newest letter of A with the newest letter of B
- up: skip a letter of A (a gap in B)
- left: skip a letter of B (a gap in A)

The top row and left column are filled first, because lining a letter up with
nothing is a string of gaps. When the grid is full, the bottom-right cell is
the score for using every letter. Walking backward from that cell rebuilds
the alignment.
"""

from __future__ import annotations

from bioalign.alignment import Alignment
from bioalign.scoring import Scoring
from bioalign.sequence import clean_sequence


def needleman_wunsch(
    sequence_a: str,
    sequence_b: str,
    scoring: Scoring | None = None,
) -> Alignment:
    """Globally align two DNA sequences. Both sequences are used completely."""
    prepared_a, prepared_b, scoring = _prepare(sequence_a, sequence_b, scoring)
    matrix = global_score_matrix(prepared_a, prepared_b, scoring)
    aligned_a, aligned_b = _traceback(prepared_a, prepared_b, matrix, scoring)
    return Alignment(
        method="global",
        sequence_a=prepared_a,
        sequence_b=prepared_b,
        aligned_a=aligned_a,
        aligned_b=aligned_b,
        score=matrix[-1][-1],
    )


def global_score_matrix(
    sequence_a: str,
    sequence_b: str,
    scoring: Scoring | None = None,
) -> list[list[int]]:
    """Fill the Needleman-Wunsch score grid. Row 0 / column 0 are pure gaps."""
    sequence_a, sequence_b, scoring = _prepare(sequence_a, sequence_b, scoring)
    rows = len(sequence_a) + 1
    cols = len(sequence_b) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        matrix[i][0] = i * scoring.gap
    for j in range(1, cols):
        matrix[0][j] = j * scoring.gap

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = matrix[i - 1][j - 1] + scoring.pair(
                sequence_a[i - 1], sequence_b[j - 1]
            )
            up = matrix[i - 1][j] + scoring.gap
            left = matrix[i][j - 1] + scoring.gap
            matrix[i][j] = max(diagonal, up, left)

    return matrix


def _traceback(
    sequence_a: str,
    sequence_b: str,
    matrix: list[list[int]],
    scoring: Scoring,
) -> tuple[str, str]:
    """Walk from the bottom-right cell back to the top-left, recording the path.

    Ties are broken in a fixed order: prefer pairing letters, then a gap in B,
    then a gap in A. The same comparison the grid used, so the walk stays on
    a real best path.
    """
    i = len(sequence_a)
    j = len(sequence_b)
    aligned_a: list[str] = []
    aligned_b: list[str] = []

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diagonal = matrix[i - 1][j - 1] + scoring.pair(
                sequence_a[i - 1], sequence_b[j - 1]
            )
            if matrix[i][j] == diagonal:
                aligned_a.append(sequence_a[i - 1])
                aligned_b.append(sequence_b[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + scoring.gap:
            aligned_a.append(sequence_a[i - 1])
            aligned_b.append("-")
            i -= 1
            continue
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + scoring.gap:
            aligned_a.append("-")
            aligned_b.append(sequence_b[j - 1])
            j -= 1
            continue
        raise RuntimeError("global traceback lost the path")

    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b)


def _prepare(
    sequence_a: str,
    sequence_b: str,
    scoring: Scoring | None,
) -> tuple[str, str, Scoring]:
    return (
        clean_sequence(sequence_a, "sequence A"),
        clean_sequence(sequence_b, "sequence B"),
        scoring if scoring is not None else Scoring(),
    )
