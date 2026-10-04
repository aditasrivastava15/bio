"""Smith-Waterman: find the best matching stretch, and ignore the rest.

The grid is the same idea as Needleman-Wunsch, with one extra choice in every
cell: score 0, which means "throw away the alignment so far and start over."

Because a fresh start is free, the first row and first column stay 0. There is
no punishment for skipping letters before or after the matching stretch.

The best local alignment is not stored in the corner. It is the highest number
anywhere in the grid. The walk starts at that cell and stops when it reaches 0.
"""

from __future__ import annotations

from bioalign.alignment import Alignment
from bioalign.scoring import Scoring
from bioalign.sequence import clean_sequence


def smith_waterman(
    sequence_a: str,
    sequence_b: str,
    scoring: Scoring | None = None,
) -> Alignment:
    """Locally align two DNA sequences. Only the best shared stretch is returned."""
    prepared_a, prepared_b, scoring = _prepare(sequence_a, sequence_b, scoring)
    matrix = local_score_matrix(prepared_a, prepared_b, scoring)
    end_i, end_j, best = _best_cell(matrix)
    if best == 0:
        aligned_a, aligned_b = "", ""
    else:
        aligned_a, aligned_b = _traceback(
            prepared_a, prepared_b, matrix, scoring, end_i, end_j
        )
    return Alignment(
        method="local",
        sequence_a=prepared_a,
        sequence_b=prepared_b,
        aligned_a=aligned_a,
        aligned_b=aligned_b,
        score=best,
    )


def local_score_matrix(
    sequence_a: str,
    sequence_b: str,
    scoring: Scoring | None = None,
) -> list[list[int]]:
    """Fill the Smith-Waterman score grid. Borders stay 0 so a start is free."""
    sequence_a, sequence_b, scoring = _prepare(sequence_a, sequence_b, scoring)
    rows = len(sequence_a) + 1
    cols = len(sequence_b) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    for i in range(1, rows):
        for j in range(1, cols):
            diagonal = matrix[i - 1][j - 1] + scoring.pair(
                sequence_a[i - 1], sequence_b[j - 1]
            )
            up = matrix[i - 1][j] + scoring.gap
            left = matrix[i][j - 1] + scoring.gap
            restart = 0
            matrix[i][j] = max(restart, diagonal, up, left)

    return matrix


def _best_cell(matrix: list[list[int]]) -> tuple[int, int, int]:
    """Return the row, column, and score of the best cell.

    If several cells tie, the earliest row and then the earliest column wins.
    That keeps the result the same every time the program runs.
    """
    best_score = 0
    best_i = 0
    best_j = 0
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if value > best_score:
                best_score = value
                best_i = i
                best_j = j
    return best_i, best_j, best_score


def _traceback(
    sequence_a: str,
    sequence_b: str,
    matrix: list[list[int]],
    scoring: Scoring,
    i: int,
    j: int,
) -> tuple[str, str]:
    """Walk backward from the best cell until the score falls to 0."""
    aligned_a: list[str] = []
    aligned_b: list[str] = []

    while i > 0 and j > 0 and matrix[i][j] > 0:
        diagonal = matrix[i - 1][j - 1] + scoring.pair(
            sequence_a[i - 1], sequence_b[j - 1]
        )
        if matrix[i][j] == diagonal:
            aligned_a.append(sequence_a[i - 1])
            aligned_b.append(sequence_b[j - 1])
            i -= 1
            j -= 1
            continue
        if matrix[i][j] == matrix[i - 1][j] + scoring.gap:
            aligned_a.append(sequence_a[i - 1])
            aligned_b.append("-")
            i -= 1
            continue
        if matrix[i][j] == matrix[i][j - 1] + scoring.gap:
            aligned_a.append("-")
            aligned_b.append(sequence_b[j - 1])
            j -= 1
            continue
        raise RuntimeError("local traceback lost the path")

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
