"""The points the aligner adds or subtracts for each decision.

An alignment is a path through two sequences. At every step the path can:

- pair two letters (reward a match, penalize a mismatch)
- skip a letter in one sequence (a gap)

Those three numbers are the whole scoring system. Change them and the
"best" alignment can change, the same way a grading rubric changes a grade.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Scoring:
    """Linear scoring: one reward for a match, one cost for a mismatch, one cost per gap letter.

    ``gap`` is zero or negative. A gap of -2 means "each skipped letter costs 2 points."
    """

    match: int = 2
    mismatch: int = -1
    gap: int = -2

    def __post_init__(self) -> None:
        if self.gap > 0:
            raise ValueError("gap must be 0 or negative, for example -2")
        if self.match <= self.mismatch:
            raise ValueError("match score must be higher than mismatch score")

    def pair(self, left: str, right: str) -> int:
        """Points for lining one base up with another."""
        if left == right:
            return self.match
        return self.mismatch
