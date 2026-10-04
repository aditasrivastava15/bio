"""Shared alignment result: the two lined-up strings, the score, and a printable view."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Alignment:
    """One finished alignment.

    ``aligned_a`` and ``aligned_b`` are the same length. A dash means "this
    sequence has no letter here" (a gap). ``score`` is the total points earned
    by the path. ``method`` is ``global`` or ``local``.
    """

    method: str
    sequence_a: str
    sequence_b: str
    aligned_a: str
    aligned_b: str
    score: int

    def __post_init__(self) -> None:
        if len(self.aligned_a) != len(self.aligned_b):
            raise ValueError("aligned sequences must be the same length")

    @property
    def matches(self) -> int:
        return sum(
            left == right and left != "-"
            for left, right in zip(self.aligned_a, self.aligned_b)
        )

    @property
    def mismatches(self) -> int:
        return sum(
            left != right and left != "-" and right != "-"
            for left, right in zip(self.aligned_a, self.aligned_b)
        )

    @property
    def gaps(self) -> int:
        return sum(
            left == "-" or right == "-"
            for left, right in zip(self.aligned_a, self.aligned_b)
        )

    @property
    def identity(self) -> float:
        """Share of columns where both letters exist and are the same.

        Empty alignments (no shared stretch worth keeping) report 0.
        """
        if not self.aligned_a:
            return 0.0
        return self.matches / len(self.aligned_a)

    def columns(self) -> str:
        """A middle row of symbols: ``|`` match, ``.`` mismatch, `` `` gap."""
        symbols: list[str] = []
        for left, right in zip(self.aligned_a, self.aligned_b):
            if left == "-" or right == "-":
                symbols.append(" ")
            elif left == right:
                symbols.append("|")
            else:
                symbols.append(".")
        return "".join(symbols)


def format_alignment(alignment: Alignment, width: int = 60) -> str:
    """Render an alignment in blocks so long sequences stay readable."""
    if width < 1:
        raise ValueError("width must be at least 1")

    title = "global" if alignment.method == "global" else "local"
    header = (
        f"{title} alignment   score {alignment.score}"
        f"   identity {alignment.identity:.1%}"
        f"   matches {alignment.matches}"
        f"   mismatches {alignment.mismatches}"
        f"   gaps {alignment.gaps}"
    )
    if not alignment.aligned_a:
        return header + "\n(no positive local alignment)"

    blocks: list[str] = [header]
    middle = alignment.columns()
    for start in range(0, len(alignment.aligned_a), width):
        stop = start + width
        blocks.append(f"A {alignment.aligned_a[start:stop]}")
        blocks.append(f"  {middle[start:stop]}")
        blocks.append(f"B {alignment.aligned_b[start:stop]}")
        blocks.append("")
    return "\n".join(blocks).rstrip()


def format_matrix(sequence_a: str, sequence_b: str, matrix: list[list[int]]) -> str:
    """Print the score grid the aligner filled in, with sequence letters as labels."""
    width = max(len(str(value)) for row in matrix for value in row)
    width = max(width, 2)

    def cell(value: str) -> str:
        return value.rjust(width)

    header = " ".join([cell("")] + [cell("")] + [cell(base) for base in sequence_b])
    rows = [header]
    for i, row in enumerate(matrix):
        label = " " if i == 0 else sequence_a[i - 1]
        body = " ".join(cell(str(value)) for value in row)
        rows.append(f"{cell(label)} {body}")
    return "\n".join(rows)
