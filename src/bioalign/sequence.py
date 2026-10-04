"""Turn raw text into a DNA string the rest of the program can trust."""

from __future__ import annotations

VALID_BASES = frozenset("ACGT")


def clean_sequence(raw: str, name: str = "sequence") -> str:
    """Uppercase a DNA string and reject anything that is not A, C, G, or T.

    Spaces and line breaks are ignored so FASTA files and pasted text both work.
    Letters such as N (unknown base) are rejected on purpose: the aligner scores
    one letter against one letter, and an unknown base has no single score.
    """
    if not isinstance(raw, str):
        raise TypeError(f"{name} must be a string of DNA letters")

    sequence = "".join(raw.split()).upper()
    if not sequence:
        raise ValueError(f"{name} is empty")

    bad = sorted({base for base in sequence if base not in VALID_BASES})
    if bad:
        shown = ", ".join(bad)
        raise ValueError(
            f"{name} has letters other than A, C, G, T: {shown}"
        )
    return sequence


def read_fasta(text: str) -> list[tuple[str, str]]:
    """Read FASTA text into a list of (name, sequence) pairs.

    A FASTA record starts with a ">" line (the name) and continues with one or
    more lines of bases. A file with no ">" header is treated as one unnamed
    sequence.
    """
    if not text.strip():
        raise ValueError("FASTA text is empty")

    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    saw_header = False

    def finish() -> None:
        if name is None:
            return
        records.append((name, clean_sequence("".join(chunks), name)))

    for line_number, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            saw_header = True
            finish()
            name = stripped[1:].strip() or f"sequence_{len(records) + 1}"
            chunks = []
            continue
        if not saw_header:
            name = "sequence"
            saw_header = True
        if name is None:
            raise ValueError(f"FASTA line {line_number} has bases before a name")
        chunks.append(stripped)

    finish()
    if not records:
        raise ValueError("FASTA text did not contain a sequence")
    return records
