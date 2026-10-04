"""Short synthetic DNA strings used by the demo and the sample FASTA file.

Nothing here is a real gene. The strings are built so each lesson has an
obvious expected result.
"""

MOTIF = "ATGCGTAAGCTT"

# Same 12-letter motif in the middle. The flanks share no letters,
# so a local alignment should keep only the motif.
FLANKED_A = "TTTTTTTT" + MOTIF + "GGGGGGGG"
FLANKED_B = "CCCCCCCC" + MOTIF + "AAAAAAAA"

# Identical except B has one extra T in the middle.
INDEL_A = "ACGTACGT"
INDEL_B = "ACGTTACGT"

TINY_A = "GATTACA"
TINY_B = "GCATGCT"


def pattern_sequence() -> str:
    """AT-rich block, then a CG island, then an A/T run, then mixed DNA."""
    at_rich = "AT" * 30
    cpg_island = "CG" * 30
    a_tract = "AAAAAATTTTT"
    mixed = "AGCT" * 15
    return at_rich + cpg_island + a_tract + mixed
