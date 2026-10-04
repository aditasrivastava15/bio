"""DNA alignment and statistical pattern prediction.

Two jobs live here:

1. Line up two DNA strings and measure how similar they are.
2. Scan one DNA string for stretches that look structurally unusual.
"""

from bioalign.needleman_wunsch import needleman_wunsch
from bioalign.scoring import Scoring
from bioalign.sequence import clean_sequence, read_fasta
from bioalign.smith_waterman import smith_waterman
from bioalign.structure import predict_patterns

__all__ = [
    "Scoring",
    "clean_sequence",
    "needleman_wunsch",
    "predict_patterns",
    "read_fasta",
    "smith_waterman",
]
