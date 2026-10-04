import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bioalign.alignment import format_alignment
from bioalign.needleman_wunsch import global_score_matrix, needleman_wunsch
from bioalign.samples import FLANKED_A, FLANKED_B, MOTIF, TINY_A, TINY_B
from bioalign.scoring import Scoring
from bioalign.sequence import clean_sequence, read_fasta
from bioalign.smith_waterman import smith_waterman


class SequenceTests(unittest.TestCase):
    def test_clean_accepts_lowercase_and_whitespace(self):
        self.assertEqual(clean_sequence("  at gc\n"), "ATGC")

    def test_clean_rejects_unknown_letters(self):
        with self.assertRaises(ValueError):
            clean_sequence("ATGCN")

    def test_fasta_round_trip_shape(self):
        text = Path("data/sample_sequences.fasta").read_text()
        records = read_fasta(text)
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0][0], "flanked_a")
        self.assertEqual(records[0][1], FLANKED_A)
        self.assertEqual(records[1][1], FLANKED_B)


class GlobalAlignmentTests(unittest.TestCase):
    def test_identical_sequences_are_a_perfect_match(self):
        result = needleman_wunsch("ACGT", "ACGT")
        self.assertEqual(result.aligned_a, "ACGT")
        self.assertEqual(result.aligned_b, "ACGT")
        self.assertEqual(result.score, 8)
        self.assertEqual(result.identity, 1.0)

    def test_all_mismatches(self):
        result = needleman_wunsch("AAAA", "CCCC")
        self.assertEqual(result.score, -4)
        self.assertEqual(result.matches, 0)

    def test_one_gap_for_a_missing_letter(self):
        result = needleman_wunsch("ACGT", "AGT")
        self.assertEqual(result.score, 4)
        self.assertEqual(result.matches, 3)
        self.assertEqual(result.gaps, 1)
        self.assertIn("-", result.aligned_a + result.aligned_b)

    def test_hand_worked_pair(self):
        result = needleman_wunsch("GA", "GC")
        self.assertEqual(result.aligned_a, "GA")
        self.assertEqual(result.aligned_b, "GC")
        self.assertEqual(result.score, 1)

    def test_border_of_the_grid_is_the_gap_cost(self):
        matrix = global_score_matrix("AC", "AGT")
        self.assertEqual(matrix[0][0], 0)
        self.assertEqual(matrix[1][0], -2)
        self.assertEqual(matrix[2][0], -4)
        self.assertEqual(matrix[0][3], -6)
        self.assertEqual(len(matrix), 3)
        self.assertEqual(len(matrix[0]), 4)

    def test_flanked_motif_keeps_the_mismatched_ends(self):
        result = needleman_wunsch(FLANKED_A, FLANKED_B)
        self.assertEqual(result.matches, len(MOTIF))
        self.assertLess(result.identity, 0.5)
        self.assertIn(MOTIF, result.aligned_a.replace("-", ""))


class LocalAlignmentTests(unittest.TestCase):
    def test_no_shared_letters_scores_zero(self):
        result = smith_waterman("AAAA", "CCCC")
        self.assertEqual(result.score, 0)
        self.assertEqual(result.aligned_a, "")
        self.assertEqual(result.identity, 0.0)

    def test_finds_an_embedded_motif(self):
        result = smith_waterman("TTTT" + "ACGTACGT" + "GGGG", "ACGTACGT")
        self.assertEqual(result.aligned_a, "ACGTACGT")
        self.assertEqual(result.aligned_b, "ACGTACGT")
        self.assertEqual(result.score, 16)

    def test_single_match_beats_a_mismatch_extension(self):
        result = smith_waterman("GA", "GC")
        self.assertEqual(result.score, 2)
        self.assertEqual(result.aligned_a, "G")
        self.assertEqual(result.aligned_b, "G")

    def test_flanked_motif_is_the_local_hit(self):
        result = smith_waterman(FLANKED_A, FLANKED_B)
        self.assertEqual(result.aligned_a, MOTIF)
        self.assertEqual(result.aligned_b, MOTIF)
        self.assertEqual(result.identity, 1.0)
        self.assertGreater(result.score, needleman_wunsch(FLANKED_A, FLANKED_B).score)

    def test_tiny_demo_pair_is_stable(self):
        result = needleman_wunsch(TINY_A, TINY_B)
        self.assertEqual(len(result.aligned_a), len(result.aligned_b))
        self.assertIn("GATTACA", result.aligned_a.replace("-", ""))
        text = format_alignment(result)
        self.assertIn("global alignment", text)


class ScoringTests(unittest.TestCase):
    def test_positive_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            Scoring(gap=2)


if __name__ == "__main__":
    unittest.main()
