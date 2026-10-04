import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bioalign.samples import pattern_sequence
from bioalign.structure import (
    cpg_observed_expected,
    find_a_tracts,
    gc_fraction,
    mean_stacking,
    predict_patterns,
    shannon_entropy,
    train_markov,
)


class MeasurementTests(unittest.TestCase):
    def test_gc_fraction(self):
        self.assertEqual(gc_fraction("GGCC"), 1.0)
        self.assertEqual(gc_fraction("AAAA"), 0.0)

    def test_cpg_enrichment(self):
        self.assertGreater(cpg_observed_expected("CG" * 20), 1.5)
        self.assertEqual(cpg_observed_expected("AT" * 20), 0.0)

    def test_gc_stacks_are_more_stable_than_at_stacks(self):
        self.assertLess(mean_stacking("CG" * 20), mean_stacking("AT" * 20))
        self.assertLess(mean_stacking("CG" * 20), -1.70)
        self.assertGreater(mean_stacking("AT" * 20), -1.10)

    def test_entropy_endpoints(self):
        self.assertEqual(shannon_entropy("AAAA"), 0.0)
        self.assertAlmostEqual(shannon_entropy("ACGT"), 2.0)

    def test_markov_rows_are_probabilities_and_favor_the_common_pair(self):
        model = train_markov("A" * 20 + "C")
        for previous in "ACGT":
            total = sum(model.transitions[previous].values())
            self.assertAlmostEqual(total, 1.0)
        self.assertAlmostEqual(sum(model.base_probability.values()), 1.0)
        self.assertGreater(model.transitions["A"]["A"], 0.8)

    def test_a_tract_coordinates(self):
        sequence = "GGAAAAAATT"
        tracts = find_a_tracts(sequence)
        self.assertEqual([(tract.start, tract.end) for tract in tracts], [(2, 8)])
        self.assertEqual(sequence[2:8], "AAAAAA")


class PatternCallTests(unittest.TestCase):
    def test_composite_sequence_calls_the_planted_patterns(self):
        sequence = pattern_sequence()
        report = predict_patterns(sequence, window=40, step=20)
        labels = {region.label for region in report.regions}
        self.assertIn("cpg_like", labels)
        self.assertIn("stable_duplex", labels)
        self.assertIn("at_flexible", labels)
        self.assertIn("a_tract", labels)

        tracts = [region for region in report.regions if region.label == "a_tract"]
        pieces = {sequence[region.start : region.end] for region in tracts}
        self.assertIn("AAAAAA", pieces)
        self.assertIn("TTTTT", pieces)

    def test_pure_cg_window_is_cpg_like_and_stable(self):
        report = predict_patterns("CG" * 30, window=40, step=40)
        labels = {label for window in report.windows for label in window.labels}
        self.assertIn("cpg_like", labels)
        self.assertIn("stable_duplex", labels)

    def test_fasta_pattern_record_matches_the_builder(self):
        from bioalign.sequence import read_fasta

        records = read_fasta(Path("data/sample_sequences.fasta").read_text())
        self.assertEqual(records[2][1], pattern_sequence())


if __name__ == "__main__":
    unittest.main()
