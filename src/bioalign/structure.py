"""Scan one DNA sequence for structural patterns using counts, not a 3D model.

Four separate measurements are combined:

1. Composition. What fraction of the letters are G or C? GC-rich DNA holds
   the two strands together more tightly.
2. A first-order Markov chain. This is a table of "if the current letter is X,
   how often is the next letter Y?" Trained on the whole sequence, the table
   is a small statistical model of that sequence's habits.
3. A sliding window. The same measurements are repeated on short stretches so
   a local pattern is not diluted by the rest of the sequence.
4. Named calls, using published rules of thumb:
   - CpG-like: GC-rich and CG pairs are more common than chance
     (Gardiner-Garden and Frommer, 1987, scaled down to the window size)
   - stable duplex / AT-flexible: average nearest-neighbor stacking energy
     (SantaLucia, 1998). More negative means a more stable helix.
   - A-tract: a run of A or a run of T. Runs of about four or more tend to bend DNA.
   - low complexity: the stretch is dominated by one letter
   - unlike the rest: the window is a poor fit to the Markov model of the whole sequence
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass

from bioalign.sequence import clean_sequence

# SantaLucia 1998 nearest-neighbor ΔG° at 37°C, kcal/mol, for one DNA step.
# More negative means the two bases stack into a more stable duplex.
# The reverse-complement step has the same energy, so AA and TT share a value.
STACKING = {
    "AA": -1.00,
    "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "TG": -1.45,
    "GT": -1.44,
    "AC": -1.44,
    "CT": -1.28,
    "AG": -1.28,
    "GA": -1.30,
    "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
    "CC": -1.84,
}

CPG_GC_MIN = 0.50
CPG_OE_MIN = 0.60
STABLE_STACKING = -1.70
FLEXIBLE_STACKING = -1.10
LOW_ENTROPY = 0.50
A_TRACT_MIN = 4

LABELS = {
    "cpg_like": "CpG-like",
    "stable_duplex": "stable duplex",
    "at_flexible": "AT-flexible",
    "low_complexity": "low complexity",
    "a_tract": "A-tract",
    "unlike_background": "unlike the rest",
}


@dataclass(frozen=True)
class Composition:
    """Letter counts for one sequence."""

    length: int
    counts: dict[str, int]

    def fraction(self, base: str) -> float:
        return self.counts[base] / self.length

    @property
    def gc(self) -> float:
        return (self.counts["G"] + self.counts["C"]) / self.length


@dataclass(frozen=True)
class MarkovModel:
    """P(letter) and P(next | current), with a little pseudocount smoothing.

    Pseudocounts stop a never-seen pair from getting probability 0, which
    would make the log-probability explode to negative infinity.
    """

    base_probability: dict[str, float]
    transitions: dict[str, dict[str, float]]

    def log2_probability(self, sequence: str) -> float:
        """Log2 of the probability of a sequence under this model.

        The first letter uses the base frequency. Every later letter uses the
        transition from the letter before it. Dividing by the length later
        makes long and short windows comparable.
        """
        total = math.log2(self.base_probability[sequence[0]])
        for previous, nxt in zip(sequence, sequence[1:]):
            total += math.log2(self.transitions[previous][nxt])
        return total


@dataclass(frozen=True)
class WindowProfile:
    """Measurements for one stretch of the sequence. ``end`` is exclusive."""

    start: int
    end: int
    gc: float
    cpg_observed_expected: float
    stacking: float
    entropy: float
    log2_per_base: float
    labels: tuple[str, ...]


@dataclass(frozen=True)
class Region:
    """A merged run of windows (or an exact A/T run) that share one label."""

    start: int
    end: int
    label: str

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PatternReport:
    sequence: str
    composition: Composition
    model: MarkovModel
    windows: tuple[WindowProfile, ...]
    regions: tuple[Region, ...]


def composition_of(sequence: str) -> Composition:
    sequence = clean_sequence(sequence)
    return Composition(
        length=len(sequence),
        counts={base: sequence.count(base) for base in "ACGT"},
    )


def train_markov(sequence: str, pseudocount: float = 0.5) -> MarkovModel:
    """Count letters and neighboring pairs, then turn the counts into probabilities."""
    if pseudocount <= 0:
        raise ValueError("pseudocount must be positive")
    sequence = clean_sequence(sequence)
    bases = "ACGT"
    base_counts = {base: pseudocount for base in bases}
    pair_counts = {base: {nxt: pseudocount for nxt in bases} for base in bases}

    for base in sequence:
        base_counts[base] += 1
    for previous, nxt in zip(sequence, sequence[1:]):
        pair_counts[previous][nxt] += 1

    base_total = sum(base_counts.values())
    base_probability = {base: base_counts[base] / base_total for base in bases}
    transitions = {
        previous: {
            nxt: pair_counts[previous][nxt] / sum(pair_counts[previous].values())
            for nxt in bases
        }
        for previous in bases
    }
    return MarkovModel(base_probability=base_probability, transitions=transitions)


def gc_fraction(sequence: str) -> float:
    if not sequence:
        return 0.0
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def cpg_observed_expected(sequence: str) -> float:
    """How many CG pairs there are, divided by how many chance would predict.

    Chance here means "C's and G's placed independently." The formula is
    (count of CG) * (length) / (count of C * count of G). A value near 1
    means CG appears about as often as the C and G counts suggest. Above 1
    means CG is enriched.
    """
    length = len(sequence)
    if length < 2:
        return 0.0
    cytosines = sequence.count("C")
    guanines = sequence.count("G")
    if cytosines == 0 or guanines == 0:
        return 0.0
    cg_pairs = sum(
        1 for index in range(length - 1) if sequence[index : index + 2] == "CG"
    )
    return (cg_pairs * length) / (cytosines * guanines)


def mean_stacking(sequence: str) -> float:
    """Average nearest-neighbor stacking energy. More negative means more stable."""
    if len(sequence) < 2:
        return 0.0
    total = sum(STACKING[sequence[index : index + 2]] for index in range(len(sequence) - 1))
    return total / (len(sequence) - 1)


def shannon_entropy(sequence: str) -> float:
    """How mixed the letters are, in bits. 0 = one letter only. 2 = all four equal."""
    if not sequence:
        return 0.0
    entropy = 0.0
    for base in "ACGT":
        probability = sequence.count(base) / len(sequence)
        if probability > 0:
            entropy -= probability * math.log2(probability)
    return entropy


def find_a_tracts(sequence: str, min_run: int = A_TRACT_MIN) -> tuple[Region, ...]:
    """Find exact runs of the same A, or the same T, at least ``min_run`` long."""
    if min_run < 1:
        raise ValueError("min_run must be at least 1")
    sequence = clean_sequence(sequence)
    regions: list[Region] = []
    index = 0
    while index < len(sequence):
        base = sequence[index]
        if base not in "AT":
            index += 1
            continue
        end = index + 1
        while end < len(sequence) and sequence[end] == base:
            end += 1
        if end - index >= min_run:
            regions.append(Region(start=index, end=end, label="a_tract"))
        index = end
    return tuple(regions)


def predict_patterns(
    sequence: str,
    window: int = 40,
    step: int = 10,
    pseudocount: float = 0.5,
) -> PatternReport:
    """Train a Markov model on the sequence and label unusual windows."""
    if window < 2:
        raise ValueError("window must be at least 2")
    if step < 1:
        raise ValueError("step must be at least 1")

    sequence = clean_sequence(sequence)
    model = train_markov(sequence, pseudocount=pseudocount)
    spans = _windows(len(sequence), window, step)

    measured: list[dict] = []
    for start, end in spans:
        chunk = sequence[start:end]
        labels: list[str] = []
        gc = gc_fraction(chunk)
        observed_expected = cpg_observed_expected(chunk)
        stacking = mean_stacking(chunk)
        entropy = shannon_entropy(chunk)
        if len(chunk) >= 20 and gc >= CPG_GC_MIN and observed_expected >= CPG_OE_MIN:
            labels.append("cpg_like")
        if stacking <= STABLE_STACKING:
            labels.append("stable_duplex")
        elif stacking >= FLEXIBLE_STACKING:
            labels.append("at_flexible")
        if entropy <= LOW_ENTROPY:
            labels.append("low_complexity")
        measured.append(
            {
                "start": start,
                "end": end,
                "gc": gc,
                "cpg_observed_expected": observed_expected,
                "stacking": stacking,
                "entropy": entropy,
                "log2_per_base": model.log2_probability(chunk) / len(chunk),
                "labels": labels,
            }
        )

    _mark_unlike_background(measured)

    windows = tuple(
        WindowProfile(
            start=item["start"],
            end=item["end"],
            gc=item["gc"],
            cpg_observed_expected=item["cpg_observed_expected"],
            stacking=item["stacking"],
            entropy=item["entropy"],
            log2_per_base=item["log2_per_base"],
            labels=tuple(item["labels"]),
        )
        for item in measured
    )
    regions = _regions_from(windows, find_a_tracts(sequence))
    return PatternReport(
        sequence=sequence,
        composition=composition_of(sequence),
        model=model,
        windows=windows,
        regions=regions,
    )


def format_report(report: PatternReport) -> str:
    """Plain-text summary of composition, the Markov table, regions, and windows."""
    composition = report.composition
    counts = "   ".join(
        f"{base} {composition.counts[base]:<4} {composition.fraction(base):6.1%}"
        for base in "ACGT"
    )
    lines = [
        f"length {composition.length}    GC content {composition.gc:.1%}",
        counts,
        "",
        "Markov model: probability that the column letter follows the row letter",
        _format_transitions(report.model),
        "",
        "Called regions",
    ]
    if not report.regions:
        lines.append("  (none)")
    for region in report.regions:
        preview = report.sequence[region.start : region.end]
        if len(preview) > 24:
            preview = preview[:24] + "..."
        lines.append(
            f"  {region.start:>4}-{region.end:<4}  {LABELS[region.label]:<16}  "
            f"length {region.length:<4}  {preview}"
        )

    lines.extend(["", "Windows", _window_header()])
    for window in report.windows:
        label_text = ",".join(LABELS[label] for label in window.labels) or "-"
        lines.append(
            f"  {window.start:>4}-{window.end:<4}  "
            f"GC {window.gc:5.1%}  "
            f"CpG {window.cpg_observed_expected:5.2f}  "
            f"stack {window.stacking:6.2f}  "
            f"H {window.entropy:4.2f}  "
            f"bits {window.log2_per_base:6.2f}  "
            f"{label_text}"
        )
    lines.extend(
        [
            "",
            "stack is average stacking energy in kcal/mol (more negative = more stable helix).",
            "H is Shannon entropy in bits (0 = one repeated letter, 2 = equal A/C/G/T).",
            "bits is the Markov log2-probability per base (closer to 0 = more expected).",
            "CpG is observed/expected CG pairs (above 0.60 with GC >= 50% calls CpG-like).",
        ]
    )
    return "\n".join(lines)


def _windows(length: int, window: int, step: int) -> list[tuple[int, int]]:
    if length <= window:
        return [(0, length)]
    spans = []
    start = 0
    while start + window <= length:
        spans.append((start, start + window))
        start += step
    if spans[-1][1] < length:
        spans.append((length - window, length))
    return spans


def _mark_unlike_background(measured: list[dict]) -> None:
    """Flag windows whose per-base probability sits more than one stdev below the mean.

    With a single window the sequence is the background, so nothing is flagged.
    """
    if len(measured) < 2:
        return
    scores = [item["log2_per_base"] for item in measured]
    cutoff = statistics.mean(scores) - statistics.pstdev(scores)
    if statistics.pstdev(scores) < 1e-9:
        return
    for item in measured:
        if item["log2_per_base"] < cutoff:
            item["labels"].append("unlike_background")


def _regions_from(
    windows: tuple[WindowProfile, ...],
    a_tracts: tuple[Region, ...],
) -> tuple[Region, ...]:
    by_label: dict[str, list[tuple[int, int]]] = {}
    for window in windows:
        for label in window.labels:
            by_label.setdefault(label, []).append((window.start, window.end))

    regions: list[Region] = []
    for label, spans in by_label.items():
        for start, end in _merge(spans):
            regions.append(Region(start=start, end=end, label=label))
    # A-tracts stay as exact runs. An A-run touching a T-run is two bends, not one.
    regions.extend(a_tracts)
    regions.sort(key=lambda region: (region.start, region.end, region.label))
    return tuple(regions)


def _merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    ordered = sorted(spans)
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def _format_transitions(model: MarkovModel) -> str:
    header = "       " + " ".join(f"{base:>6}" for base in "ACGT")
    rows = [header]
    for previous in "ACGT":
        cells = " ".join(
            f"{model.transitions[previous][nxt]:6.2f}" for nxt in "ACGT"
        )
        rows.append(f"    {previous}  {cells}")
    return "\n".join(rows)


def _window_header() -> str:
    return (
        "  start-end     GC      CpG     stack      H       bits   labels"
    )
