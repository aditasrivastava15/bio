# Bioinformatics Sequence Analysis

This program compares DNA sequences and looks for simple structural patterns in them.

Give it two DNA strings and it lines them up so you can see where they match, mismatch, or have an extra/missing letter. It does this two ways:

- **Needleman–Wunsch (global):** uses the whole sequence, even the ends that do not match
- **Smith–Waterman (local):** keeps only the best matching stretch and ignores the rest

Give it one DNA string and it scans that string with a sliding window. It does not fold DNA into 3D. It uses counts and a small probability model to flag stretches that look GC-rich, unusually stable, AT-flexible, or like an A-tract.

Input is DNA letters (`A`, `C`, `G`, `T`), typed in or read from a FASTA file. Output is text: an alignment picture with a score, or a table of labeled regions.

```bash
python3 demo.py
python3 demo.py align GATTACA GCATGCT
python3 demo.py patterns ATGCGTACGT
python3 demo.py fasta data/sample_sequences.fasta
```
