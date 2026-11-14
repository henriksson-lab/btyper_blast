#!/usr/bin/env python3
"""Generate deterministic synthetic search fixtures (not a performance dataset)."""
import argparse
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPO / "tests/data/smoke")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    gc_source = REPO / "vendor/ncbi-blast-2.17.0+-src/c++/src/objects/seqfeat/gc.prt"
    code = re.search(r'id\s+1\s*,\s*ncbieaa\s+"([A-Z*]+)"', gc_source.read_text()).group(1)
    codons = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
    choices = {}
    for codon, aa in zip(codons, code):
        choices.setdefault(aa, []).append(codon)
    state = 0x3BD91E07

    def letters(alphabet, length):
        nonlocal state
        result = []
        for _ in range(length):
            state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
            result.append(alphabet[(state >> 16) % len(alphabet)])
        return "".join(result)

    def reverse(sequence):
        return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

    def fasta(name, records):
        with (args.output / name).open("w") as output:
            for identifier, sequence in records:
                output.write(">" + identifier + "\n")
                for i in range(0, len(sequence), 80):
                    output.write(sequence[i:i + 80] + "\n")

    peptide = letters("ACDEFGHIKLMNPQRSTVWY", 80)
    cds = "".join(choices[aa][i % len(choices[aa])] for i, aa in enumerate(peptide))
    nucleotide_query = letters("ACGT", 180)
    subjects = []
    for i in range(6):
        sequence = letters("ACGT", 4096)
        positive_cds = cds if i < 3 else reverse(cds)
        # Positive frames +1,+2,+3 and negative frames -1,-2,-3.
        position = 600 + i if i < 3 else 4096 - 240 - (600 + i - 3)
        sequence = sequence[:position] + positive_cds + sequence[position + len(cds):]
        insert = nucleotide_query if i % 2 == 0 else reverse(nucleotide_query)
        if i == 2:
            insert = insert[:70] + "GAC" + insert[70:]
        if i == 4:
            insert = insert[:40] + "NRY" + insert[43:]
        sequence = sequence[:1800] + insert + sequence[1800 + len(insert):]
        subjects.append((f"subject{i}", sequence))
    fasta("subjects.fna", subjects)
    fasta("nucleotide-queries.fna", [
        ("nucleotide_match", nucleotide_query),
        ("nucleotide_no_hit", letters("ACGT", 150)),
        ("low_complexity", "AT" * 80),
    ])
    fasta("protein-queries.faa", [
        ("protein_match", peptide),
        ("protein_no_hit", letters("ACDEFGHIKLMNPQRSTVWY", 70)),
    ])


if __name__ == "__main__":
    main()
