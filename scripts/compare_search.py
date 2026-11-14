#!/usr/bin/env python3
"""Compare native in-memory search/display with the existing pinned BLAST smoke.

All non-floating fields and hit order must match exactly. Native raw floats are
compared after the original tabular precision rules (not an arbitrary tolerance).
This small fixture is an integration check, not a complete parity or speed proof.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess

from reference_search import FIELDS

REPO = Path(__file__).resolve().parents[1]
PROGRAMS = {"megablast": "blastn_default_megablast", "blastn": "blastn_traditional",
            "tblastn": "tblastn_default"}


def displayed_number(field, value):
    """BLAST 2.17.0+ tabular.cpp SetScores and GetScoreString formatting."""
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"invalid native {field}: {value}")
    if field == "bitscore":
        if value > 99999:
            return float(f"{value:.3e}")
        if value > 99.9:
            return float(int(value))
        return float(f"{value:.1f}")
    if value < 1e-180:
        return 0.0
    if value < 0.0009:
        return float(f"{value:.2e}")
    if value < 0.1:
        return float(f"{value:.3f}")
    if value < 1:
        return float(f"{value:.2f}")
    if value < 10:
        return float(f"{value:.1f}")
    return float(f"{value:.0f}")


def compare(native, expected, program):
    if len(native) != len(expected):
        raise ValueError(f"{program}: hit count {len(native)} != original {len(expected)}")
    differences = []
    for index, (actual, original) in enumerate(zip(native, expected)):
        if len(actual) != len(FIELDS):
            raise ValueError(f"{program} row {index}: {len(actual)} fields, expected {len(FIELDS)}")
        for field, value in zip(FIELDS, actual):
            same = (displayed_number(field, float(value)) == float(original[field])
                    if field in ("bitscore", "evalue") else value == str(original[field]))
            if not same:
                differences.append(f"row {index} {field}: native {value!r}, original {original[field]!r}")
    if differences:
        raise ValueError(program + ":\n" + "\n".join(differences[:30]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--program", action="append", choices=PROGRAMS,
                        help="Repeat to select programs; defaults to all three presets")
    args = parser.parse_args()
    binary = args.native.resolve()
    if not binary.exists() and binary.with_suffix(".exe").exists():
        binary = binary.with_suffix(".exe")
    fixture = json.loads((REPO / "tests/fixtures/search-smoke-2.17.0.json").read_text())
    version = json.loads((REPO / "upstream/blast-source.json").read_text())["version"]
    if fixture["reference_version"] != version:
        raise ValueError("fixture and translated source versions differ")
    for program in args.program or PROGRAMS:
        query = "protein-queries.faa" if program == "tblastn" else "nucleotide-queries.fna"
        result = subprocess.run([str(binary), program,
                                 str(REPO / "tests/data/smoke/subjects.fna"),
                                 str(REPO / "tests/data/smoke" / query)],
                                text=True, capture_output=True, check=True)
        rows = [line.split("\t") for line in result.stdout.splitlines()]
        compare(rows, fixture["searches"][PROGRAMS[program]]["hits"], program)
        print(f"{program}: {len(rows)} hits match BLAST+ {version}, including order, alignments and printed statistics")


if __name__ == "__main__":
    main()
