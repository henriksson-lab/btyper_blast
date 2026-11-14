#!/usr/bin/env python3
"""Capture default query encodings and DUST masks from official BLAST archives.

The archive request LCaseMask options retain whole-query masks and duplicate
strand entries; the reported result masks can omit whole-query intervals.
"""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from source_ranges import definition

REPO = Path(__file__).resolve().parents[1]


def cases():
    sequences = []
    for motif in ("A", "C", "G", "T", "N", "AT", "ACG", "ACGT", "RYSWKMBDHVN"):
        for length in (4, 8, 31, 62, 63, 64, 65, 127, 512):
            sequence = (motif * (length // len(motif) + 1))[:length]
            sequences.append(sequence)
            sequences.append("ACGTGCACTG" + sequence + "TGCACTGAGC")
    state = 0x2b3d15a9
    for i in range(40):
        sequence = []
        for _ in range(100 + i * 5):
            state = (state * 1664525 + 1013904223) & 0xffffffff
            sequence.append("ACGTNRYSWKMBVDH"[(state >> 16) % 15])
        sequence = "".join(sequence)
        sequences.append(sequence[:30] + "A" * 80 + sequence[30:])
    sequences.extend(["acgtu" + "a" * 80 + "RYSWKMBDHVN", "ACGT" + "-" * 70 + "A" * 80 + "CAGT",
                      "NNNACGT" + "A" * 100 + "N" * 100 + "C" * 100,
                      "N" * 4096, "A" * 4096])
    return sequences


def parse_archive(archive):
    sequences = {}
    for match in re.finditer(r"^          seq\s*(?=\{)", archive, re.M):
        start, end = definition(archive, r"^          seq\s*", match.start())
        block = archive[start:end]
        ordinal = int(re.search(r'id\s*\{\s*local str "Query_(\d+)"', block)[1]) - 1
        length = int(re.search(r"\blength (\d+)", block)[1])
        coding = re.search(r"\bseq-data (ncbi2na|ncbi4na)\s*'([0-9A-Fa-f\s]+)'H", block)
        if not coding:
            raise ValueError("Unexpected raw nucleotide archive encoding")
        packed = bytes.fromhex(coding[2])
        if coding[1] == "ncbi2na":
            decoded = bytes(1 << ((packed[i // 4] >> (6 - 2 * (i % 4))) & 3) for i in range(length))
        else:
            decoded = bytes((packed[i // 2] >> (4 if i % 2 == 0 else 0)) & 15 for i in range(length))
        assert ordinal not in sequences
        sequences[ordinal] = decoded
    masks = defaultdict(set)
    for match in re.finditer(r'name "LCaseMask"', archive):
        start, end = definition(archive, r"\bvalue query-mask\s*", match.end())
        for left, right, ordinal in re.findall(r'from (\d+),\s*to (\d+),\s*id local str "Query_(\d+)"', archive[start:end]):
            masks[int(ordinal) - 1].add((int(left), int(right)))
    return sequences, {ordinal: sorted(intervals) for ordinal, intervals in masks.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    tools = REPO / pin["root"]
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1", NCBI_DONT_USE_LOCAL_CONFIG="1")
    env.pop("BLASTINPUT_GEN_DELTA_SEQ", None)
    originals = cases()
    with tempfile.TemporaryDirectory(prefix="btyper-query-oracle-") as temporary:
        directory = Path(temporary)
        query_path = directory / "queries.fa"
        query_path.write_text("".join(f">case{i}\n{sequence}\n" for i, sequence in enumerate(originals)))
        subject_path = directory / "subject.fa"
        # A modest database is enough to exercise full setup; hits are not used.
        state = 45179
        subject = []
        for _ in range(3000):
            state = (state * 1664525 + 1013904223) & 0xffffffff
            subject.append("ACGT"[state >> 30])
        subject_path.write_text(">subject\n" + "".join(subject) + "\n")
        db = directory / "subject"
        subprocess.run([str(tools / "makeblastdb"), "-in", str(subject_path), "-dbtype", "nucl", "-out", str(db)],
                       env=env, capture_output=True, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w") as output:
            for task, task_args in (("megablast", []), ("blastn", ["-task", "blastn"])):
                search = subprocess.run([str(tools / "blastn"), "-query", str(query_path), "-db", str(db),
                                         "-num_threads", "1", "-outfmt", "11"] + task_args,
                                        env=env, capture_output=True, text=True, check=True)
                assert 'ident "2.17.0+"' in search.stdout
                assert re.search(r'name "MaskAtHash",\s*value boolean TRUE', search.stdout)
                sequences, masks = parse_archive(search.stdout)
                assert set(sequences) == set(range(len(originals)))
                assert set(masks) <= set(sequences)
                for ordinal, original in enumerate(originals):
                    encoded_masks = ",".join(f"{start}:{stop}" for start, stop in masks.get(ordinal, [])) or "-"
                    output.write(f"query {task} case{ordinal} {original.encode().hex()} {sequences[ordinal].hex()} {encoded_masks}\n")
    print(f"Original BLAST query reference written to {args.output}")


if __name__ == "__main__":
    main()
