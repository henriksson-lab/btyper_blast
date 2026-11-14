#!/usr/bin/env python3
"""Record canonical subjects and exact .nsq payloads from pinned makeblastdb v5."""
import argparse
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]


def cases():
    result = []
    for length in range(1, 17):
        result.append(("ACGT" * 4)[:length])
    for base in "ACMGRSVTWYHKDBN-":
        for length in (1, 2, 3, 4, 5, 14, 15, 16, 17, 127, 4094, 4095, 4096, 4097):
            # A record made entirely of parsed '-' gaps is rejected by NCBI's
            # FASTA-to-Bioseq adapter. Keep gaps internal to a valid subject.
            result.append("A" + base * length + "C" if base == "-" else base * length)
    for length in range(1, 130):
        state = length
        sequence = []
        for _ in range(length):
            state = (1664525 * state + 1013904223) & 0xffffffff
            sequence.append("ACMGRSVTWYHKDBN"[(state >> 16) % 15])
        result.append("".join(sequence))
    result.extend(["acgtuX-RYSWKMBDHVN", "AuC", "AXC", "A-C", "A?C", "A.C",
                   "ACGT" + "N" * 4096 + "ACGT", "N" * 4095 + "R" + "N" * 17])
    return result


def read_index(path):
    # Layout documented in upstream writedb_files.cpp: version 5, nucleotide,
    # volume ordinal, three length-prefixed strings, count, LE64 length, BE32 max.
    data = path.read_bytes()
    position = 0

    def integer():
        nonlocal position
        value, = struct.unpack_from(">I", data, position)
        position += 4
        return value

    assert (integer(), integer(), integer()) == (5, 0, 0)
    for _ in range(3):
        length = integer()
        position += length
    count = integer()
    total, = struct.unpack_from("<Q", data, position)
    position += 8
    maximum = integer()
    arrays = []
    for _ in range(3):
        arrays.append([integer() for _ in range(count + 1)])
    assert position == len(data)
    return count, total, maximum, arrays[1], arrays[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    tools = REPO / pin["root"]
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
               NCBI_DONT_USE_LOCAL_CONFIG="1")
    sequences = cases()
    with tempfile.TemporaryDirectory(prefix="btyper-makeblastdb-oracle-") as temporary:
        directory = Path(temporary)
        fasta = directory / "subjects.fa"
        fasta.write_text("".join(f">subject{i}\n{sequence}\n" for i, sequence in enumerate(sequences)))
        database = directory / "reference"
        build = subprocess.run([str(tools / "makeblastdb"), "-in", str(fasta), "-dbtype", "nucl",
                                "-out", str(database)], env=env, capture_output=True, text=True)
        if build.returncode:
            raise RuntimeError(build.stderr)
        count, total, maximum, seq_offsets, amb_offsets = read_index(database.with_suffix(".nin"))
        assert count == len(sequences)
        canonical = subprocess.run([str(tools / "blastdbcmd"), "-db", str(database), "-entry", "all",
                                    "-outfmt", "%s"], env=env, capture_output=True, text=True, check=True).stdout.splitlines()
        assert len(canonical) == count
        data = database.with_suffix(".nsq").read_bytes()
        assert data[0] == 0 and seq_offsets[0] == 1 and seq_offsets[-1] == len(data)
        alphabet = {letter: value for value, letter in enumerate("-ACMGRSVTWYHKDBN")}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w") as output:
            output.write(f"stats {count} {total} {maximum}\n")
            for ordinal, (original, normalized) in enumerate(zip(sequences, canonical)):
                packed = data[seq_offsets[ordinal]:amb_offsets[ordinal]]
                ambiguity = data[amb_offsets[ordinal]:seq_offsets[ordinal + 1]]
                assert len(packed) == len(normalized) // 4 + 1
                encoded = bytes(alphabet[letter] for letter in normalized)
                output.write(f"sequence subject{ordinal} {original.encode().hex()} {encoded.hex()} "
                             f"{packed.hex()} {ambiguity.hex() or '-'}\n")
    print(f"Original makeblastdb v5 payloads written to {args.output}")


if __name__ == "__main__":
    main()
