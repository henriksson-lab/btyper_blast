#!/usr/bin/env python3
"""Execute unchanged BLAST nucleotide word-finder kernels."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

REPO = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_source.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    root = REPO / pin["root"]
    core = root / "src/algo/blast/core"
    with tempfile.TemporaryDirectory(prefix="btyper-word-finder-oracle-") as temporary:
        executable = Path(temporary) / "ungapped"
        command = shlex.split(os.environ.get("CC", "cc")) + [
            "-std=c11", "-D_DEFAULT_SOURCE", "-O2", "-fwrapv", "-ffp-contract=off",
            "-ffunction-sections", "-fdata-sections", "-I" + str(REPO / "tests/reference/include"),
            "-I" + str(root / "include"), "-I" + str(core), str(REPO / "tests/reference/word_finder.c"),
        ]
        command += [str(core / (name + ".c")) for name in
                    ("blast_diagnostics", "blast_stat", "ncbi_math", "blast_program", "blast_encoding", "blast_query_info",
                     "blast_filter", "blast_options", "blast_seqsrc", "lookup_util")]
        command += ["-Wl,-dead_strip" if sys.platform == "darwin" else "-Wl,--gc-sections", "-lm", "-o", str(executable)]
        subprocess.run(command, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    print(f"Original word-finder reference written to {args.output}")


if __name__ == "__main__":
    main()
