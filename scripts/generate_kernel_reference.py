#!/usr/bin/env python3
"""Build unchanged pinned C kernels and emit differential fixtures (Linux/macOS)."""
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
    with tempfile.TemporaryDirectory(prefix="btyper-blast-oracle-") as temporary:
        executable = Path(temporary) / "kernels"
        command = shlex.split(os.environ.get("CC", "cc")) + [
            "-std=c11", "-D_DEFAULT_SOURCE", "-O2", "-ffp-contract=off",
            "-ffunction-sections", "-fdata-sections",
            "-I" + str(REPO / "tests/reference/include"),
            "-I" + str(root / "include"), "-I" + str(core),
            str(REPO / "tests/reference/kernels.c"),
        ]
        command += [str(core / (name + ".c")) for name in
                    ("blast_encoding", "blast_program", "gapinfo", "ncbi_math")]
        command += ["-Wl,-dead_strip" if sys.platform == "darwin" else "-Wl,--gc-sections",
                    "-lm", "-o", str(executable)]
        subprocess.run(command, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    print(f"C reference fixtures written to {args.output}")


if __name__ == "__main__":
    main()
