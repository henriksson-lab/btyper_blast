#!/usr/bin/env python3
"""Compile unchanged pinned C lookup kernels and record their defined outputs."""
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
    parser.add_argument("--kind", choices=("lookup", "ordinary_lookup", "mb_lookup", "seed_lookup"), default="lookup")
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_source.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    root = REPO / pin["root"]
    core = root / "src/algo/blast/core"
    with tempfile.TemporaryDirectory(prefix="btyper-lookup-oracle-") as temporary:
        executable = Path(temporary) / "lookup"
        command = shlex.split(os.environ.get("CC", "cc")) + [
            "-std=c11", "-D_DEFAULT_SOURCE", "-O2", "-fwrapv",
            "-ffunction-sections", "-fdata-sections",
            "-I" + str(REPO / "tests/reference/include"),
            "-I" + str(root / "include"), "-I" + str(core),
            str(REPO / "tests/reference" / (args.kind + ".c")), str(core / "blast_program.c"),
            str(core / "blast_filter.c"),
            str(core / "blast_options.c"),
            str(core / "blast_seqsrc.c"),
            "-Wl,-dead_strip" if sys.platform == "darwin" else "-Wl,--gc-sections",
            "-lm", "-o", str(executable),
        ]
        subprocess.run(command, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    print(f"Original C lookup reference written to {args.output}")


if __name__ == "__main__":
    main()
