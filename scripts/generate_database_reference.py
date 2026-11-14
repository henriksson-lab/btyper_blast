#!/usr/bin/env python3
"""Compile unchanged selected C++ writer/RNG definitions into a reference driver."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

from source_ranges import definition

REPO = Path(__file__).resolve().parents[1]
CPP = "src/objtools/blast/seqdb_writer/writedb_convert.cpp"
RNG_CPP = "src/util/random_gen.cpp"
RNG_HPP = "include/util/random_gen.hpp"


def fragments(root):
    selections = [
        (RNG_CPP, r"^CRandom::CRandom\(TValue seed\)", False),
        (RNG_CPP, r"^void CRandom::SetSeed\(TValue seed\)", False),
        (RNG_HPP, r"^inline CRandom::TValue CRandom::x_GetRand32Bits\(void\)", False),
        (RNG_HPP, r"^inline CRandom::TValue CRandom::GetRand\(void\)", False),
        ("include/objtools/blast/seqdb_writer/writedb_general.hpp", r"^inline int s_DivideRoundUp\(", False),
        ("include/objtools/blast/seqdb_writer/writedb_convert.hpp", r"^inline void s_AppendInt4\(", False),
        (CPP, r"^class CAmbiguousRegion\s*", True),
        (CPP, r"^class CAmbigDataBuilder\s*", True),
        (CPP, r"^inline vector<unsigned char> s_BuildNa4ToNa2Table\(\)", False),
        (CPP, r"^void WriteDB_Ncbi4naToBinary\(const char \* ncbi4na,", False),
    ]
    result = []
    for path, pattern, semicolon in selections:
        text = (root / path).read_text()
        start, end = definition(text, pattern, semicolon=semicolon)
        data = text[start:end]
        result.append({"source": path, "byte_start": len(text[:start].encode()),
                       "byte_end": len(text[:end].encode()), "sha256": hashlib.sha256(data.encode()).hexdigest(),
                       "text": data})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fragments", type=Path, help="Write extracted source provenance")
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_source.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    root = REPO / pin["root"]
    selected = fragments(root)
    shim = (REPO / "tests/reference/database_shim.hpp").read_text()
    # The state offset comes from its original declaration, not a rewritten value.
    offset = next(line for line in (root / RNG_CPP).read_text().splitlines()
                  if line.startswith("static const size_t kStateOffset ="))
    with tempfile.TemporaryDirectory(prefix="btyper-database-oracle-") as temporary:
        directory = Path(temporary)
        (directory / "reference_database.hpp").write_text(
            shim + "\n" + offset + "\n\n" + "\n\n".join(entry["text"] for entry in selected))
        executable = directory / "database"
        command = shlex.split(os.environ.get("CXX", "c++")) + [
            "-std=c++17", "-O2", "-I" + str(directory),
            str(REPO / "tests/reference/database.cpp"), "-o", str(executable),
        ]
        subprocess.run(command, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    if args.fragments:
        args.fragments.write_text(json.dumps([{k: v for k, v in entry.items() if k != "text"}
                                             for entry in selected], indent=2) + "\n")
    print(f"Original C++ database reference written to {args.output}")


if __name__ == "__main__":
    main()
