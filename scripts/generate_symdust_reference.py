#!/usr/bin/env python3
"""Compile the original symmetric DUST class without object-manager serialization."""
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
RNG_CPP = "src/util/random_gen.cpp"
RNG_HPP = "include/util/random_gen.hpp"
DUST_CPP = "src/algo/dustmask/symdust.cpp"
DUST_HPP = "include/algo/dustmask/symdust.hpp"


def fragments(root):
    selections = [
        (RNG_CPP, r"^CRandom::CRandom\(EGetRandMethod method\)", False),
        (RNG_CPP, r"^CRandom::CRandom\(TValue seed\)", False),
        (RNG_CPP, r"^void CRandom::SetSeed\(TValue seed\)", False),
        (RNG_CPP, r"^void CRandom::Reset\(void\)", False),
        (RNG_HPP, r"^inline CRandom::TValue CRandom::x_GetRand32Bits\(void\)", False),
        (RNG_HPP, r"^inline CRandom::TValue CRandom::GetRand\(void\)", False),
        (DUST_HPP, r"^class NCBI_XALGODUSTMASK_EXPORT CSymDustMasker\s*", True),
        (DUST_CPP, r"^CSymDustMasker::triplets::triplets\(", False),
        (DUST_CPP, r"^bool CSymDustMasker::triplets::shift_high\(", False),
        (DUST_CPP, r"^bool CSymDustMasker::triplets::shift_window\(", False),
        (DUST_CPP, r"^inline void CSymDustMasker::triplets::find_perfect\(", False),
        (DUST_CPP, r"^CSymDustMasker::CSymDustMasker\(", False),
        (DUST_CPP, r"^inline void CSymDustMasker::save_masked_regions\(", False),
        (DUST_CPP, r"^std::unique_ptr< CSymDustMasker::TMaskList >\s*\nCSymDustMasker::operator\(\)\( const sequence_type & seq,", False),
        (DUST_CPP, r"^std::unique_ptr< CSymDustMasker::TMaskList >\s*\nCSymDustMasker::operator\(\)\( const sequence_type & seq \)", False),
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
    parser.add_argument("--fragments", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, str(REPO / "scripts/blast_source.py"), "verify"], check=True)
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    root = REPO / pin["root"]
    selected = fragments(root)
    offset = next(line for line in (root / RNG_CPP).read_text().splitlines()
                  if line.startswith("static const size_t kStateOffset ="))
    with tempfile.TemporaryDirectory(prefix="btyper-dust-oracle-") as temporary:
        directory = Path(temporary)
        shim = (REPO / "tests/reference/database_shim.hpp").read_text() + "\n" + (REPO / "tests/reference/symdust_shim.hpp").read_text()
        (directory / "reference_symdust.hpp").write_text(shim + "\n" + offset + "\n\n" + "\n\n".join(entry["text"] for entry in selected))
        executable = directory / "symdust"
        command = shlex.split(os.environ.get("CXX", "c++")) + ["-std=c++17", "-O2", "-I" + str(directory),
                    str(REPO / "tests/reference/symdust.cpp"), "-o", str(executable)]
        subprocess.run(command, check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    if args.fragments:
        args.fragments.write_text(json.dumps([{k: v for k, v in entry.items() if k != "text"}
                                             for entry in selected], indent=2) + "\n")
    print(f"Original symmetric DUST reference written to {args.output}")


if __name__ == "__main__":
    main()
