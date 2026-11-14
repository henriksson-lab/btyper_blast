#!/usr/bin/env python3
"""Build unchanged C/C++ query setup, sequence stitching and mask complement."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

from source_ranges import definition
from generate_symdust_reference import fragments as dust_fragments, RNG_CPP

REPO = Path(__file__).resolve().parents[1]
SETUP = "src/algo/blast/api/blast_setup_cxx.cpp"


def fragments(root):
    result = []
    for path, pattern, semicolon in [
        ("src/algo/blast/api/blast_setup.hpp", r"^enum ESentinelType\s*", True),
        (SETUP, r"^static void\ns_QueryInfo_SetContext\(", False),
        (SETUP, r"^TSeqPos CalculateSeqBufferLength\(", False),
        (SETUP, r"^Uint1 GetSentinelByte\(", False),
        (SETUP, r"^static SBlastSequence\nGetSequenceSingleNucleotideStrand\(", False),
        (SETUP, r"^static SBlastSequence\nGetSequenceNucleotideBothStrands\(", False),
        ("src/objects/seqloc/seqloc.asn", r"^Na-strand ::= ENUMERATED\s*", False),
    ]:
        text = (root / path).read_text()
        start, end = definition(text, pattern, semicolon=semicolon)
        data = text[start:end]
        result.append(dict(source=path, byte_start=len(text[:start].encode()),
                           byte_end=len(text[:end].encode()), sha256=hashlib.sha256(data.encode()).hexdigest(), text=data))
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
    dust = dust_fragments(root)
    core = root / "src/algo/blast/core"
    with tempfile.TemporaryDirectory(prefix="btyper-context-oracle-") as temporary:
        directory = Path(temporary)
        asn = next(item["text"] for item in selected if item["source"].endswith("seqloc.asn"))
        variants = re.findall(r"([a-z]+(?:-[a-z]+)?)\s*\((\d+)\)", asn)
        assert len(variants) == 6
        (directory / "reference_na_strand.hpp").write_text(
            "namespace objects { enum ENa_strand {" +
            ",".join("eNa_strand_" + name.replace("-", "_") + "=" + value for name, value in variants) + "}; }\n")
        offset = next(line for line in (root / RNG_CPP).read_text().splitlines() if line.startswith("static const size_t kStateOffset ="))
        (directory / "reference_symdust.hpp").write_text(
            (REPO / "tests/reference/database_shim.hpp").read_text() + "\n" +
            (REPO / "tests/reference/symdust_shim.hpp").read_text() + "\n" + offset + "\n" +
            "\n\n".join(entry["text"] for entry in dust))
        (directory / "reference_context.hpp").write_text(
            (REPO / "tests/reference/context_shim.hpp").read_text() + "\n" +
            "\n\n".join(entry["text"] for entry in selected if not entry["source"].endswith(".asn")))
        common = ["-O2", "-fwrapv", "-ffunction-sections", "-fdata-sections",
                  "-I" + str(REPO / "tests/reference/include"), "-I" + str(root / "include")]
        objects = []
        for name in ("blast_query_info", "blast_program", "blast_encoding", "blast_util", "blast_filter", "pattern", "blast_nalookup", "lookup_util"):
            target = directory / (name + ".o")
            subprocess.run(shlex.split(os.environ.get("CC", "cc")) + ["-std=c11", "-D_DEFAULT_SOURCE"] + common +
                           ["-c", str(core / (name + ".c")), "-o", str(target)], check=True)
            objects.append(str(target))
        executable = directory / "context"
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) + ["-std=c++17"] + common +
                       ["-I" + str(directory), str(REPO / "tests/reference/context.cpp")] + objects +
                       ["-Wl,-dead_strip" if sys.platform == "darwin" else "-Wl,--gc-sections", "-lm", "-o", str(executable)], check=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("wb") as output:
            subprocess.run([str(executable)], stdout=output, check=True)
    if args.fragments:
        args.fragments.write_text(json.dumps([{k: v for k, v in item.items() if k != "text"} for item in selected + dust], indent=2) + "\n")
    print(f"Original query context reference written to {args.output}")


if __name__ == "__main__":
    main()
