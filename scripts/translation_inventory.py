#!/usr/bin/env python3
"""Generate a scoped function map from a comparator C report (not a scope closure)."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[1]
CORE = "src/algo/blast/core/"

# A reviewed starting shortlist. Do not expand this to every parsed function.
CANDIDATES = {
    "blast_program.c": (
        "Program classification for nucleotide and translated-subject paths",
        "Blast_QueryIsProtein Blast_QueryIsNucleotide Blast_SubjectIsProtein "
        "Blast_SubjectIsNucleotide Blast_QueryIsTranslated Blast_SubjectIsTranslated "
        "Blast_ProgramIsNucleotide Blast_ProgramIsValid",
    ),
    "gapinfo.c": (
        "Traceback edit-script storage and operation accumulation",
        "GapEditScriptNew GapEditScriptDelete s_GapPrelimEditBlockRealloc "
        "s_GapPrelimEditBlockAddNew GapPrelimEditBlockAdd GapPrelimEditBlockNew "
        "GapPrelimEditBlockFree GapPrelimEditBlockReset",
    ),
    "blast_util.c": (
        "Sequence encoding, ambiguity-aware translation, and frame bookkeeping",
        "s_CodonToAA BLAST_GetTranslation GetReverseNuclSequence BLAST_PackDNA "
        "BLAST_ContextToFrame BLAST_FrameToContext BLAST_GetNumberOfContexts "
        "BLAST_GetTranslatedProteinLength s_BlastGetTranslationTable "
        "BLAST_TranslateCompressedSequence BLAST_GetAllTranslations Blast_GetPartialTranslation",
    ),
    "ncbi_math.c": (
        "Small numerical kernels used by statistical calculations; verify each caller before porting",
        "BLAST_Expm1 BLAST_Log1p BLAST_Gcd BLAST_Gdb3 BLAST_Nint BLAST_Powi",
    ),
    "boost_erf.c": (
        "Error functions used by the statistical path; retain Boost license and recursive implementation",
        "ErfImpl Erf ErfC",
    ),
}


def snake_case(name):
    name = re.sub(r"([A-Z]+)([A-Z][a-z])", r"\1_\2", name)
    return re.sub(r"([a-z0-9])([A-Z])", r"\1_\2", name).lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("report", type=Path)
    parser.add_argument("--order", type=Path, help="Comparator order CSV; advisory only")
    args = parser.parse_args()
    pin = json.loads((REPO / "upstream/blast-source.json").read_text())
    root = REPO / pin["root"]
    report = json.loads(args.report.read_text())
    ledger_path = REPO / "upstream/translations.json"
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else []
    translated = {(entry["source"], entry["symbol"]): entry for entry in ledger}
    order = {}
    if args.order:
        with args.order.open() as stream:
            for rank, entry in enumerate(csv.DictReader(stream), 1):
                order[(Path(entry["file"]).name, entry["name"])] = rank
    rows, seen = [], set()
    for filename, (reason, names) in CANDIDATES.items():
        for symbol in names.split():
            matches = [f for f in report["functions"]
                       if f["name"] == symbol and
                       Path(f["location"]["file"]).as_posix().endswith(CORE + filename)]
            if len(matches) != 1:
                raise SystemExit(f"Expected one definition of {filename}:{symbol}, got {len(matches)}")
            function = matches[0]
            loc = function["location"]
            raw = (root / CORE / filename).read_bytes()
            rust_name = snake_case(symbol)
            rust_path = "src/core/" + Path(filename).stem + ".rs"
            key = (rust_path, rust_name)
            if key in seen:
                raise SystemExit(f"Naming collision: {key}")
            seen.add(key)
            entry = translated.get((CORE + filename, symbol))
            rows.append({
                "version": pin["version"], "source": CORE + filename,
                "symbol": symbol, "line_start": loc["line_start"], "line_end": loc["line_end"],
                "byte_start": loc["byte_start"], "byte_end": loc["byte_end"],
                "source_sha256": hashlib.sha256(raw[loc["byte_start"]:loc["byte_end"]]).hexdigest(),
                "rust_path": rust_path, "rust_type": "", "rust_symbol": rust_name,
                "status": entry["status"] if entry else "candidate",
                "reason": reason, "adaptations": entry["adaptations"] if entry else "",
                "advisory_order": order.get((filename, symbol), ""),
                "direct_calls": ";".join(c["callee"] for c in function["calls"]),
            })
    included = {(row["source"], row["symbol"]) for row in rows}
    for entry in ledger:
        if (entry["source"], entry["symbol"]) not in included:
            rows.append({key: entry.get(key, "") for key in rows[0]})
    with (REPO / "upstream/translation-candidates.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Recorded {len(rows)} scoped entries; {sum(r['status'] == 'translated' for r in rows)} translated")


if __name__ == "__main__":
    main()
