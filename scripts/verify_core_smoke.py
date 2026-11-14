#!/usr/bin/env python3
"""Run the existing smoke against raw original core results (Linux x86-64).

Unlike tabular comparison, this checks exact IEEE754 E-value/bit-score bits,
including values that original Seq-align formatting clamps to zero.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from compare_core_results import logical_results, compare

REPO = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    binary = args.native.resolve()
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix="btyper-raw-smoke-"))
    work.mkdir(parents=True,exist_ok=True)
    report = {"native_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "precision": "exact raw core HSP fields and IEEE754 bits before formatting/clamping",
              "status": "running", "programs": {}}
    report_path = args.report or work/"report.json"
    try:
        for mode in ("megablast","blastn","tblastn"):
            query = "protein-queries.faa" if mode == "tblastn" else "nucleotide-queries.fna"
            subject = REPO/"tests/data/smoke/subjects.fna"
            query = REPO/"tests/data/smoke"/query
            original, native = work/(mode+"-original.jsonl"), work/(mode+"-native.jsonl")
            subprocess.run([sys.executable,str(REPO/"scripts/core_result_oracle.py"),
                "--program",mode,"--database",str(subject),"--queries",str(query),
                "--output",str(original),"--work-dir",str(work/(mode+"-oracle"))],check=True)
            output = subprocess.check_output([str(binary),"--raw-json",mode,str(subject),str(query)])
            native.write_bytes(output)
            expected,actual = logical_results(original),logical_results(native)
            differences = compare(expected,actual)
            report["programs"][mode] = {"hsp_count":len(actual["hsps"]),
                "original_sha256":hashlib.sha256(original.read_bytes()).hexdigest(),
                "native_sha256":hashlib.sha256(output).hexdigest(),"differences":differences}
            if differences:
                raise ValueError(f"{mode}: {len(differences)} raw HSP differences")
        report["status"] = "match"
    except Exception as error:
        report["status"],report["error"] = "failed",str(error)
        raise
    finally:
        report_path.parent.mkdir(parents=True,exist_ok=True)
        report_path.write_text(json.dumps(report,indent=2)+"\n")
    print("All25 smoke HSPs match the original core, including exact float bits")


if __name__ == "__main__":
    main()
