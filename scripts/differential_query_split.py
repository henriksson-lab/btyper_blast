#!/usr/bin/env python3
"""Reproduce two default long-query boundary comparisons with pinned BLAST.

Generates deterministic FASTA inputs instead of storing large fixtures. Both
classic blastn (2,000,800nt) and tblastn (40,200aa) cross the original query
split boundary and contain a later hit. This is a focused parity check, not a
benchmark or exhaustive corpus. All inputs, outputs and reports are retained.
Optional --core uses the existing Linux x86-64 GDB child-process oracle; run it
only where the caller permits child debugging/ptrace.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

from compare_search import compare
from reference_search import FIELDS

REPO = Path(__file__).resolve().parents[1]
SEED = 21717
ENV = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
           NCBI_DONT_USE_LOCAL_CONFIG="1")
# Exercise default source geometry even if a caller previously experimented
# with the original toolkit's environment overrides.
for _variable in ("CHUNK_SIZE", "OVERLAP_CHUNK_SIZE"):
    ENV.pop(_variable, None)


def workloads():
    """Exact original boundary-guard generator; preserve RNG consumption order."""
    rng = random.Random(SEED)
    nucleotide = "".join(rng.choices("ACGT", k=2_000_800))
    mid = len(nucleotide) // 2
    subject = ("".join(rng.choices("ACGT", k=300))
               + nucleotide[mid - 160:mid + 180]
               + "".join(rng.choices("ACGT", k=300))
               + nucleotide[-500:-250])
    letters = "ARNDCQEGHILKMFPSTWYV"
    protein = "".join(rng.choices(letters, k=40_200))
    mid = len(protein) // 2
    codons = dict(zip(letters, (
        "GCT", "CGT", "AAT", "GAT", "TGT", "CAA", "GAA", "GGT", "CAT", "ATT",
        "CTG", "AAA", "ATG", "TTT", "CCT", "TCT", "ACT", "TGG", "TAT", "GTT")))
    translated = protein[mid - 180:mid + 180] + "GGSGGS" + protein[-420:-220]
    coding = "".join(codons[residue] for residue in translated)
    return (
        ("blastn", "nucleotide-query.fna", "long_nt", nucleotide,
         "nucleotide-subject.fna", "subject", subject),
        ("tblastn", "protein-query.faa", "long_aa", protein,
         "protein-subject.fna", "coding", coding),
    )


def fasta_bytes(identifier, sequence):
    return (f">{identifier}\n" + "".join(sequence[i:i + 80] + "\n"
                                        for i in range(0, len(sequence), 80))).encode("ascii")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command, directory, label):
    result = subprocess.run(list(map(str, command)), cwd=directory, env=ENV,
                            capture_output=True, text=True)
    (directory / f"{label}.stdout").write_text(result.stdout)
    (directory / f"{label}.stderr").write_text(result.stderr)
    if result.returncode:
        raise ValueError(f"{label} exited {result.returncode}; see {directory / (label + '.stderr')}")
    return result


def core_check(native, program, subject, query, directory):
    original_json = directory / "original.jsonl"
    native_json = directory / "native.jsonl"
    core_report = directory / "core-report.json"
    oracle = [sys.executable, REPO / "scripts/core_result_oracle.py",
              "--program", program, "--database", subject, "--queries", query,
              "--output", original_json, "--work-dir", directory / "oracle"]
    run(oracle, directory, "core-oracle")
    result = run([native, "--raw-json", program, subject, query], directory, "native-core")
    native_json.write_text(result.stdout)
    run([sys.executable, REPO / "scripts/compare_core_results.py",
         "--original", original_json, "--native-json", native_json,
         "--report", core_report], directory, "core-compare")
    return json.loads(core_report.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--work-dir", type=Path,
                        help="retain generated inputs and outputs here; otherwise create a retained temporary directory")
    parser.add_argument("--core", action="store_true",
                        help="also compare all logical core fields/float bits using the Linux GDB oracle")
    args = parser.parse_args()
    native = args.native.resolve()
    if not native.is_file():
        parser.error(f"native executable does not exist: {native}")
    subprocess.run([sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"],
                   check=True, stdout=subprocess.DEVNULL, env=ENV)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    reference = REPO / pin["root"]
    work = (args.work_dir.resolve() if args.work_dir else
            Path(tempfile.mkdtemp(prefix="btyper-query-split-differential-")))
    work.mkdir(parents=True, exist_ok=True)
    report = dict(version=pin["version"], archive_sha256=pin["sha256"], seed=SEED,
                  purpose="two default query-split boundary cases; parity only",
                  work_dir=str(work), native=str(native), native_sha256=sha256(native),
                  geometry_overrides_removed=["CHUNK_SIZE", "OVERLAP_CHUNK_SIZE"], runs=[])
    for program, qname, qid, sequence, sname, sid, target in workloads():
        directory = work / program
        directory.mkdir(exist_ok=True)
        query, subject = directory / qname, directory / sname
        query.write_bytes(fasta_bytes(qid, sequence))
        subject.write_bytes(fasta_bytes(sid, target))
        entry = dict(program=program, query_length=len(sequence), subject_length=len(target),
                     query=str(query), subject=str(subject), query_sha256=sha256(query),
                     subject_sha256=sha256(subject), reference_sha256=sha256(reference / program))
        report["runs"].append(entry)
        try:
            version = run([reference / program, "-version"], directory, "reference-version")
            if f"{program}: {pin['version']}" not in version.stdout:
                raise ValueError("reference executable version differs from pin")
            run([reference / "makeblastdb", "-in", subject, "-dbtype", "nucl",
                 "-out", directory / "db"], directory, "makeblastdb")
            original_command = [reference / program, "-query", query, "-db", directory / "db",
                                "-num_threads", "1", "-outfmt", "6 " + " ".join(FIELDS)]
            if program == "blastn":
                original_command += ["-task", "blastn"]
            native_command = [native, program, subject, query]
            entry.update(reference_command=list(map(str, original_command)),
                         native_command=list(map(str, native_command)))
            original = run(original_command, directory, "original")
            result = run(native_command, directory, "native")
            (directory / "original.tsv").write_text(original.stdout)
            (directory / "native.tsv").write_text(result.stdout)
            expected = [dict(zip(FIELDS, line.split("\t"))) for line in original.stdout.splitlines()]
            actual = [line.split("\t") for line in result.stdout.splitlines()]
            compare(actual, expected, program + " default query split")
            entry.update(status="match", hits=len(actual), fields=FIELDS)
            if args.core:
                entry["core"] = core_check(native, program, subject, query, directory)
            print(f"{program}: {len(actual)} hits match, including order and all {len(FIELDS)} fields", flush=True)
        except (ValueError, OSError) as error:
            entry.update(status="different", difference=str(error))
            print(f"{program}: DIFFERENT: {error}", flush=True)
        (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    failures = sum(entry["status"] != "match" for entry in report["runs"])
    report["status"] = "different" if failures else "match"
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"{len(report['runs'])} boundary cases, {failures} differences; retained inputs/report: {work}")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
