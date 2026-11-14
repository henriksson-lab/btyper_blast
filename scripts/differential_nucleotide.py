#!/usr/bin/env python3
"""Focused deterministic whole-pipeline nucleotide comparisons with pinned BLAST.

Inputs and both outputs remain in a temporary work directory for reproduction.
This is a small integration workload, not exhaustive parity or a speed benchmark.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile

from compare_search import compare
from reference_search import FIELDS

REPO = Path(__file__).resolve().parents[1]
ENV = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
           NCBI_DONT_USE_LOCAL_CONFIG="1")


def rc(sequence):
    return sequence.translate(str.maketrans("ACGTRYMKBDHVN", "TGCAYRKMVHDBN"))[::-1]


def workloads():
    rng = random.Random(217001)

    def dna(length):
        return "".join(rng.choice("ACGT") for _ in range(length))

    def mutated(sequence, step):
        return "".join(rng.choice([base for base in "ACGT" if base != value])
                       if i % step == 0 else value for i, value in enumerate(sequence))

    a, b, c = dna(240), dna(170), dna(300)
    return {
        "exact_mutated": ([("a", a), ("b", b)],
                            [("exact", a), ("changed", mutated(a, 29)), ("other", b), ("none", dna(350))]),
        "gaps": ([("query", c)], [("ins1", c[:80] + "A" + c[80:]),
                                      ("ins7", c[:140] + "ATGCTTA" + c[140:]),
                                      ("del3", c[:170] + c[173:])]),
        "reverse": ([("forward", a), ("reverse", rc(b))],
                      [("reverse_a", dna(31) + rc(a) + dna(40)), ("forward_b", b)]),
        "short_queries": ([(f"q{n}", a[:n]) for n in (10, 14, 20, 28, 40, 55)],
                            [("target", a), ("short", b[:9])]),
        "all_short_queries": ([("ten", a[:10]), ("seven", b[:7])], [("target", a), ("other", b)]),
        "all_masked_queries": ([("poly_a", "A" * 120), ("poly_ac", "AC" * 100)], [("target", a), ("other", b)]),
        "ambiguous_subject": ([("query", c)],
                                [("n_gap", c[:90] + "N" * 8 + c[98:]),
                                 ("iupac", c[:140] + "RYSWKMBDHVN" + c[151:]), ("exact", c)]),
        "ambiguous_query": ([("ambig", c[:100] + "RYSWKMBDHVN" + c[111:]),
                              ("exact", a), ("all_n", "N" * 80)], [("a", a), ("c", c)]),
        "dust_flanks": ([("low_flanks", "A" * 90 + b + "AT" * 60), ("simple", "AC" * 100)],
                          [("target", "A" * 90 + b + "AT" * 60), ("embedded", dna(40) + b + dna(40))]),
        "repeated_blocks": ([("a", a[:100]), ("b", b[:80])],
                              [("copies", a[:100] + dna(45) + a[:100] + dna(55) + rc(a[:100])),
                               ("pair", b[:80] + dna(40) + a[:100])]),
        "ordering": ([("one", a), ("two", b)],
                       [(f"tie_{i}", a) for i in range(5)] + [("two", b), ("changed", mutated(a, 37))]),
        "ends_and_case": ([("mixed", a[:30].lower() + a[30:]), ("partial", b[20:140])],
                           [("at_start", a + dna(75)), ("at_end", dna(83) + a),
                            ("reverse_partial", rc(b[20:140]))]),
    }


def write_fasta(path, records):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def run(command, directory):
    return subprocess.run([str(value) for value in command], cwd=directory, env=ENV,
                          capture_output=True, text=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--case", action="append", choices=workloads())
    parser.add_argument("--program", action="append", choices=("blastn", "megablast"))
    args = parser.parse_args()
    subprocess.run([os.sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"],
                   check=True, stdout=subprocess.DEVNULL)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    binary = REPO / pin["root"]
    version = run([binary / "blastn", "-version"], REPO)
    if version.returncode or "blastn: " + pin["version"] not in version.stdout:
        raise ValueError("reference executable version differs from pin")
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix="btyper-nucleotide-differential-"))
    work.mkdir(parents=True, exist_ok=True)
    report = {"version": pin["version"], "archive_sha256": pin["sha256"], "work_dir": str(work), "native": str(args.native.resolve()),
              "native_sha256": hashlib.sha256(args.native.read_bytes()).hexdigest(), "runs": []}
    for name, (queries, subjects) in workloads().items():
        if args.case and name not in args.case:
            continue
        directory = work / name
        directory.mkdir(exist_ok=True)
        query, subject = directory / "query.fna", directory / "subject.fna"
        write_fasta(query, queries)
        write_fasta(subject, subjects)
        built = run([binary / "makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", directory / "db"], directory)
        if built.returncode:
            raise RuntimeError(built.stderr)
        for program in args.program or ("blastn", "megablast"):
            original_command = [binary / "blastn", "-query", query, "-db", directory / "db",
                                "-num_threads", "1", "-outfmt", "6 " + " ".join(FIELDS)]
            if program == "blastn":
                original_command += ["-task", "blastn"]
            native_command = [args.native.resolve(), program, subject, query]
            original, native = run(original_command, directory), run(native_command, directory)
            for label, result in (("original", original), ("native", native)):
                (directory / f"{program}-{label}.tsv").write_text(result.stdout)
                (directory / f"{program}-{label}.stderr").write_text(result.stderr)
            entry = {"case": name, "program": program,
                     "reference_command": list(map(str, original_command)), "native_command": list(map(str, native_command))}
            try:
                if original.returncode or native.returncode:
                    raise ValueError(f"exit status original={original.returncode}, native={native.returncode}; see stderr files")
                expected = [dict(zip(FIELDS, line.split("\t"))) for line in original.stdout.splitlines()]
                actual = [line.split("\t") for line in native.stdout.splitlines()]
                compare(actual, expected, program)
                entry.update(status="match", hits=len(actual))
                print(f"{name}/{program}: {len(actual)} hits match")
            except ValueError as error:
                entry.update(status="different", difference=str(error))
                print(f"{name}/{program}: DIFFERENT: {error}")
            report["runs"].append(entry)
    (work / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    failures = sum(item["status"] != "match" for item in report["runs"])
    print(f"{len(report['runs'])} runs, {failures} differences; reproduction inputs: {work}")
    raise SystemExit(bool(failures))


if __name__ == "__main__":
    main()
