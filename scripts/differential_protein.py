#!/usr/bin/env python3
"""Small deterministic whole-default-tblastn differential workloads.

Execute the native in-memory API driver and pinned BLAST+ 2.17.0+ against the
same FASTA/database inputs. Compare every field and hit order using the exact
original tabular float formatting rules; never weaken comparisons on failure.
Failed inputs/outputs are retained for source investigation, not added as a
large tracked fixture corpus. This is neither a complete parity nor speed gate.
"""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

from compare_search import compare
from reference_search import FIELDS

REPO = Path(__file__).resolve().parents[1]
AMINO = "ARNDCQEGHILKMFPSTWYV"
CODONS = dict(zip(AMINO, ("GCT", "CGT", "AAT", "GAT", "TGT", "CAA", "GAA", "GGT",
                          "CAT", "ATT", "CTG", "AAA", "ATG", "TTT", "CCT", "TCT",
                          "ACT", "TGG", "TAT", "GTT")))


def protein(length, seed):
    rng = random.Random(seed)
    return "".join(rng.choice(AMINO) for _ in range(length))


def encode(sequence):
    return "".join(CODONS[aa] for aa in sequence)


def reverse(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mutate(sequence, step):
    result = list(sequence)
    for pos in range(7, len(result), step):
        result[pos] = AMINO[(AMINO.index(result[pos]) + 7) % len(AMINO)]
    return "".join(result)


def workloads():
    p = protein(105, 217001)
    q = protein(92, 217002)
    six = [(f"frame{sign}{frame+1}", dna) for sign in ("p", "m") for frame in range(3)
           for dna in [("A" * frame + encode(p) + "TAA")]]
    six = [(name, reverse(dna) if name.startswith("framem") else dna) for name, dna in six]
    yield "six_frames", [("q", p)], six
    yield "ranking", [("q", p)], [("exact", encode(p)), ("other", encode(q))] + [
        (f"mut{step}", encode(mutate(p, step))) for step in (7, 9, 12, 15, 21)] + [
        ("equal_last", encode(p))]
    yield "insertions", [("q", p)], [(f"ins{n}", encode(p[:47] + q[:n] + p[47:]))
                                       for n in (1, 3, 7, 15)]
    yield "deletions", [("q", p)], [(f"del{n}", encode(p[:43] + p[43+n:]))
                                      for n in (1, 3, 7, 15)]
    yield "multiple_hsps", [("p", p), ("q", q)], [
        ("repeat", encode(p) + "TAATAGTGA" + encode(q) + "TAATAGTGA" + encode(p)),
        ("reverse", reverse(encode(q) + "TAATAGTGA" + encode(p)))]
    ambiguous = encode(p)
    yield "ambiguous_codons", [("q", p)], [
        ("internal_N", ambiguous[:117] + "NNN" + ambiguous[120:]),
        ("ambiguous_run", ambiguous[:117] + "N" * 18 + ambiguous[135:]),
        ("end_N", "NNN" + ambiguous[3:-3] + "NNN")]
    masked = p[:40] + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + p[40:]
    yield "masked_query", [("masked", masked), ("original", p)], [
        ("masked_exact", encode(masked)), ("diverse", encode(p))]
    biased = "".join("EK" if i % 3 else "DE" for i in range(45))
    mixed = p[:35] + biased + p[35:]
    yield "composition_bias", [("mixed", mixed), ("diverse", p)], [
        ("mixed", encode(mixed)), ("mutated", encode(mutate(mixed, 13))),
        ("diverse", encode(p))]
    yield "multiple_queries", [("p", p), ("short", q[:39]), ("q", q)], [
        ("p", encode(p)), ("q", encode(q)), ("all", encode(q)+"TAATAGTGA"+encode(p))]
    rare = p[:20] + "XBZJOU" + p[26:]
    yield "protein_ambiguities", [("rare", rare), ("plain", p)], [
        ("plain", encode(p)), ("substituted", encode(p[:20]+"ADELIC"+p[26:]))]
    yield "short_queries", [("eighteen", p[:18]), ("thirty", p[:30]), ("p", p)], [
        ("exact", encode(p)), ("terminal", "TAA"+encode(p[-30:])+"TGA")]
    yield "negative_ambiguous", [("q", p)], [
        ("minus2", reverse("A"+ambiguous[:180]+"NNN"+ambiguous[183:]+"TAA")),
        ("minus3", reverse("AA"+encode(mutate(p, 17))+"TAG"))]


def fasta(path, records):
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--case", action="append", help="Repeat to select named workloads")
    parser.add_argument("--workdir", type=Path, help="Retain all generated inputs/output here")
    args = parser.parse_args()
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    source = json.loads((REPO / "upstream/blast-source.json").read_text())
    if pin["version"] != source["version"]:
        raise ValueError("reference tools and translated source pins differ")
    subprocess.run([os.sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"],
                   check=True, stdout=subprocess.DEVNULL)
    reference = REPO / pin["root"]
    native = args.native.resolve()
    work = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="btyper-protein-diff-"))
    work.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
               NCBI_DONT_USE_LOCAL_CONFIG="1")
    failures = []
    selected = [(name, queries, subjects) for name, queries, subjects in workloads()
                if not args.case or name in args.case]
    if args.case and set(args.case) - {name for name, _, _ in selected}:
        raise ValueError("unknown case requested")
    for name, queries, subjects in selected:
        case = work / name
        case.mkdir(exist_ok=True)
        query, subject = case / "queries.faa", case / "subjects.fna"
        fasta(query, queries)
        fasta(subject, subjects)
        try:
            build = subprocess.run([str(reference / "makeblastdb"), "-in", str(subject),
                                    "-dbtype", "nucl", "-out", str(case / "db")],
                                   capture_output=True, text=True, check=True, env=env)
            (case / "makeblastdb.txt").write_text(build.stdout + build.stderr)
            original = subprocess.run([str(reference / "tblastn"), "-query", str(query),
                                       "-db", str(case / "db"), "-num_threads", "1",
                                       "-outfmt", "6 " + " ".join(FIELDS)],
                                      capture_output=True, text=True, check=True, env=env)
            (case / "original.tsv").write_text(original.stdout)
            (case / "original.stderr").write_text(original.stderr)
            result = subprocess.run([str(native), "tblastn", str(subject), str(query)],
                                    capture_output=True, text=True, env=env)
            (case / "native.tsv").write_text(result.stdout)
            (case / "native.stderr").write_text(result.stderr)
            result.check_returncode()
            expected = [dict(zip(FIELDS, row.split("\t"))) for row in original.stdout.splitlines()]
            rows = [row.split("\t") for row in result.stdout.splitlines()]
            compare(rows, expected, name)
            print(f"{name}: {len(rows)} hits match BLAST+ {pin['version']}", flush=True)
        except (ValueError, subprocess.CalledProcessError) as error:
            failures.append(name)
            print(f"FAIL {name}: {error}\n  retained input/output: {case}", flush=True)
    if failures:
        raise SystemExit(f"{len(failures)}/{len(selected)} workloads failed: {', '.join(failures)}; inputs retained at {work}")
    print(f"All {len(selected)} small default-tblastn workloads agree; this is not a complete parity/speed claim.")
    if not args.workdir:
        shutil.rmtree(work)


if __name__ == "__main__":
    main()
