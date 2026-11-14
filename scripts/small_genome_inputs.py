#!/usr/bin/env python3
"""Prepare pinned public small-genome inputs and deterministic derived queries.

Raw genomes are kept outside Git. Supply local FASTA files or --fetch; canonical
sequence hashes must match upstream/benchmark-inputs.json in either case.
Derived peptides are stop-delimited standard-code ORF fragments, not a claim of
curated annotation. This generator does not plant sequence in the genomes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import gzip

REPO = Path(__file__).resolve().parents[1]
CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = dict(zip((a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG"), CODE))
COMPLEMENT = str.maketrans("ACGTRYMKBDHVSWN", "TGCAYRKMVHDBSWN")


def records(text):
    result, identifier, sequence = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if identifier is not None:
                result.append((identifier, "".join(sequence).upper()))
            identifier, sequence = line[1:].split()[0], []
        elif line.strip():
            if identifier is None:
                raise ValueError("residues before FASTA header")
            sequence.append("".join(line.split()))
    if identifier is not None:
        result.append((identifier, "".join(sequence).upper()))
    if not result or any(not seq for _, seq in result):
        raise ValueError("empty genome or record")
    return result


def fasta(items):
    return "".join(">"+name+"\n"+"\n".join(seq[i:i+80] for i in range(0,len(seq),80))+"\n"
                   for name,seq in items).encode("ascii")


def reverse(seq):
    return seq.translate(COMPLEMENT)[::-1]


def derive(items):
    ranked = sorted(enumerate(items), key=lambda x: (-len(x[1][1]), x[0]))
    nt, peptide, details = [], [], {"nucleotide_queries": [], "protein_queries": []}
    # Eight realistic query fragments, with exact, substituted and gapped cases.
    for index in range(8):
        oid,(name,seq) = ranked[index % len(ranked)]
        length = (250,500,1000,2000)[index % 4]
        start = min(len(seq)-length, (len(seq)*(index+1))//10)
        if start < 0:
            raise ValueError("selected contig too short")
        fragment = seq[start:start+length]
        kind = ("exact", "substitutions", "insertion", "exact")[index % 4]
        if kind == "substitutions":
            chars = list(fragment)
            for at in range(13,len(chars),23):
                chars[at] = {"A":"C","C":"G","G":"T","T":"A"}.get(chars[at],"A")
            fragment = "".join(chars)
        elif kind == "insertion":
            fragment = fragment[:length//2]+"GACTGAC"+fragment[length//2:]
        strand = -1 if index % 2 else 1
        if strand < 0:
            fragment = reverse(fragment)
        qid = f"n{index}_{kind}_{'minus' if strand<0 else 'plus'}"
        nt.append((qid,fragment))
        details["nucleotide_queries"].append(dict(id=qid,subject=name,oid=oid,
            start_zero_based=start,source_length=length,strand=strand,change=kind))
    # Select actual genomic ORFs in all six translation frames. Search a fixed
    # prefix on each strand; longest qualifying runs avoid arbitrary random AA.
    for index in range(8):
        oid,(name,seq) = ranked[index % len(ranked)]
        frame = (1,2,3,-1,-2,-3,1,-1)[index]
        oriented = seq if frame > 0 else reverse(seq)
        phase = abs(frame)-1
        stop = min(len(oriented),150_000)
        aa = "".join(CODONS.get(oriented[i:i+3],"X") for i in range(phase,stop-2,3))
        candidates = []
        begin = 0
        for part in aa.split("*"):
            if len(part)>=120 and "X" not in part:
                candidates.append((len(part), begin, part))
            begin += len(part)+1
        if not candidates:
            raise ValueError(f"no qualifying ORF for {name} frame{frame}")
        _,begin,part = sorted(candidates,key=lambda x:(-x[0],x[1]))[0]
        length = min((150,220,350,450)[index % 4],len(part))
        part = part[:length]
        kind = "exact"
        if index in (2,5):
            chars = list(part)
            for at in range(17,len(chars),31):
                chars[at] = "A" if chars[at] != "A" else "G"
            part, kind = "".join(chars), "substitutions"
        qid = f"p{index}_{kind}_frame{frame}"
        peptide.append((qid,part))
        details["protein_queries"].append(dict(id=qid,subject=name,oid=oid,
            frame=frame,oriented_nucleotide_start=phase+3*begin,source_aa_length=length,
            change=kind,genetic_code=1,annotation="derived stop-delimited ORF fragment"))
    return nt, peptide, details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", action="append", help="Default: all pinned datasets")
    parser.add_argument("--local", action="append", default=[], metavar="NAME=FASTA")
    parser.add_argument("--fetch", action="store_true")
    parser.add_argument("--output", type=Path, default=REPO/"vendor/benchmark-inputs")
    args = parser.parse_args()
    pin = json.loads((REPO/"upstream/benchmark-inputs.json").read_text())
    local = dict(value.split("=",1) for value in args.local)
    for name in args.dataset or pin["datasets"]:
        entry = pin["datasets"][name]
        output = args.output/name
        output.mkdir(parents=True,exist_ok=True)
        path = output/"subjects.fna"
        if name in local:
            text = Path(local[name]).read_text()
        elif path.exists():
            text = path.read_text()
        elif args.fetch:
            with urllib.request.urlopen(entry["url"],timeout=120) as stream:
                blob = stream.read()
            text = (gzip.decompress(blob) if entry.get("gzip") else blob).decode("ascii")
        else:
            raise SystemExit(f"{name}: provide --local {name}=FASTA or --fetch")
        items = records(text)
        canonical = fasta(items)
        sha = hashlib.sha256(canonical).hexdigest()
        if sha != entry["canonical_sha256"]:
            raise SystemExit(f"{name}: sequence hash differs from pin: {sha}")
        path.write_bytes(canonical)
        nt,aa,details = derive(items)
        (output/"nucleotide-queries.fna").write_bytes(fasta(nt))
        (output/"protein-queries.faa").write_bytes(fasta(aa))
        details.update(dataset=name,provenance=entry,subject_count=len(items),
            total_bases=sum(len(seq) for _,seq in items),generator_version=pin["generator_version"])
        details["files"] = {p.name:{"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
            "bytes":p.stat().st_size} for p in (path,output/"nucleotide-queries.fna",output/"protein-queries.faa")}
        (output/"workload.json").write_text(json.dumps(details,indent=2)+"\n")
        print(f"{name}: {len(items)} sequences, {details['total_bases']} bases, {len(nt)} nt/{len(aa)} protein queries")


if __name__ == "__main__":
    main()
