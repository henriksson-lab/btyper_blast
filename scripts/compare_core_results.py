#!/usr/bin/env python3
"""Compare logical raw HSP results, including exact floating-point bit patterns.

Original input comes from core_result_oracle.py before Seq-align formatting;
native input comes from native_search --raw-json. Container allocation/storage
metadata is diagnostic and is not an alignment-result contract. No numeric
rounding, E-value clamping, alternate alignment normalization or sorting occurs.
For CLI nucleotide batches, the oracle retains raw local indexes and supplies
an explicit global_list_query_index shadow; only that identifier is rebased.
"""
import argparse
import json
from pathlib import Path

HSP_FIELDS = ("query_index", "subject_list_index", "oid", "hsp_index", "score",
              "num_ident", "num_positives", "context", "num", "comp_adjustment_method",
              "bit_score_bits", "evalue_bits", "query", "subject", "gap_info")
SUBJECT_FIELDS = ("query_index", "subject_list_index", "oid", "list_query_index",
                  "hspcnt", "best_evalue_bits")


def logical_results(path):
    records = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    headers = [r for r in records if r.get("kind") == "results"]
    if len(headers) != 1:
        raise ValueError(f"{path}: expected one results record")
    count = headers[0]["num_queries"]
    queries = {r["query_index"]: r["hsplist_count"] for r in records if r.get("kind") == "query"}
    if any(index not in range(count) for index in queries):
        raise ValueError(f"{path}: invalid query index")
    subjects = []
    for row in records:
        if row.get("kind") != "subject":
            continue
        subject = {key: row[key] for key in SUBJECT_FIELDS}
        if "global_list_query_index" in row:
            subject["list_query_index"] = row["global_list_query_index"]
        subjects.append(subject)
    hsps = [{key: row[key] for key in HSP_FIELDS}
            for row in records if row.get("kind") == "hsp"]
    return {"num_queries": count, "queries": [queries.get(i,0) for i in range(count)],
            "subjects": subjects, "hsps": hsps}


def compare(original, native):
    differences = []
    for field in ("num_queries", "queries"):
        if original[field] != native[field]:
            differences.append(dict(field=field,original=original[field],native=native[field]))
    for kind in ("subjects", "hsps"):
        if len(original[kind]) != len(native[kind]):
            differences.append(dict(field=kind+".length",original=len(original[kind]),native=len(native[kind])))
        for index,(expected,actual) in enumerate(zip(original[kind],native[kind])):
            for key in expected:
                if expected[key] != actual[key]:
                    differences.append(dict(field=f"{kind}[{index}].{key}",original=expected[key],native=actual[key]))
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", required=True, type=Path)
    parser.add_argument("--native-json", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    original, native = logical_results(args.original), logical_results(args.native_json)
    differences = compare(original,native)
    report = dict(status="mismatch" if differences else "match", hsp_count=len(native["hsps"]),
                  original=str(args.original.resolve()),native=str(args.native_json.resolve()),
                  precision="exact IEEE754 bit patterns; all logical HSP fields and ordered lists",
                  differences=differences)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+"\n")
    if differences:
        for error in differences[:20]:
            print(json.dumps(error))
        raise SystemExit(f"{len(differences)} raw core result differences")
    print(f"{len(native['hsps'])} raw HSPs match, including full-precision float bits, frames, seeds and edit scripts")


if __name__ == "__main__":
    main()
