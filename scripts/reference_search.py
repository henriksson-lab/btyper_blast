#!/usr/bin/env python3
"""Run the pinned original BLAST defaults and capture search results in memory.

This is a correctness oracle, not a speed comparison. Tabular bit scores and
E-values retain BLAST's printed precision; full precision needs the core oracle.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[1]
FIELDS = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qframe", "sframe",
          "length", "nident", "mismatch", "gaps", "gapopen", "score", "bitscore", "evalue",
          "btop", "qseq", "sseq"]
INTEGER_FIELDS = set(FIELDS[2:14])


def run_reference(database, nucleotide_queries, protein_queries):
    subprocess.run([os.sys.executable, str(REPO / "scripts/blast_tools.py"), "verify"],
                   check=True, stdout=subprocess.DEVNULL)
    pin = json.loads((REPO / "upstream/blast-tools.json").read_text())
    binary = REPO / pin["root"]
    env = dict(os.environ, BLAST_USAGE_REPORT="false", NCBI_DONT_USE_NCBIRC="1",
               NCBI_DONT_USE_LOCAL_CONFIG="1")
    result = {
        "reference_version": pin["version"], "executable_archive_sha256": pin["sha256"],
        "inputs": {name: {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                   for name, path in (("database", database), ("nucleotide_queries", nucleotide_queries),
                                      ("protein_queries", protein_queries))},
        "precision": "Original CLI tabular precision for bitscore/evalue; integers, BTOP and alignment strings preserved",
        "searches": {},
    }
    with tempfile.TemporaryDirectory(prefix="btyper-blast-search-") as temporary:
        db = Path(temporary) / "nucleotide"
        build = subprocess.run([str(binary / "makeblastdb"), "-in", str(database), "-dbtype", "nucl",
                                "-out", str(db)], capture_output=True, text=True, check=True,
                               cwd=temporary, env=env)
        result["database_format"] = 5
        # Default executable task, explicitly selected traditional task, and tblastn.
        for name, program, query, task_args in (
            ("blastn_default_megablast", "blastn", nucleotide_queries, []),
            ("blastn_traditional", "blastn", nucleotide_queries, ["-task", "blastn"]),
            ("tblastn_default", "tblastn", protein_queries, []),
        ):
            command = [str(binary / program), "-query", str(query), "-db", str(db),
                       "-num_threads", "1", "-outfmt", "6 " + " ".join(FIELDS)] + task_args
            search = subprocess.run(command, capture_output=True, text=True, check=True,
                                    cwd=temporary, env=env)
            hits = []
            for line in search.stdout.splitlines():
                values = line.split("\t")
                if len(values) != len(FIELDS):
                    raise ValueError(f"Unexpected BLAST result field count: {line}")
                hits.append({field: int(value) if field in INTEGER_FIELDS else value
                             for field, value in zip(FIELDS, values)})
            result["searches"][name] = {"hits": hits, "diagnostics": search.stderr}
        # Check the builder really emitted the requested modern/default format.
        info = subprocess.run([str(binary / "blastdbcmd"), "-db", str(db), "-info"],
                              capture_output=True, text=True, check=True, cwd=temporary, env=env)
        if "BLASTDB Version: 5" not in info.stdout:
            raise ValueError("Reference makeblastdb did not create a v5 database")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--nucleotide-queries", required=True, type=Path)
    parser.add_argument("--protein-queries", required=True, type=Path)
    args = parser.parse_args()
    result = run_reference(args.database.resolve(), args.nucleotide_queries.resolve(),
                           args.protein_queries.resolve())
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
