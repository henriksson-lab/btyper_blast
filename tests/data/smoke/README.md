# Synthetic BLAST search smoke cases

These six 4096-base subjects contain an 80-residue coding sequence placed in all six subject translation frames. A 180-base nucleotide query is planted on both strands, with an insertion in one subject and ambiguity symbols in another. Independent queries without planted hits and a low complexity query exercise reporting and masking behavior.

Regenerate the FASTA inputs with `python3 scripts/generate_smoke_sequences.py`. This requires the pinned BLAST sources for the standard genetic code table. Regenerate the original tool results with:

```sh
python3 scripts/blast_tools.py fetch
python3 scripts/reference_search.py --database tests/data/smoke/subjects.fna --nucleotide-queries tests/data/smoke/nucleotide-queries.fna --protein-queries tests/data/smoke/protein-queries.faa
```

The recorded output in `tests/fixtures/search-smoke-2.17.0.json` contains six default megablast hits, eight traditional BLASTN hits, and eleven default tblastn hits, including all six subject frames. It preserves original result order, alignment strings, BTOP, and printed numeric scores. BLAST's tabular bit scores and E-values have limited precision. The native database-wide search driver now matches all 25 recorded hits, including order, alignment strings, BTOP and numerical values at original CLI precision. Run `python3 scripts/compare_search.py --native target/release/examples/native_search` after building the example. This is a focused end-to-end integration check; broader parity is still required.

These are synthetic correctness inputs. They do not constitute the representative small-genome performance benchmark required to complete the project.
