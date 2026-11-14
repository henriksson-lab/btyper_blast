# BLAST Rust implementation progress

The goal remains open: standard BLAST+ 2.17.0+ results from a usable Rust library, memory databases/results, verified Linux/Windows/macOS CI and a repeatable speed improvement. The pinned algorithms and systematic function correspondence remain authoritative. See [goal.md](goal.md).

Use parallel workers for independent work. While a required translation is missing, prioritize implementation and pipeline assembly, with focused checks. Complete-result and performance validation are completion gates.

## Implemented scope

- **596 upstream functions/methods** have source ranges, definition hashes, Rust mappings and adaptations in [translations.json](translations.json). Data/types/constants are tracked separately in [translation-data.json](translation-data.json).
- `DatabaseBuilder` and `NucleotideDatabase::from_fasta` create memory subjects with original packed bytes, ambiguity restoration and database statistics. Packed payloads match official makeblastdb fixtures. Disk/NCBI v5 persistence is outside scope by user instruction.
- `PreparedNucleotideBatch` retains DUST soft masks, both query strands and both nucleotide presets. `PreparedProteinBatch` retains SEG masks and original unmasked copies.
- `search_nucleotide_database` includes original seed/initial-HSP search, gapped preliminary evaluation, candidate collector limits, selected-list traceback, ambiguity/endpoint filtering, statistics and final ranking/limits. `search_protein_database` includes six-frame translation, AA two-hit search, gapped frame merging, the collector, actual CBS2 optimization/redo, linked HSP statistics and composition heaps.
- Chunking uses the original five-million-residue limit and 100-residue overlap. Nucleotide chunks borrow original packed bytes; translated chunks retain physical offsets within complete frames. There is no artificial five-million-nucleotide tblastn cap.
- The source partial target translation cache and reusable protein/nucleotide scratch owners retain original lifetimes. Large query batches use original chunk geometry, inherited DUST/local DUST union or fresh local SEG, full effective-space prefix semantics and preliminary overlap merging before one full-query traceback/CBS. Full-query final setup skips the unused full lookup. Nucleotide traceback restores each OID once; CBS matrix rows persist across alignments. Safe AA scanner changes bind complete input/output spans while preserving loops, order and restart behavior. No unchecked indexing was added.
- Owned HSP results retain core coordinates and full-precision fields. `AlignmentReport` reconstructs oriented endpoints, alignment strings and BTOP. Original-core tblastn identities/positives use masked query residues; formatter identities use restored aligned residues, so these counts can differ.

## Validation evidence

The frozen combined release artifact passes the synthetic **25-HSP** smoke and **1,259 real-genome HSPs** against unmodified official core results. All logical ordered result fields match, including raw scores, identities/positives, contexts/frames, seed coordinates, edit scripts and exact floating-point bits. Source query batching is mapped explicitly; it does not discard early batches or normalize scores, alignments or statistics.

| Pinned workload | Megablast | Blastn | Tblastn |
| --- | ---: | ---: | ---: |
| *S. epidermidis* ST14 | 6 | 481 | 52 |
| *S. cerevisiae* S288C | 16 | 573 | 131 |

Current local proof records are `/tmp/btyper-final-real-core-parity.json` and `/tmp/btyper-final-core-smoke/report.json`; they record the tested binary and reference output hashes. These temporary records are evidence for that artifact, not a substitute for rerunning validation after changes.

[small_genome_inputs.py](../scripts/small_genome_inputs.py) fetches pinned public genomes or accepts hash-checked local inputs and derives queries without planting sequence. [differential_nucleotide.py](../scripts/differential_nucleotide.py) covers 12 deterministic workloads with both presets; [differential_protein.py](../scripts/differential_protein.py) covers 12 workloads/61 hits. Cases include short/masked queries, ambiguity, gaps, reverse frames and ordering.

[core_result_oracle.py](../scripts/core_result_oracle.py) reads official core results through child GDB, using compiled pinned-header layouts. It captures each completed final traceback batch before Seq-align E-value clamping. [compare_core_results.py](../scripts/compare_core_results.py) checks exact logical fields and float bits. [verify_core_smoke.py](../scripts/verify_core_smoke.py) runs the smoke automatically, including in configured Linux CI. CLI/display comparisons remain separate.

Focused checks also cover the preliminary 550/final 500 gate on 551 matching subjects, a 5,000,400-nt nucleotide boundary match, and a translated 5M-AA boundary match in a 15,000,300-nt subject. Earlier source-C/C++ component fixtures remain intact; integration adaptations and focused checks are recorded under [integrations/](integrations/).

## Reproduce

```sh
cargo check --all-targets --locked
cargo test --locked
cargo test --release --locked
cargo build --release --locked --example native_search
python3 scripts/verify_translations.py
python3 scripts/compare_search.py --native target/release/examples/native_search
python3 scripts/verify_core_smoke.py --native target/release/examples/native_search
```

Fetch the source and executable pins using `scripts/blast_source.py fetch` and `scripts/blast_tools.py fetch`. The full-precision executable oracle requires Linux x86-64, GDB and permission to trace its child. The Rust library has no original BLAST C/C++ build dependency.

The repository release profile is ThinLTO with one codegen unit. A consumer application's root Cargo profile controls dependency code generation and must reproduce these settings for equivalent build measurements.

## Recorded benchmark result

The [corrected SEG-default benchmark](benchmarks/seg-default-c144e261/README.md) retains all raw reports: 54 measured pairs plus 12 warmup pairs. Rust/original full-workflow ratio is 0.8497, 95% interval [0.8420, 0.8574]: 15.0% less wall time. Five cases improve and one is inconclusive, with no supported regression. The workflow includes database construction, search and output: disk storage/two processes for original, memory storage/one process for Rust.

The supplementary view excludes original database creation while retaining Rust's full workflow: ratio 1.2108 [1.1980, 1.2237], or 21.1% more time for Rust. This favors original by excluding its build cost. Neither view measures engine-only speed. All source/build/input identities passed integrity checks; earlier reports remain separately available.

## Remaining completion gates

Actual Linux/Windows/macOS CI remains unverified. The workflow is configured, now fetches source before checking source hashes on Linux/macOS, and includes the default query-split and result-limit comparisons on Linux. This checkout has no Git remote or commits. Authenticated reads confirm that neither henriksson-lab/btyper_blast nor mahogny/btyper_blast exists. A repository destination was requested and is still needed before CI can run.

The [known isolated last-bit math difference](numerical-compatibility/README.md) is caused by reference glibc's extended exponential. The user's choice between documented platform-math tolerance and exact GNU bits remains pending. No global numerical equivalence claim or input-specific patch was made.

The corrected source passes all 1,259 real HSP raw comparisons, both default query-split boundary raw cases, the printed smoke and 12 protein differential workloads. Seven result-limit cases also match all 4,000 ordered HSPs at CLI precision, including original/Rust agreement on 550 multi-query CBS early terminations. These finite cases do not prove universal input parity. Memory database creation is available; disk/v5 persistence is excluded from the goal.
