# BLAST Rust library

A scoped Rust translation of **NCBI BLAST+ 2.17.0+** for standard megablast, traditional blastn and tblastn searches on nucleotide databases. Database construction and search results stay in memory. No original BLAST executable or C/C++ library is invoked by the Rust search code.

The search pipelines include DUST/SEG filtering, the original lookup/scanning and gapped algorithms, traceback, ambiguity handling, CBS2 composition adjustment, linked statistics, preliminary candidate selection, ranking and result limits. Long subjects and translated frames retain the original chunking and HSP merging. Only dependencies needed by these presets are translated; the toolkit's other programs and optional modes are outside this scope.

## Use as a library

```rust
use btyper_blast::{
    DatabaseBuilder, NucleotideTask, PreparedNucleotideBatch,
    PreparedNucleotideQuery, search_nucleotide_database,
};

let sequence = b"ACGTTGCAGATCGTACGATTCGAGCTAGCTGA";
let mut builder = DatabaseBuilder::new();
builder.add_sequence("contig_1", sequence)?;
let database = builder.build()?;
let query = PreparedNucleotideQuery::new("query_1", sequence, NucleotideTask::Blastn)?;
let batch = PreparedNucleotideBatch::new(&[query])?;
let results = search_nucleotide_database(&batch, &database)?;
```

`NucleotideTask::Megablast` selects the plain `blastn` command's default task. `NucleotideTask::Blastn` selects its traditional task. For tblastn, construct `PreparedProteinQuery` values and a `PreparedProteinBatch`, then call `search_protein_database(&batch, &database)`. The default uses BLOSUM62, gaps 11/1, SEG, composition mode 2 and genetic code 1.

`NucleotideDatabase::from_fasta(reader)` accepts ordinary nucleotide FASTA and exposes omitted invalid-residue diagnostics through `warnings()`. Direct builders reject invalid residues. Database U becomes T and '-' becomes N; query hyphens are ignored. Identifiers, insertion order, duplicate identifiers and database statistics are retained. Empty database records are rejected. Full NCBI Bioseq/FASTA extensions are not implemented. Database persistence is outside scope; the library creates databases in memory.

Results contain owned HSPs, native context coordinates, edit scripts and full-precision statistics. `report::AlignmentReport` provides oriented endpoints, original alignment strings and BTOP. For tblastn, core identities/positives use the SEG-masked query as in upstream C; displayed identities are recomputed from restored aligned residues as in the original formatter. Those counts can legitimately differ.

## Correctness evidence

There are **596 mapped upstream functions/methods**, with pinned source ranges, hashes and necessary adaptations. The synthetic smoke's **25 HSPs** and **1,259 HSPs on two public genomes** match the unmodified official executable's core results: ordered query/subject/HSP lists, scores, identities/positives, contexts, frames, seeds, edit scripts and exact IEEE754 bit patterns for E-values and bit scores.

| Public genome workload | Megablast | Traditional blastn | Tblastn |
| --- | ---: | ---: | ---: |
| *Staphylococcus epidermidis* ST14, 2.48 Mb | 6 | 481 | 52 |
| *Saccharomyces cerevisiae* S288C, 12.16 Mb | 16 | 573 | 131 |

Seven default-limit cases also match 4,000 ordered HSPs at CLI precision, including coupled multi-query CBS early termination; see [the corrected query-default evidence](upstream/validation/seg-default-c144e261/README.md). A [known isolated ErfC difference](upstream/numerical-compatibility/README.md) of one final floating-point ULP is caused by the reference math library; universal bitwise numerical parity is not claimed.

The genome inputs and deterministic derived queries are pinned in [benchmark-inputs.json](upstream/benchmark-inputs.json). These are finite workload results, not a claim of equivalence for every possible input. Focused differential suites also cover 12 nucleotide workloads with both presets and 12 protein workloads. Earlier kernel and database conversion checks remain available.

The core oracle reads the official child process through GDB before Seq-align clamps E-values below `1e-180`. It derives ABI layouts from the pinned headers and captures every final query batch. Global indexes have explicit batch mappings; original local fields, numeric values, seeds and scripts are preserved. Tabular comparison separately checks original CLI formatting.

## Build and reproduce

```sh
cargo check --all-targets --locked
cargo test --locked
cargo test --release --locked
cargo build --release --locked --example native_search
python3 scripts/verify_translations.py
python3 scripts/compare_search.py --native target/release/examples/native_search
```

Fetch references with `python3 scripts/blast_source.py fetch` and `python3 scripts/blast_tools.py fetch`. The system-installed BLAST is not the oracle. Python 3.12+ is used for safe source extraction. The original executables and full-precision debugger oracle run on Linux x86-64; child-process tracing must be available.

```sh
python3 scripts/verify_core_smoke.py --native target/release/examples/native_search
python3 scripts/differential_nucleotide.py --native target/release/examples/native_search
python3 scripts/differential_protein.py --native target/release/examples/native_search
python3 scripts/differential_query_split.py --native target/release/examples/native_search
python3 scripts/differential_result_limits.py --native target/release/examples/native_search
python3 scripts/small_genome_inputs.py --fetch
```

`native_search [--raw-json] <megablast|blastn|tblastn> <database.fna> <queries.fasta>` is a validation/display example. Its results come from the native memory API; it does not invoke original BLAST. `--raw-json` exposes core values and float bits. [compare_core_results.py](scripts/compare_core_results.py) compares them with [core_result_oracle.py](scripts/core_result_oracle.py).

## Performance and portability

The [recorded six-case workflow benchmark](upstream/benchmarks/seg-default-c144e261/README.md) uses **15.0% less wall time** than original BLAST (95% paired interval: 14.3–15.8% less). Five cases improve; yeast tblastn is inconclusive. It includes database construction and output, with disk storage/two processes for original and memory storage/one process for Rust. All paired observations and artifact/source/build/input identities are retained.

A supplementary view excludes original database creation while retaining Rust's full workflow. Rust then takes **21.1% more time** (interval: 19.8–22.4% more). The original database is already created and warm, while Rust still parses and builds its database. Neither comparison isolates engine speed.

Repository release builds use ThinLTO and one codegen unit. Cargo applications consuming this library must set their own root release profile to reproduce those settings; dependency profiles do not propagate.

CI is configured for Linux, Windows and macOS, including native smoke checks, with an additional full-precision original-core comparison on Linux. Actual remote CI runs remain unverified. The goal remains open: actual cross-platform CI and the numerical compatibility policy are unresolved. The finite parity and workflow-speed results are documented above.

- [Current progress and remaining gates](upstream/progress.md)
- [Source investigation and generated-code boundaries](upstream/investigation.md)
- [Translation policy](upstream/translation-policy.md), [function ledger](upstream/translations.json) and [data provenance](upstream/translation-data.json)
- [Pinned upstream source](upstream/blast-source.json) and [CI configuration](.github/workflows/build.yml)
