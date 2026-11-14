# Active goal and execution policy

## Objective

Complete a production-usable Rust library for blastn with standard settings and tblastn with standard settings, plus in-memory nucleotide database creation for small genomes. Use the pinned NCBI BLAST+ 2.17.0+ source and preserve one-to-one upstream function correspondence with systematic Rust snake_case names initially; translate only functions and dependencies required by the supported paths, recording necessary adaptations and retaining optimized upstream behavior. Return search results in memory. Resolve and document the distinction between the blastn executable's default megablast task and the traditional blastn task; make the supported presets explicit. Build and test on Linux, Windows, and macOS through CI. Establish reproducible differential tests against original BLAST+ 2.17.0+ covering hits, coordinates, strands, translation frames, edit scripts, raw scores, bit scores, E-values, ranking, filtering, and result limits, including realistic small genomes and adversarial edge cases. Establish a reproducible representative small-genome benchmark suite and fair original-BLAST baseline, compare equivalent search/database workloads with matched settings, resources, and output work, and optimize only after correctness parity is established. The goal is complete only when working library APIs and runnable validation/benchmark tools demonstrate equivalent results to original BLAST for the supported scope, passing cross-platform CI, and a repeatable statistically supported speed improvement over original BLAST on the documented representative benchmark suite, with regressions and limits reported honestly. Scaffolding, partial translations, unverified parity, or unmeasured speed claims do not satisfy completion.

## Execution policy (user update)

Use parallel workers whenever independent tasks make that useful. While required translation is unfinished, focus on getting the missing code written and the pipeline assembled, rather than adding large test corpora. Retain focused correctness checks and source correspondence; perform comprehensive differential validation and performance measurements before claiming completion.

The goal remains active. This file records the user's execution-policy amendment because the goal API exposes status updates but does not expose an operation to edit an active objective's wording. The original completion criteria are unchanged.

## User scope and fidelity amendments

The user explicitly omitted disk/NCBI v5 persistence. Database creation and search results remain in memory. Original makeblastdb v5 use in validation/benchmark tools is a reference workload, not a disk-storage feature of the Rust library.

The performance comparison must first establish that upstream optimizations have been transferred: memory lifetimes/layouts, original arenas and pools, scan unrolling, bit tricks, prefetching, SIMD and architecture dispatch where selected by the supported paths. Restore missing upstream mechanisms before considering bumpalo or a new arena design. Any later allocator change requires profiling evidence and equivalent results.

These amendments supersede the corresponding storage/execution wording in the original goal API objective; that API exposes status changes but cannot edit an active objective. The result-parity, cross-platform CI and measured-speed completion gates remain active.
