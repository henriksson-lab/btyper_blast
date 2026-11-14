# Working instructions

- Use parallel workers whenever independent translation or integration tasks can make useful progress concurrently. Give each worker clear ownership of files to avoid conflicting edits.
- While required translation remains unfinished, prioritize writing the missing code and assembling the search pipeline. Use focused checks for translation errors; defer broad fixture expansion and exhaustive test matrices until the pipeline is assembled.
- Preserve the pinned BLAST+ 2.17.0+ version, one-to-one function correspondence, systematic Rust snake_case names, optimized upstream algorithms, and the narrow dependency scope.
- Keep source provenance and necessary adaptations recorded. Do not claim parity or speed from compilation or isolated kernel tests.
- The completion criteria remain full supported-search result parity, Linux/Windows/macOS CI, and a demonstrated speed improvement. See upstream/goal.md for the objective and execution policy.

- Database creation is in memory. The user explicitly omitted disk/NCBI v5 persistence; do not implement it as a goal requirement.
- Before introducing bumpalo or a new arena design, audit and transfer upstream optimizations on the supported paths, including memory lifetimes/layouts, existing pools/arenas, unrolling, bit tricks, prefetching and SIMD/CPU dispatch. Restore those mechanisms first; consider new allocation strategies only after that audit and profiling evidence.
