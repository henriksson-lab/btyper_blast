# Source lifetime restoration: recorded workflow benchmark

The artifact and exact source/build fingerprints are in build.json. Source is pinned BLAST+2.17.0+; Rust release uses ThinLTO/one codegen unit, single thread, CPU2 affinity, warm caches. Every measured output passed original parity checks. No source/executable/input integrity change occurred during measurement.

Nine measured pairs and two declared warmup pairs per case give54measured and12warmup pairs. Each pair randomizes original/native order. All measured trials, including outliers, are retained in the individual reports and embedded in summary.json. The primary metric sums six case arithmetic means and takes the native/reference ratio. A deterministic10,000-draw joint paired bootstrap gives its95% interval.

| Dataset | Preset | Reference mean s | Native mean s | Native/reference [95% interval] | Classification |
| --- | --- | ---: | ---: | --- | --- |
| staph_st14 | megablast | 0.0862 | 0.0479 | 0.5560 [0.5392, 0.5735] | win |
| staph_st14 | blastn | 0.2384 | 0.1554 | 0.6517 [0.6371, 0.6658] | win |
| staph_st14 | tblastn | 0.3100 | 0.2801 | 0.9035 [0.8850, 0.9257] | win |
| yeast_s288c | megablast | 0.2150 | 0.1360 | 0.6326 [0.6247, 0.6404] | win |
| yeast_s288c | blastn | 0.5079 | 0.3708 | 0.7300 [0.7252, 0.7352] | win |
| yeast_s288c | tblastn | 0.8920 | 0.8943 | 1.0027 [0.9857, 1.0205] | inconclusive |

Suite ratio **0.83776 [0.82988,0.84682]**: **16.2% less wall time**, interval15.3–17.0% less, or1.194x throughput at the fixed workload. Five cases improve and one is inconclusive; no supported regression. The yeast tblastn interval includes both a small improvement and a small regression.

This is FASTA database creation plus search and complete19-field output capture. Original uses two processes and a fresh disk v5 database; Rust uses one process and an in-memory database. Temporary-directory creation/removal, executable verification, parity comparison and reporting are outside the clock. No engine-only or universal speed claim follows. These public genomes/derived queries are the fixed suite in upstream/benchmark-inputs.json, not all possible bacterial/yeast workloads.

The source restore includes query working-set splitting, original partial translation caching, nucleotide/protein reusable scratch, OID-batch nucleotide expansion, persistent CBS matrix rows and checked AA span hoists. This combined comparison does not isolate each change's causal contribution. Earlier reports are retained under baseline-5ae87386. Remote Linux/Windows/macOS CI remains unverified.

## Original search-phase comparison

[reference-search-phase.json](reference-search-phase.json) retains a supplementary calculation from the same paired observations. Excluding the original database-build phase while retaining Rust's complete workflow gives a suite ratio of **1.21954 [1.20889, 1.23104]**: Rust takes about 22% more time in this view. Five cases regress; staph traditional blastn improves. The original database is already created and warm, while Rust still pays parsing, encoding and construction. This is neither an independently measured prebuilt-database benchmark nor an engine-only comparison. It prevents interpreting the full-workflow improvement as proof of faster search kernels.
