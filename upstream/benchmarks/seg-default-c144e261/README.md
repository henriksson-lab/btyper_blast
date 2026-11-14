# Corrected SEG default: workflow benchmark

[build.json](build.json) identifies artifact `c144e261` and its exact source/build fingerprints. The source is pinned BLAST+ 2.17.0+; the Rust release uses ThinLTO and one codegen unit. Both programs use one thread on CPU 2, with warm caches. Every measured output passed parity checks. Source, executable and input identity checks passed throughout measurement.

Nine measured pairs and two declared warmup pairs per case give 54 measured and 12 warmup pairs. Each pair randomizes original/Rust order. All trials, including outliers, remain in the reports. The primary metric sums six case arithmetic means and takes the Rust/original ratio; a deterministic 10,000-draw joint paired bootstrap gives its 95% interval.

| Dataset | Preset | Original mean s | Rust mean s | Rust/original [95% interval] | Classification |
| --- | --- | ---: | ---: | --- | --- |
| staph_st14 | megablast | 0.0730 | 0.0446 | 0.6105 [0.5966, 0.6263] | win |
| staph_st14 | blastn | 0.2280 | 0.1595 | 0.6997 [0.6873, 0.7146] | win |
| staph_st14 | tblastn | 0.3112 | 0.2998 | 0.9633 [0.9424, 0.9876] | win |
| yeast_s288c | megablast | 0.2053 | 0.1297 | 0.6318 [0.6228, 0.6389] | win |
| yeast_s288c | blastn | 0.5138 | 0.3649 | 0.7103 [0.6941, 0.7247] | win |
| yeast_s288c | tblastn | 0.8640 | 0.8668 | 1.0033 [0.9896, 1.0167] | inconclusive |

The suite ratio is **0.8497 [0.8420, 0.8574]**: **15.0% less wall time**, with an interval of 14.3–15.8% less. Five cases improve and yeast tblastn is inconclusive. No case has a supported regression in this full-workflow view. Earlier measurements remain in their original directories; these corrected-source results are not selected from multiple timing attempts.

The clock covers FASTA database construction, search and complete 19-field output capture. Original BLAST uses two processes and fresh disk v5 storage; Rust uses one process and memory storage. Directory creation/removal, executable verification, parity comparisons and reporting are outside the clock. These are the fixed public-genome workloads in [benchmark-inputs.json](../../benchmark-inputs.json), with deterministic derived queries. This does not establish engine-only or universal speed.

## Original search phase with database construction excluded

[reference-search-phase.json](reference-search-phase.json) retains a supplementary calculation from the same observations. Rust's complete workflow divided by the original search phase is **1.2108 [1.1980, 1.2237]**: about **21.1% more wall time**, interval 19.8–22.4% more. Five cases regress in this view; staph traditional blastn improves. The original database has already been created and is warm, while Rust still pays parsing, encoding and construction. This is not an independently measured prebuilt-database benchmark or an engine-only comparison. It shows why the primary improvement must not be described as proof of faster search kernels.

The combined source includes original query splitting and buffer lifetimes, partial translation caching, checked AA span hoists and the corrected source query-SEG overlap setting. This measurement does not isolate the causal contribution of each change. Actual remote Linux/Windows/macOS CI remains unverified.
