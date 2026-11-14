# Recorded workflow baseline

These reports measure artifact `5ae87386bbbd1db8191c5a9ccb8880eff0b650f96193f0ea623da85b8c5b9366`, before subsequent nucleotide buffer-lifetime and query-splitting restoration. They are retained as a historical baseline, not measurements of the evolving working tree.

The six-case suite used 54 measured pairs and 12 warmup pairs on two pinned public genomes. Native/reference total mean wall time ratio was 0.87178 (95% paired bootstrap interval 0.86434–0.87920): 12.8% less workflow wall time. Five cases improved; yeast tblastn was inconclusive. No case had a statistically supported regression. This compares complete database creation plus search/output workflows: original uses two processes and disk v5 storage, Rust one process and memory storage. It does not isolate engine speed.

`summary.json` retains all raw reports, paired samples, identities and methods; individual reports and build metadata are also retained. New source changes require new measurements before attributing this result to them. Cross-platform CI remains unverified.
