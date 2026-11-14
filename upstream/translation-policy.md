# BLAST Rust translation policy

The reference is BLAST+ 2.17.0+, pinned by archive SHA256 in [blast-source.json](blast-source.json). The initial scope is nucleotide databases, `blastn` with standard settings, and `tblastn` with standard settings. Search results belong in memory as Rust values. Database persistence is a separate concern.

## Function correspondence

Translate only functions needed by these paths, together with their required callees, types, constants, macros, and data tables. A file appearing in a dependency report does not authorize translating all of it. Candidate records are investigation results, not claims that the dependency closure is complete. Unsupported programs must be rejected at the public boundary.

Initially preserve a one-to-one relationship between upstream and Rust functions. Do not combine, split, simplify, replace, or redesign algorithm functions. Retain upstream fast paths, specialized scanning variants selected by supported inputs, table values, loop and branch order, tie breaking, sentinels, and score calculations. Standard settings still choose different implementations according to sequence length and content.

Map function names deterministically to snake case:

1. Preserve prefixes, including `BLAST`, `Blast`, `s_`, and `WriteDB`; do not strip them.
2. Insert an underscore between an acronym and a following capitalized word, and between a lowercase letter or digit and a following uppercase letter.
3. Lowercase the result. Preserve digits and existing underscores.
4. Use source module identity to distinguish equal names in different files. Detect collisions within a module; resolve them explicitly in the mapping rather than silently renaming.
5. For C++ methods retain the class identity as the Rust type, map the method name using the same rule, and record the qualified original name and overload signature. Constructors, destructors, and overloads need explicit mappings when selected.

Examples: `BLAST_GetTranslation` becomes `blast_get_translation`; `GapEditScriptNew` becomes `gap_edit_script_new`; `s_CodonToAA` becomes `s_codon_to_aa`; `WriteDB_Ncbi4naToBinary` becomes `write_db_ncbi4na_to_binary`.

Keep upstream fixed numeric widths and signedness. Make C promotions, casts, and defined wrapping explicit when needed. Preserve floating point expression order; do not substitute other math routines or enable reassociation without equivalence evidence. C undefined behavior is not a behavior to reproduce: record any correction. Layout, pointer ownership, output parameters, allocation, callbacks, and synchronization may require Rust adaptations; record each departure beside its mapping. Rust naming alone does not justify changes in algorithm structure.

Every translated function needs its release, original path and symbol, original source range or hash, Rust path and symbol, translation status, inclusion reason, and any adaptation. Keep NCBI notices and author attribution; retain the separate Boost license for `boost_erf.c` and licenses of any other selected third party code.

## Dependency guidance

Use `/home/mahogny/github/claude/code-complexity-comparator` for suggested bottom-up order and later comparisons. Restrict its output to the reviewed function scope. Its name-based graph does not establish macro expansion, generated ASN.1 definitions, function-pointer dispatch, or runtime selection. An unresolved edge is an investigation item, not proof of a leaf. Recursive groups must be translated together when needed.

## Verification before optimization

Compare each translated kernel against the pinned C implementation. Then compare end-to-end default searches against BLAST+ 2.17.0+: hits, coordinates, strands and frames, edit scripts, raw scores, bit scores, E-values, ordering, and truncation. Include ambiguity symbols, both strands, all subject frames, low complexity, and inputs selecting alternate scan kernels. Performance work comes after behavioral agreement; a close source translation alone does not guarantee equivalent machine code or speed.


The nucleotide lookup translations preserve thin-chain allocation thresholds, signed 16-bit final-table encoding, overflow limits, scan unrolling, alignment-phase entry points, strict buffer comparisons, and exact restart offsets. Linked query intervals use ordered slices; owned arrays replace allocator/destructor calls. Rolling signed shifts use explicit wrapping machine arithmetic, with a matching `-fwrapv` C oracle. Pre-allocation C pointers used by the stride-4 remainder entry become signed offsets so Rust only dereferences valid array indices. All deviations are per-function ledger entries. Table selection retains unsupported-mode branches without treating those branches as supported search APIs.


Mode-specific constructors may translate only branches required by the explicit supported presets, provided unsupported modes are rejected before construction and every omitted branch is recorded in the function ledger. `BlastMBLookupTableNew` is such a scoped translation: contiguous, nonmapping, no database word-count filtering, widths 9..12. Its supported allocation/PV/filling path retains the original formulas and call order. Native ownership/dispatch adapters (database, query, seed) remain separate from the per-function kernel ledger; they are compared to original tooling or core orchestration, and do not imply complete search parity.
