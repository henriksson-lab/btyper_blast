# BLAST translation starting points

Start with selected kernels from `algo/blast/core`, preserving one upstream function per Rust function. The C core already separates search computation from result formatting and database storage. Its large shared files also contain unrelated programs, so port functions as they become necessary rather than entire files.

The reference is **NCBI BLAST+ 2.17.0+**, downloaded from NCBI and checked against its published MD5, with SHA256 pinned in [blast-source.json](blast-source.json). The complete extracted reference is under [the local toolkit root](../vendor/ncbi-blast-2.17.0+-src/c++/). Paths below are relative to that root. There are 556 mapped upstream functions/methods. Native memory database construction and database-wide standard search APIs are assembled. The 25-HSP smoke and 1,259 real-genome HSPs match the unmodified original core, including logical field order, seeds/scripts and exact E-value/bit-score bits. These finite workloads do not establish universal parity; actual cross-platform CI and final speed evidence remain pending. Disk/NCBI v5 persistence is outside scope. See [progress.md](progress.md) for current state; the dependency investigation below records the original starting scope.

## Current assembled boundary and evidence

The selected kernels now feed complete native memory searches. The current real-workload counts are ST14: 6 megablast, 481 blastn, 52 tblastn; S288C: 16, 573, 131. [benchmark-inputs.json](benchmark-inputs.json) pins the public genomes and derived-query provenance. Full-precision comparison uses the official unmodified executable, not a reimplemented statistics oracle: GDB captures `BLAST_ComputeTraceback_MT` returns after final sorting/pruning and before Seq-align clamps very small E-values. Every adaptive CLI batch is retained, with raw local indexes and explicit global mappings. Layouts come from compiled pinned headers.

Implementation checks exposed required API-layer details as well as kernel paths: C++ nucleotide diagonal-separation defaults affect which seed survives containment; the outer final result sort applies the original tiny-E-value equivalence and OID tie break; tblastn identity/positive counts use the masked core query, while displayed identities come from restored aligned columns. The three source target-translation-cache functions, scratch owner lifetimes and safe AA scan span/hit-run bounds are preserved with recorded adaptations. No unchecked indexing or alternate scan algorithm was introduced.

The repository release profile uses ThinLTO/one codegen unit; consuming Cargo applications choose their own root profile. Final paired benchmark reports are pending, with no speed claim yet. See [progress.md](progress.md) for current evidence and open gates. The remaining sections retain the dependency investigation and source boundaries that guided translation.

## Supported behavior

Both searches need a **nucleotide subject database**. `tblastn` accepts a protein query and translates nucleotide subjects in six frames. Start with a shared in-memory nucleotide database for small genomes; keep storage in memory; disk persistence is excluded. Search results always return as Rust values.

The executable named `blastn` defaults to task `megablast`. Task `blastn` is a different preset. Both presets are exposed explicitly because their lookup and gapped-extension paths differ. Do not accidentally implement traditional BLASTN and call it the executable default.

| Setting | Plain blastn command | blastn task | Default tblastn |
| --- | --- | --- | --- |
| Task | megablast | blastn | tblastn |
| Word size | 28 | 11 | 3 |
| Nucleotide reward and penalty | 1 and -2 | 2 and -3 | BLOSUM62 |
| Gap open and extension options | 0 and 0, invoking default greedy handling | 5 and 2 | 11 and 1 |
| Extension and traceback | Greedy | Dynamic programming | Dynamic programming with composition adjustment |
| Query filtering | DUST, soft masking | DUST, soft masking | SEG, hard masking |
| Expect threshold | 10 | 10 | 10 |
| Subject genetic code | Not applicable | Not applicable | 1 |

Sources: `src/algo/blast/api/blast_nucl_options.cpp`, `blast_prot_options.cpp`, `tblastn_options.cpp`; `include/algo/blast/core/blast_options.h`; and `src/algo/blast/blastinput/{blastn_args,tblastn_args,blast_args,cmdline_flags}.cpp`. CLI argument extraction can override API defaults, so both layers matter. These are distinguishing settings, not an exhaustive frozen option set. Pin all remaining defaults when implementing the option structs.

## First translation batch

Each `.c` normally has its public declarations in `include/algo/blast/core/<stem>.h`. Bring across only required declarations and constants, with private headers where needed. Line counts describe whole source files, not proposed translation size.

| Source under src/algo/blast/core | Whole file lines | Starting subset and dependency order |
| --- | ---: | --- |
| `blast_encoding.c` | 121 | Constant alphabet conversion tables and sentinels. Strongest first candidate: no runtime function dependencies. Retain exact byte values and alphabet distinctions. |
| `blast_program.c` | 108 | Predicates needed by supported program paths. Expand `SAFE_CAST_INT_TO_BOOLEAN` explicitly; preserve program bit masks. |
| `gapinfo.c` | 232 | Edit-script allocation and reset, then reallocation helper, add-new helper, and operation accumulation. Preserve operation ordering and run merging. Record allocation adaptations. |
| `blast_util.c` | 1413 | Select frame bookkeeping and packing helpers; translate `s_CodonToAA` before `BLAST_GetTranslation`. Reverse sequence handling and ambiguity rules are shared foundations. Do not translate this entire utility file. |
| `ncbi_math.c` | 494 | Small numerical functions only when a required statistics caller reaches them. Preserve formulas and evaluation order. |
| `boost_erf.c` | 256 | `ErfImpl`, then `Erf` and `ErfC` as required by statistics. `ErfImpl` is self-recursive. This is Boost licensed code. |

[translation-candidates.csv](translation-candidates.csv) records 556 scoped functions and methods with original source ranges, function hashes, and Rust identities. All currently listed definitions are translated; this is an inventory of selected implemented paths, not the complete toolkit. Direct calls and advisory order are available for the initial C core analysis; later C++ rows are recorded from selected source ranges and do not claim a parsed call graph. Packed/partial translation, nucleotide packing, ambiguity encoding, required RNG modes, and symmetric DUST are implemented. Unsupported program/option branches remain outside this inventory. [translations.json](translations.json) records implemented definitions and adaptations, including C++ class, constructor, and operator mappings. `blast_encoding.c` contains data rather than functions; its tables are also translated and checked against the C reference. Translated candidate rows identify their implemented Rust counterparts.

`blast_query_info.c` is a subsequent structural candidate, after the necessary pieces of `blast_util.c`. Its PHI pattern-copy and pattern-free paths do not belong to the supported searches. The genetic-code singleton is storage, not a translation algorithm: trace how code 1 is loaded before deciding which storage functions are necessary. Do not port generic dynamic arrays solely because the singleton uses them.

## Required search layers after the foundations

| Layer | Relevant sources | Scope constraint |
| --- | --- | --- |
| Default construction and thresholds | `blast_options.c`, `blast_parameters.c`, `blast_setup.c`, API preset files | Port required function paths and option structures; omit generic CLI and remote-service infrastructure. |
| Nucleotide masking | `src/algo/dustmask/symdust.cpp`, `include/algo/dustmask/symdust.hpp`, API `dust_filter.cpp` | Symmetric DUST is outside the C core. Translate the masker methods one-to-one; replace only object-manager input/output adaptation. |
| Protein masking | `blast_seg.c`, selected `blast_filter.c` functions | Needed for default protein query filtering and relevant subject masking during redo alignment. |
| Nucleotide seeds and scans | `blast_nalookup.c`, `blast_nascan.c`, `na_ungapped.c` | Follow dispatch for the supported preset and all accepted input sizes. Do not assume one scanner suffices. |
| Protein seeds and scans | `blast_aalookup.c`, `blast_aascan.c`, `aa_ungapped.c` | Needed by tblastn despite its nucleotide database. Keep protein neighborhood generation and threshold behavior. |
| Gapped extension | `blast_extend.c`, `blast_gapalign.c`, `greedy_align.c`, `gapinfo.c` | Megablast needs its greedy path; traditional BLASTN and tblastn use DP paths. |
| Scoring and statistics | `blast_stat.c`, `ncbi_math.c`, `boost_erf.c`, required matrix data | E-values depend on sequence composition and effective lengths. A fixed BLOSUM62 table alone is insufficient. |
| Default tblastn score adjustment | `blast_kappa.c`, required `src/algo/blast/composition_adjustment/*.c` functions | Composition mode 2 is enabled by default; do not disable it to shrink the port. Trace its actual callees rather than porting the directory wholesale. |
| Hits and traceback | `blast_hits.c`, `blast_itree.c`, `link_hsps.c`, `blast_traceback.c` | Preserve pruning, HSP linking, sorting, tie handling, and final score calculation. |
| Search orchestration | `blast_engine.c`, `lookup_wrap.c`, selected stream and collector functions | Assemble last. Preserve selected engine function structure; adapt library input and result ownership at the boundary. |

Do not initially implement PSI iterations, PHI pattern searches, RPS databases, translated queries, mapper/spliced alignment, remote BLAST, optional masking databases, culling/best-hit options, CLI argument parsing, or report/XML/ASN.1 result formatting. Shared support code remains eligible if a required default path actually calls it. Retaining enum discriminants and checking an unsupported-program condition does not require implementing that program.

Symmetric DUST and its default RNG path are now implemented, with 22,473 original-class comparisons. The native query adapter also matches 414 encoding/mask setups captured from actual BLAST archive requests. Query hyphens are ignored (the BLAST-specific FASTA reader uses fHyphensIgnoreAndWarn); database parsed gaps become N. The standalone dustmasker application's long-N preprocessing differs from BLAST's direct use of CSymDustMasker, so its output is not substituted for the query-mask oracle.

## In-memory results and databases

The C engine offers `Blast_RunPreliminarySearch`, `Blast_RunTracebackSearch`, and `Blast_RunFullSearch`. Its `BlastHSPResults` and edit scripts are an existing in-memory result boundary. Follow `src/algo/blast/api/local_blast.cpp` to understand preliminary search and traceback orchestration, without translating the whole C++ object-manager wrapper.

`blast_seqsrc.c` and `include/algo/blast/core/blast_seqsrc_impl.h` describe the sequence-source abstraction. A Rust storage adapter can provide sequence count, total nucleotide length, subject retrieval, iteration, and required encodings from memory. Track its callback dispatch manually: the comparator cannot prove those dependencies. Preserve database statistics for E-values and database-wide ranking when searching in-memory subjects; independent pairwise searches are not an equivalent replacement.

For small genomes, a database builder can keep identifiers, sequence offsets, packed nucleotide data, ambiguity information, and database statistics in memory. Retrieval needs both packed scanning data and restored ambiguity data for scoring, traceback, and translation. This storage boundary now feeds the native search APIs; upstream sequence-source and encoding semantics still govern the kernels.

## Modern on-disk database creation

NCBI's current builder is **makeblastdb**, replacing the old `formatdb`. This release defaults to **version 5**. A compatible on-disk writer involves more than LMDB:

1. `src/app/blastdb/makeblastdb.cpp` selects version 5 and calls `CBuildDatabase`.
2. `src/objtools/blast/seqdb_writer/build_db.cpp` handles input conversion and calls `CWriteDB`.
3. `writedb.cpp` and `writedb_impl.cpp` assemble records and volumes.
4. `writedb_volume.cpp` and `writedb_files.cpp` serialize sequence, index, and header files. `CWriteDB_File::Create` opens an actual file stream.
5. `writedb_lmdb.cpp` manages LMDB identifier and taxonomy indexes. LMDB is file backed; it does not provide a drop-in heap-only replacement for this writer.

`writedb_convert.cpp` is the best early database-code candidate, especially the raw-buffer overload of `WriteDB_Ncbi4naToBinary` and its ambiguity-region builder. Translate only these selected methods and helpers if NCBI-compatible packing is needed. The packed sequence has a length remainder in its final byte; ambiguity runs have a separate representation. Ambiguous bases are filled using `CRandom`, seeded from sequence length. Reproducing byte output requires the same random generator and call sequence, not an arbitrary Rust RNG. Its `CSeq_inst` overloads and FASTA/object-manager wrappers can wait.

That raw-buffer conversion, required ambiguity-builder methods, and seeded LFG methods are now translated. `DatabaseBuilder` uses them for in-memory subjects and retains the unambiguous input packing choice, including short-sequence ties. Exact packed data and ambiguity bytes match 377 subjects from official makeblastdb v5. The native FASTA adapter supports ordinary sequence input and returns diagnostics in memory; it does not implement the toolkit's complete Bioseq/gap metadata layer. No compatible v5 files are written by the library yet.

The user explicitly excluded NCBI v5 disk compatibility from this goal. The first in-memory database can feed the same search kernels without implementing ASN.1 headers, taxonomy indexes, aliases, volumes, or LMDB. If v5 compatibility is implemented, retain serialization routines one-to-one and isolate necessary I/O adaptations. Do not describe a custom in-memory database as an NCBI v5 database.

## Generated code boundaries

The source archive is a toolkit source subset, not a fully generated build tree. For example, `src/objects/blastdb/Makefile.in` declares `ASN_PROJ = blastdb`; `blastdb.asn` and `blastdb.def` drive generated metadata classes. The hand-maintained `Blast_def_line*.cpp` extensions alone do not define those classes. Database writer dependencies on `CBlast_def_line_set` and `CSeq_inst` therefore cannot be resolved by looking only for ordinary function definitions.

Keep generated serialization classes outside the first memory-database and result boundary. The excluded v5 header compatibility path would require ASN.1 serialization; it is not part of this implementation. The C core still needs toolkit configuration headers and macros for building the C comparison harness, even though its algorithm implementations are directly available.

## Comparator guidance and reproducibility

The local comparator's C analysis of `src/algo/blast/core` reported **1273 functions, 2125 edges, 0 ambiguous call sites, and 1731 unresolved call sites**. This includes headers and unsupported code. Missing edges include ordinary external calls and macros as well as runtime dispatch; this count is not a generated-code defect count.

To recreate the advisory analysis from the repository root:

```sh
CCC=/home/mahogny/github/claude/code-complexity-comparator/target/release/ccc-rs
CORE=vendor/ncbi-blast-2.17.0+-src/c++/src/algo/blast/core
"$CCC" analyze "$CORE" -l c --recurse -o /tmp/btyper-blast-core-analysis.json
"$CCC" order /tmp/btyper-blast-core-analysis.json -o /tmp/btyper-blast-core-order.csv
python3 scripts/translation_inventory.py /tmp/btyper-blast-core-analysis.json --order /tmp/btyper-blast-core-order.csv
```

Advisory order numbers refer to the full parsed core; they are not a required translation sequence. Review direct calls, macro bodies, callbacks, generated dependencies, and the actual preset before moving a candidate into the implementation scope. Hashes in [reference-files.json](reference-files.json) cover the reference files recorded for this investigation; the archive hash pins the complete release.

## Validation and portability

The crate must support Linux, Windows, and macOS. `.github/workflows/build.yml` checks formatting, all Cargo targets, and kernel tests in debug and release mode on all three systems. Linux and macOS additionally compile the pinned C sources and run the differential tests against fresh reference output. Native smoke searches run on all three configured platforms; Linux additionally compares full-precision original core results through GDB. Remote execution remains unverified. Preserve target ABI semantics for C `long`, use explicit byte order in serialized values, and guard architecture intrinsics. `BLAST_Nint`, for example, returns the target's `std::os::raw::c_long`. Record target-dependent behavior explicitly.

Use the pinned C code as the oracle for kernels and the pinned executable behavior for complete searches. Relevant upstream fixtures are in `src/algo/blast/unit_tests/api/`, including `gapinfo_unit_test.cpp`, `queryinfo_unit_test.cpp`, `ntlookup_unit_test.cpp`, `aalookup_unit_test.cpp`, `stat_unit_test.cpp`, `traceback_unit_test.cpp`, and `redoalignment_unit_test.cpp`. Translate only relevant test cases. NCBI v5 file compatibility also needs reader/writer round trips against NCBI tools; an in-memory search test cannot establish it.


## Nucleotide lookup dispatch traced for standard presets

`BlastChooseNaLookupTable` chooses small width 8 for traditional word 11 below 12,000 estimated entries and default megablast word 28 below 8,500 entries. It can instead choose the ordinary nucleotide table when the query offset or entry count exceeds the small table's 15-bit limits. Traditional word 11 then selects megablast width 10 below 180,000 entries or width 11 above that; default word 28 selects megablast width 11 below 300,000 entries or width 12 above that. `EstimateNumTableEntries` sums each included interval's `right - left`, without adding one or subtracting word length. Preserve that estimate rather than substituting a precise count.

For small width 8, `scan_step = word_length - lut_word_length + 1`: traditional blastn selects `s_BlastSmallNaScanSubject_8_4`, and default megablast selects `s_BlastSmallNaScanSubject_8_1Mod4` with stride 21. Both scanners and small-table construction are translated and compared against unchanged original C. The generic `Any` scanner is not a substitute for these selected paths. The native seed adapter now assembles this dispatch and the required larger-table alternatives.

Small-table finalization reserves overflow slots 0 and 1; backbone -1 means empty and values <= -2 are negative overflow cursors. Chains end in -1. The overflow limit includes these terminators: two chains may fail even if their individual hit counts fit in 15 bits. `s_SeqLocListInvert` deliberately omits leading excluded intervals of fewer than four bases and retains subsequent inverted intervals even when their left bound exceeds the right bound. These details are differential-tested.

`BlastCompressBlastnaSequence` is rolling compression for query extension, distinct from packed database NCBI2na. It allocates `length + 3` bytes, sets the sequence pointer three bytes into the allocation, and includes right-justified prefix and left-justified suffix bytes. Its signed C rolling accumulator can overflow on long queries; Rust uses wrapping shifts to retain the original low-byte machine behavior. The lookup C oracle uses `-fwrapv`. The unrolled stride-4 C scanner forms pointers before the allocation for its remainder entry; Rust preserves the entry/unroll using signed offsets and only dereferences valid indices. These are recorded adaptations, not algorithm substitutions.


The ordinary nucleotide table fallback is now translated, including `s_BlastNaLookupFinalize`, `BlastNaLookupTableNew`, mask-at-hash predicates, `s_BlastLookupGetNumHits`, `s_BlastLookupRetrieve`, both scanners, and `s_NaChooseScanSubject`. Its cells hold three direct offsets, or a cursor into 32-bit overflow storage; the presence vector gates access. Preserve this representation rather than replacing it with the small table's negative cursors. The ordinary scanner refuses a whole hit chain if it exceeds remaining capacity; unlike the small scanner, it does not reserve `longest_chain` up front. A zero or undersized output buffer can therefore return with the same scan offset.

The generic ordinary scanner has separate loops for aligned width >5 without subject masking, unaligned/masked width >5, and widths <=5. The no-masking aligned loop rounds start/end offsets down to byte boundaries, including unaligned starts; the other loops retain residue offsets. All paths are compared to original C with readable tail padding in subject fixtures. Standard presets only select ordinary width 8; width 4..7 cases verify the shared function's remaining branches without adding public search presets. Query input uses inclusive ordered intervals and can retain high offsets despite a low unmasked entry count. The 100,128-base test query verifies standard lookup selection and retains offsets beyond 65,535 without truncating them to 16 bits.


## Contiguous megablast scanners selected by the presets

For traditional blastn word 11, the larger-query table uses width 10/stride 2 (`s_MBScanSubject_10_2`) or width 11/stride 1 (`s_MBScanSubject_11_1Mod4`). For default megablast word 28, it uses width 11/stride 18 (`s_MBScanSubject_11_2Mod4`) or width 12/stride 17 (`s_MBScanSubject_Any`). These four scanners plus `s_BlastMBLookupHasHits` and `s_BlastMBLookupRetrieve` are translated and checked with unchanged original C. The generic function's inherited width-16 branches are outside supported scope; no width-16 validation is claimed. In particular, original aligned width-16 inputs can reach an invalid negative shift, so they are not a supported domain.

Megablast's presence vector may compress several table entries into each bit: the vector word uses `index >> pv_array_bts`, while the bit still uses `index & 31`. A positive presence test can therefore retrieve zero hits. The hit buffer check is `total_hits >= max_hits - longest_chain`, performed on a positive presence test before retrieval. Retain this even for a false positive; substituting ordinary-table checks changes restart positions. Heads and next-position links use one-based query offsets, with zero as terminator; returned offsets subtract one and traverse latest occurrence first. The scanner fixture tests both compressed-PV false positives and original chain order.

The MB fixture now calls original `BlastMBLookupTableNew` directly and compares full Rust construction: table heads, links, PV sizes/bits, masks, dimensions and longest-chain estimate. It then compares scans on the Rust-constructed table. The source dependency `blast_seqsrc.c` is compiled only to resolve constructor references; standard options never take its database-counting branch. No generated object-manager/toolkit classes are substituted. The Rust constructor explicitly guards unsupported discontiguous, mapper, database-count-filtered and width >12 paths rather than translating their dependencies.


`LookupTableWrapInit_MT` estimates included entries, chooses the table, then retries ordinary-table construction at the same lookup width if small-table finalization returns nonzero. The native lookup adapter retains this overflow fallback alongside the direct max-offset fallback selected by `BlastChooseNaLookupTable`. Default `LookupTableOptionsNew(eBlastTypeBlastn)` starts with word 28 and zero-initialized `db_filter` and `mb_template_length`, and explicitly sets `stride = 0`. The standard traditional task changes word length to 11; it still uses contiguous lookup without database word-count filtering. API assembly also retains the C++ preset overrides, including minimum diagonal separation 50 for traditional blastn and 6 for megablast.


## Native standard-preset seed adapter

`seed::NucleotideSeedIndex::from_segments` accepts assembled BLASTna contexts and included regions, retaining full-context offsets. It follows the nucleotide factory branch, selecting all required table widths and preserving the small-table finalization fallback. It accepts empty `[right+1,right]` intervals produced by mask setup and rejects overlapping/out-of-range inputs. The public query batch adapter now performs context assembly and exact mask complement wiring before constructing the seed index.

The seed fixture executes original lookup defaults, selection, constructors, callback dispatch and scanning; its test driver supplies the nucleotide wrapper orchestration without executing the full multi-program `LookupTableWrapInit_MT`. Raw seed pairs and each restart offset are compared to the Rust adapter. Query sizes hit the 8,500/12,000/180,000/300,000 selection thresholds and the 800,000 PV boundary; high-offset masked queries select ordinary tables. Subjects of 1..129 bases exercise short subjects, remainder bytes and packed-window tail reads. The memory database stores private zero tail padding while its public packed getter exposes the exact original bytes.

A seed pair is not a reported alignment. The native search APIs now connect it to full-word validation, ungapped/gapped extension, traceback, score/E-value calculation, filtering and result collection. The included-region input is explicit at this stage; the native seed adapter does not recompute or invent query masks.


## Query-context assembly and mask complement

`PreparedNucleotideBatch` follows the both-strand nucleotide branches of `SetupQueryInfo_OMF` and `SetupQueries`, using translated `BlastQueryInfoNew`, `s_QueryInfo_SetContext`, `QueryInfo_GetSeqBufLen`, single-strand conversion and both-strand stitching. Positive contexts consume their length plus one separator; zero-length contexts in the general core setter consume no separator and invalidate the context. Reassigning a positive length does not restore validity. The final sequence allocation includes leading/trailing sentinels plus an extra zero byte; the sequence-block view excludes the endpoint sentinels and retains all internal offsets.

Mask coordinates for both nucleotide contexts enter in forward-query coordinates, even though the negative context's residues are reverse complemented. `BLAST_ComplementMaskLocations` reverses the negative mask list in place and then computes `end_offset - right` / `end_offset - left`. Passing the query adapter's already reversed masks would reverse coordinates twice. The Rust batch therefore copies forward DUST masks into both context entries and invokes the original complement algorithm once. Fully masked contexts retain the original `[end+1,end]` markers. The native seed adapter accepts those markers and preserves their effect on entry estimates/maximum offsets.

The mixed C/C++ reference compiles actual C query-info/filter functions and extracts unchanged C++ setup helper definitions. Test-only expanded-strand data/ownership interfaces avoid generated object-manager dependencies; the Na-strand enum is generated from the pinned ASN.1 values and ESentinelType is extracted from the original header. Original symmetric DUST computes batch masks. The fixture compares all context fields, strand bytes, defined allocation prefixes, mask complements and their list mutation, plus both standard lookup selections. It does not replace full official-program result comparisons.


## Nucleotide initial-hit extension

`BlastChooseNaExtend` selects the one-byte aligned small-table extension for traditional word 11 / lookup width 8 / stride 4, and the general small-table extension for default word 28 / width 8 / stride 21. Larger MB width 10 / word 11 and width 11 / word 28 use the general nucleotide extender; equal word/lookup widths select the direct path. Ordinary width 8 / traditional word 11 uses the aligned ordinary extender. These mini-extension callbacks and their selector are translated, together with the default unmasked-database word-finder loop. Original-C differential tests cover the assembled loop and state across subjects; score cutoff inputs are still explicit.

`BlastScoreBlkNuclMatrixCreate` computes ambiguity scores using the original ordered degeneracy formula and BLAST_Nint, not a simplified match/mismatch matrix. Its sentinel row/column is INT4_MIN/2. `s_NuclUngappedExtend` first scores four-base blocks through a 256-entry table, then calls `s_NuclUngappedExtendExact` when the approximate score is >= the reduced cutoff. Keep both loops, the raw BLASTna byte packing, strict improvement/tie rules, right-side dynamic dropoff and the approximate known-match extent. Running exact extension unconditionally would change both behavior and cost.

`s_IsSeedMasked` deliberately checks lookup membership through the stored callback instead of traversing masked intervals. Its packed fetch switch avoids unnecessary tail reads at certain alignments. `s_TypeOfWord` can reposition q_off/s_off before rejecting a word; it checks the right edge, searches a valid left edge, checks interior seeds and optionally tests a double word seed-first and then base-first. All output mutations, including rejection cases, are differential-tested. The query sequence argument is unused in the original body and omitted in the Rust binding.

The diagonal hash uses 512 buckets and 16-byte cells, with a signed 31-bit level and one-bit saved flag. Insertion checks a matching diagonal first, then reuses the first stale collision, and allocates a new head only when no cell was reused. Do not replace it with a generic Rust map, perform a full-chain matching search before stale reuse, or widen the packed fields into separate words. The Rust binding retains the row size and logical bitfield conversions. Original `BlastInitialWordParametersNew` selects this hash when total concatenated context length exceeds 8,000; `Blast_ExtendWordExit` clears/rebases state at INT4_MAX/4. The native setup now wires those constructors/lifetime calls and initial-hit collection.

The extension fixture compiles unchanged C files and directly calls the original file-local functions. It supplies four-base score tables as inputs using the block in the original parameter constructor; that fixture alone does not validate the complete parameter constructor or default cutoff assembly. Those paths are now assembled and exercised by whole-search comparisons. A separate word-finder fixture now directly verifies default unmasked-database orchestration with explicit cutoff inputs. Primitive extension equality is insufficient evidence for final alignment/result parity.
