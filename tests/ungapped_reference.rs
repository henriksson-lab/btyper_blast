use btyper_blast::core::blast_diagnostics::*;
use btyper_blast::core::blast_extend::*;
use btyper_blast::core::blast_lookup::SSeqRange;
use btyper_blast::core::blast_nalookup::*;
use btyper_blast::core::blast_nascan::*;
use btyper_blast::core::blast_options::BlastInitialWordOptions;
use btyper_blast::core::blast_options::{QuerySetUpOptions, SBlastFilterOptions};
use btyper_blast::core::blast_parameters::*;
use btyper_blast::core::blast_program::EBlastProgramType;
use btyper_blast::core::blast_query_info::{BlastContextInfo, BlastQueryInfo};
use btyper_blast::core::blast_stat::*;
use btyper_blast::core::masksubj::*;
use btyper_blast::core::na_ungapped::*;

fn number<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> i32 {
    tokens.next().unwrap().parse().unwrap()
}
fn bytes(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}
fn data<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> BlastUngappedData {
    BlastUngappedData {
        q_start: number(tokens),
        s_start: number(tokens),
        length: number(tokens),
        score: number(tokens),
    }
}
enum Table {
    Small(BlastSmallNaLookupTable),
    Ordinary(BlastNaLookupTable),
    Megablast(BlastMBLookupTable),
}
impl Table {
    fn wrap(&self) -> NucleotideLookupWrap<'_> {
        match self {
            Self::Small(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Small(t),
                lookup_callback: s_small_na_lookup,
            },
            Self::Ordinary(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Ordinary(t),
                lookup_callback: s_na_lookup,
            },
            Self::Megablast(t) => NucleotideLookupWrap {
                lut: NucleotideLookup::Megablast(t),
                lookup_callback: s_mb_lookup,
            },
        }
    }
    fn masked(&self) -> &[SSeqRange] {
        match self {
            Self::Small(t) => &t.masked_locations,
            Self::Ordinary(t) => &t.masked_locations,
            Self::Megablast(t) => &t.masked_locations,
        }
    }
}
fn verify(fixture: &str) {
    let mut hash = None;
    let mut extension = BlastExtendWord::default();
    let mut hitlist = blast_init_hit_list_new();
    let mut params = None;
    let mut seeds = Vec::new();
    let mut mini_id = 0usize;
    let mut finder_stats = BlastUngappedStats::default();
    let mini_callbacks: [TNaExtendFunction; 5] = [
        s_blast_na_extend_direct,
        s_blast_na_extend,
        s_blast_na_extend_aligned,
        s_blast_small_na_extend_aligned_one_byte,
        s_blast_small_na_extend,
    ];
    let mut matrix = [[0; 16]; 16];
    let mut score_table = [0; 256];
    let mut query = Vec::new();
    let mut subject = Vec::new();
    let mut subject_length = 0;
    let mut lookup = None;
    let mut word = 0;
    let mut width = 0;
    let info = BlastQueryInfo {
        last_context: 1,
        num_queries: 1,
        min_length: 128,
        max_length: 128,
        contexts: vec![
            BlastContextInfo {
                query_length: 128,
                ..Default::default()
            },
            BlastContextInfo {
                query_offset: 129,
                query_length: 128,
                ..Default::default()
            },
        ],
        ..Default::default()
    };
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "matrix" => {
                let reward = number(&mut tokens);
                let penalty = number(&mut tokens);
                let mut actual = [[12345; 16]; 16];
                assert_eq!(
                    blast_score_blk_nucl_matrix_create(reward, penalty, &mut actual),
                    0
                );
                for r in &mut matrix {
                    for value in r {
                        *value = number(&mut tokens);
                    }
                }
                assert_eq!(actual, matrix, "line {line}");
            }
            "score_table" => {
                for score in &mut score_table {
                    *score = number(&mut tokens);
                }
            }
            "sequences" => {
                query = bytes(tokens.next().unwrap());
                subject_length = number(&mut tokens);
                subject = bytes(tokens.next().unwrap());
            }
            "exact" => {
                let qoff = number(&mut tokens);
                let soff = number(&mut tokens);
                let x = number(&mut tokens);
                let expected = data(&mut tokens);
                let mut actual = BlastUngappedData {
                    q_start: -997,
                    s_start: -998,
                    length: -999,
                    score: -1000,
                };
                s_nucl_ungapped_extend_exact(
                    &query,
                    &subject,
                    subject_length,
                    &matrix,
                    qoff,
                    soff,
                    x,
                    &mut actual,
                );
                assert_eq!(actual, expected, "line {line}");
            }
            "approx" => {
                let qoff = number(&mut tokens);
                let match_end = number(&mut tokens);
                let soff = number(&mut tokens);
                let x = number(&mut tokens);
                let cutoff = number(&mut tokens);
                let expected = data(&mut tokens);
                let mut actual = BlastUngappedData {
                    q_start: -997,
                    s_start: -998,
                    length: -999,
                    score: -1000,
                };
                s_nucl_ungapped_extend(
                    &query,
                    &subject,
                    subject_length,
                    &matrix,
                    qoff,
                    match_end,
                    soff,
                    x,
                    &mut actual,
                    &score_table,
                    cutoff,
                );
                assert_eq!(actual, expected, "line {line}");
            }
            "lookup" => {
                let kind = number(&mut tokens);
                word = number(&mut tokens);
                width = number(&mut tokens);
                query = bytes(tokens.next().unwrap());
                subject = bytes(tokens.next().unwrap());
                subject_length = 257;
                let count = number(&mut tokens);
                let locations: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let masked = number(&mut tokens) != 0;
                let options = LookupTableOptions {
                    word_size: word,
                    mb_template_length: 0,
                    program_number: EBlastProgramType::Blastn,
                    db_filter: false,
                    stride: 0,
                    max_db_word_count: 0,
                };
                let query_options = QuerySetUpOptions {
                    filtering_options: Some(SBlastFilterOptions { mask_at_hash: true }),
                    filter_string: None,
                };
                let table = match kind {
                    1 => Table::Small(
                        blast_small_na_lookup_table_new(&query, &locations, &options, width)
                            .unwrap()
                            .0,
                    ),
                    2 => Table::Ordinary(blast_na_lookup_table_new(
                        &query,
                        &locations,
                        &options,
                        Some(&query_options),
                        width,
                    )),
                    0 => Table::Megablast(
                        blast_mb_lookup_table_new(
                            &query,
                            &locations,
                            &options,
                            Some(&query_options),
                            20000,
                            width,
                        )
                        .unwrap(),
                    ),
                    _ => unreachable!(),
                };
                assert_eq!(!table.masked().is_empty(), masked, "line {line}");
                lookup = Some(table);
            }
            "initial_new" => {
                let hash = number(&mut tokens) != 0;
                let window = number(&mut tokens);
                let scan = number(&mut tokens);
                let ungapped = number(&mut tokens) != 0;
                let cutoff = number(&mut tokens);
                let reduced = number(&mut tokens);
                let matrix_only = number(&mut tokens) != 0;
                extension = blast_extend_word_new(257, hash, window);
                hitlist = blast_init_hit_list_new();
                params = Some(BlastInitialWordParameters {
                    options: BlastInitialWordOptions {
                        gap_trigger: 0.0,
                        window_size: window,
                        scan_range: scan,
                        x_dropoff: 0.0,
                        program_number: EBlastProgramType::Blastn,
                    },
                    cutoffs: (0..2)
                        .map(|c| BlastUngappedCutoffs {
                            x_dropoff_init: if c == 0 { 20 } else { 7 },
                            x_dropoff: if c == 0 { 20 } else { 7 },
                            cutoff_score: cutoff + c,
                            reduced_nucl_cutoff_score: reduced + c,
                        })
                        .collect(),
                    nucl_score_table: score_table,
                    matrix_only_scoring: matrix_only,
                    ungapped_extension: ungapped,
                    x_dropoff_max: 0,
                    cutoff_score_min: 0,
                    container_type: if hash {
                        ESeedContainerType::DiagHash
                    } else {
                        ESeedContainerType::DiagArray
                    },
                });
            }
            "finder_stats_case" => {
                let present = number(&mut tokens) != 0;
                let total = number(&mut tokens);
                let extended = number(&mut tokens);
                let saved = number(&mut tokens);
                let mut stats = BlastUngappedStats {
                    lookup_hits: i64::MAX - 1,
                    num_seqs_lookup_hits: i32::MAX,
                    init_extends: i32::MAX - 1,
                    good_init_extends: i32::MAX,
                    num_seqs_passed: i32::MAX,
                };
                blast_ungapped_stats_update(
                    if present { Some(&mut stats) } else { None },
                    total,
                    extended,
                    saved,
                );
                assert_eq!(
                    stats.lookup_hits,
                    tokens.next().unwrap().parse::<i64>().unwrap()
                );
                assert_eq!(
                    [
                        stats.num_seqs_lookup_hits,
                        stats.init_extends,
                        stats.good_init_extends,
                        stats.num_seqs_passed
                    ],
                    [
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens)
                    ]
                );
            }
            "finder_ranges" => {
                let word = number(&mut tokens);
                let width = number(&mut tokens);
                let mut range = [0, number(&mut tokens), number(&mut tokens)];
                let ranges = [
                    SSeqRange { left: 0, right: 5 },
                    SSeqRange { left: 7, right: 9 },
                    SSeqRange {
                        left: 30,
                        right: 60,
                    },
                ];
                let result = s_determine_scanning_offsets(&ranges, word, width, &mut range);
                assert_eq!(result as i32, number(&mut tokens), "line {line}");
                assert_eq!(
                    range,
                    [
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens)
                    ],
                    "line {line}"
                );
            }
            "finder_stats_new" => finder_stats = BlastUngappedStats::default(),
            "finder_rollover" => {
                if let Some(d) = &mut extension.diag_table {
                    d.offset = i32::MAX / 4;
                } else {
                    extension.hash_table.as_mut().unwrap().offset = i32::MAX / 4;
                }
            }
            "finder_offset" => {
                let offset = extension.diag_table.as_ref().map_or_else(
                    || extension.hash_table.as_ref().unwrap().offset,
                    |d| d.offset,
                );
                assert_eq!(offset, number(&mut tokens), "line {line}");
            }
            "finder_subject" => {
                let length = number(&mut tokens);
                let packed = bytes(tokens.next().unwrap());
                let capacity = number(&mut tokens) as usize;
                let with_stats = number(&mut tokens) != 0;
                let mut pairs = vec![BlastOffsetPair::default(); capacity];
                let mut wrap = lookup.as_ref().unwrap().wrap();
                let extend = blast_choose_na_extend(&mut wrap);
                let scansub = match wrap.lut {
                    NucleotideLookup::Small(t) => {
                        NucleotideScanCallback::Small(if t.scan_step == 4 {
                            s_blast_small_na_scan_subject_8_4
                        } else {
                            assert_eq!(t.scan_step, 21);
                            s_blast_small_na_scan_subject_8_1_mod4
                        })
                    }
                    NucleotideLookup::Ordinary(t) => {
                        NucleotideScanCallback::Ordinary(if t.scan_step == 4 {
                            s_blast_na_scan_subject_8_4
                        } else {
                            s_blast_na_scan_subject_any
                        })
                    }
                    NucleotideLookup::Megablast(t) => {
                        NucleotideScanCallback::Megablast(match (t.lut_word_length, t.scan_step) {
                            (10, 2) => s_mb_scan_subject_10_2,
                            (11, 1) => s_mb_scan_subject_11_1_mod4,
                            (11, 18) => s_mb_scan_subject_11_2_mod4,
                            (12, 17) => s_mb_scan_subject_any,
                            _ => panic!("unsupported fixture scanner"),
                        })
                    }
                };
                let mut sentinels = vec![15];
                sentinels.extend_from_slice(&query);
                sentinels.push(15);
                let compressed =
                    btyper_blast::core::blast_util::blast_compress_blastna_sequence(&query);
                let q = NucleotideExtensionQuery {
                    sequence_start: &sentinels,
                    compressed_nuc_seq_start: Some(&compressed),
                };
                blast_init_hit_list_reset(&mut hitlist);
                let result = blast_na_word_finder(
                    &PackedNucleotideSubject {
                        sequence: &packed,
                        mask_type: ESubjectMaskingType::NoSubjMasking,
                    },
                    length,
                    &[SSeqRange {
                        left: 0,
                        right: length,
                    }],
                    &q,
                    &info,
                    &wrap,
                    &matrix,
                    params.as_ref().unwrap(),
                    &mut extension,
                    &mut pairs,
                    capacity as i32,
                    &mut hitlist,
                    if with_stats {
                        Some(&mut finder_stats)
                    } else {
                        None
                    },
                    scansub,
                    extend,
                );
                assert_eq!(result as i32, number(&mut tokens), "line {line}");
                assert_eq!(hitlist.total, number(&mut tokens), "line {line}");
                let mut expected = Vec::new();
                let mut actual = Vec::new();
                for h in &hitlist.init_hsp_array {
                    let q = number(&mut tokens);
                    let s = number(&mut tokens);
                    let present = number(&mut tokens) != 0;
                    let d = if present {
                        let d = data(&mut tokens);
                        [d.q_start, d.s_start, d.length, d.score]
                    } else {
                        [0; 4]
                    };
                    expected.push((q, s, present, d));
                    actual.push((
                        h.offsets.q_off,
                        h.offsets.s_off,
                        h.ungapped_data.is_some(),
                        h.ungapped_data
                            .map_or([0; 4], |d| [d.q_start, d.s_start, d.length, d.score]),
                    ));
                }
                if params.as_ref().unwrap().ungapped_extension {
                    for (a, b) in actual.iter().zip(&expected) {
                        assert_eq!((a.2, a.3), (b.2, b.3), "sorted keys line {line}");
                    }
                    actual.sort();
                    expected.sort();
                }
                assert_eq!(actual, expected, "line {line}");
                assert_eq!(
                    finder_stats.lookup_hits,
                    tokens.next().unwrap().parse::<i64>().unwrap()
                );
                assert_eq!(
                    [
                        finder_stats.num_seqs_lookup_hits,
                        finder_stats.init_extends,
                        finder_stats.good_init_extends,
                        finder_stats.num_seqs_passed
                    ],
                    [
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens)
                    ],
                    "line {line}"
                );
            }
            "mini_tables" => {
                for i in 0..256 {
                    assert_eq!(
                        [
                            S_EXACT_MATCH_EXTEND_LEFT[i] as i32,
                            S_EXACT_MATCH_EXTEND_RIGHT[i] as i32
                        ],
                        [number(&mut tokens), number(&mut tokens)]
                    );
                }
            }
            "mini_seeds" => {
                let count = number(&mut tokens);
                seeds = (0..count)
                    .map(|_| btyper_blast::core::blast_nascan::BlastOffsetPair {
                        q_off: number(&mut tokens),
                        s_off: number(&mut tokens),
                    })
                    .collect();
            }
            "mini_mode" => {
                mini_id = number(&mut tokens) as usize;
                let selected = number(&mut tokens) as usize;
                let mut wrap = lookup.as_ref().unwrap().wrap();
                let callback = blast_choose_na_extend(&mut wrap);
                assert!(
                    std::ptr::fn_addr_eq(callback, mini_callbacks[selected]),
                    "line {line}"
                );
                // The returned callback and stored membership callback must agree.
                let expected_member = match wrap.lut {
                    NucleotideLookup::Small(_) => s_small_na_lookup as TLookupCallback,
                    NucleotideLookup::Ordinary(_) => s_na_lookup as TLookupCallback,
                    NucleotideLookup::Megablast(_) => s_mb_lookup as TLookupCallback,
                };
                assert!(std::ptr::fn_addr_eq(wrap.lookup_callback, expected_member));
            }
            "mini_batch" => {
                let start = number(&mut tokens) as usize;
                let count = number(&mut tokens) as usize;
                let previous = hitlist.total as usize;
                let mut sentinels = Vec::with_capacity(query.len() + 2);
                sentinels.push(15);
                sentinels.extend_from_slice(&query);
                sentinels.push(15);
                let compressed =
                    btyper_blast::core::blast_util::blast_compress_blastna_sequence(&query);
                let q = NucleotideExtensionQuery {
                    sequence_start: &sentinels,
                    compressed_nuc_seq_start: Some(&compressed),
                };
                let ready = mini_callbacks[mini_id](
                    &seeds[start..start + count],
                    params.as_ref().unwrap(),
                    &lookup.as_ref().unwrap().wrap(),
                    &q,
                    &subject,
                    subject_length,
                    &matrix,
                    &info,
                    &mut extension,
                    &mut hitlist,
                    257,
                );
                assert_eq!(ready, number(&mut tokens), "line {line}");
                assert_eq!(hitlist.total, number(&mut tokens), "line {line}");
                for h in &hitlist.init_hsp_array[previous..] {
                    assert_eq!(
                        [
                            h.offsets.q_off,
                            h.offsets.s_off,
                            h.ungapped_data.is_some() as i32
                        ],
                        [
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens)
                        ],
                        "line {line}"
                    );
                    if let Some(d) = h.ungapped_data {
                        assert_eq!(d, data(&mut tokens), "line {line}");
                    }
                }
            }
            "initial_prime" => {
                let hash = number(&mut tokens) != 0;
                let side = number(&mut tokens);
                let level = number(&mut tokens);
                extension =
                    blast_extend_word_new(257, hash, params.as_ref().unwrap().options.window_size);
                blast_init_hit_list_reset(&mut hitlist);
                if hash {
                    s_blast_diag_hash_insert(
                        extension.hash_table.as_mut().unwrap(),
                        side,
                        level,
                        word,
                        0,
                        0,
                        100000,
                    );
                } else {
                    let d = extension.diag_table.as_mut().unwrap();
                    let slot = ((d.diag_array_length + side) & d.diag_mask) as usize;
                    d.hit_level_array[slot].last_hit_and_flag = level as u32;
                    d.hit_len_array.as_mut().unwrap()[slot] = word as u8;
                }
            }
            "initial_restart" => {
                let hash = extension.hash_table.is_some();
                extension =
                    blast_extend_word_new(257, hash, params.as_ref().unwrap().options.window_size);
                blast_init_hit_list_reset(&mut hitlist);
            }
            "initial_hit" => {
                let q = number(&mut tokens);
                let s = number(&mut tokens);
                let masks = number(&mut tokens) != 0;
                let range = number(&mut tokens);
                let table = lookup.as_ref().unwrap();
                let ready = if let Some(d) = &mut extension.diag_table {
                    s_blastn_diag_table_extend_initial_hit(
                        &query,
                        &subject,
                        subject_length,
                        q,
                        s,
                        table.masked(),
                        &info,
                        range,
                        word,
                        width,
                        &table.wrap(),
                        params.as_ref().unwrap(),
                        &matrix,
                        d,
                        &mut hitlist,
                        masks,
                    )
                } else {
                    s_blastn_diag_hash_extend_initial_hit(
                        &query,
                        &subject,
                        subject_length,
                        q,
                        s,
                        table.masked(),
                        &info,
                        range,
                        word,
                        width,
                        &table.wrap(),
                        params.as_ref().unwrap(),
                        &matrix,
                        extension.hash_table.as_mut().unwrap(),
                        &mut hitlist,
                        masks,
                    )
                };
                assert_eq!(ready, number(&mut tokens), "line {line}");
                assert_eq!(hitlist.total, number(&mut tokens), "line {line}");
                if let Some(h) = hitlist.init_hsp_array.last() {
                    assert_eq!(
                        [
                            h.offsets.q_off,
                            h.offsets.s_off,
                            h.ungapped_data.is_some() as i32
                        ],
                        [
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens)
                        ],
                        "line {line}"
                    );
                    if let Some(d) = h.ungapped_data {
                        assert_eq!(d, data(&mut tokens), "line {line}");
                    }
                }
            }
            "initial_array" => {
                let d = extension.diag_table.as_ref().unwrap();
                assert_eq!(d.diag_array_length, number(&mut tokens));
                for (i, c) in d.hit_level_array.iter().enumerate() {
                    let last = (c.last_hit_and_flag.wrapping_shl(1) as i32) >> 1;
                    let flag = (c.last_hit_and_flag >> 31) as i32;
                    let len = d.hit_len_array.as_ref().map_or(-1, |v| v[i] as i32);
                    assert_eq!(
                        [last, flag, len],
                        [
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens)
                        ],
                        "line {line} cell {i}"
                    );
                }
            }
            "initial_hash" => {
                let h = extension.hash_table.as_ref().unwrap();
                assert_eq!(
                    [h.capacity, h.occupancy],
                    [number(&mut tokens) as u32, number(&mut tokens) as u32]
                );
                for v in &h.backbone {
                    assert_eq!(*v, number(&mut tokens) as u32, "line {line}");
                }
                for c in &h.chain[..h.occupancy as usize] {
                    assert_eq!(
                        [
                            c.diag,
                            (c.level_and_saved.wrapping_shl(1) as i32) >> 1,
                            c.hit_len,
                            (c.level_and_saved >> 31) as i32,
                            c.next as i32
                        ],
                        [
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens),
                            number(&mut tokens)
                        ],
                        "line {line}"
                    );
                }
            }
            "member" => {
                let index = number(&mut tokens);
                let pos = number(&mut tokens);
                let expected = number(&mut tokens) != 0;
                let wrap = lookup.as_ref().unwrap().wrap();
                assert_eq!(
                    (wrap.lookup_callback)(&wrap, index, pos),
                    expected,
                    "line {line}"
                );
            }
            "masked" => {
                let soff = number(&mut tokens);
                let pos = number(&mut tokens);
                let expected = number(&mut tokens) != 0;
                assert_eq!(
                    s_is_seed_masked(&lookup.as_ref().unwrap().wrap(), &subject, soff, width, pos),
                    expected,
                    "line {line}"
                );
            }
            "word" => {
                let mut qoff = number(&mut tokens);
                let mut soff = number(&mut tokens);
                let has_locations = number(&mut tokens) != 0;
                let double = number(&mut tokens) != 0;
                let expected_type = number(&mut tokens);
                let expected_q = number(&mut tokens);
                let expected_s = number(&mut tokens);
                let expected_ext = number(&mut tokens);
                let expected_mask = number(&mut tokens) != 0;
                let table = lookup.as_ref().unwrap();
                assert_eq!(!table.masked().is_empty(), expected_mask);
                let mut extended = -999;
                let actual = s_type_of_word(
                    &subject,
                    &mut qoff,
                    &mut soff,
                    if has_locations { table.masked() } else { &[] },
                    &info,
                    subject_length as u32,
                    word as u32,
                    width as u32,
                    &table.wrap(),
                    double,
                    &mut extended,
                );
                assert_eq!(
                    (actual, qoff, soff, extended),
                    (expected_type, expected_q, expected_s, expected_ext),
                    "line {line}"
                );
            }
            "hash_new" => {
                let num_buckets = number(&mut tokens) as u32;
                let occupancy = number(&mut tokens) as u32;
                let capacity = number(&mut tokens) as u32;
                hash = Some(BlastDiagHash {
                    num_buckets,
                    occupancy,
                    capacity,
                    backbone: vec![0; num_buckets as usize],
                    chain: vec![DiagHashCell::default(); capacity as usize],
                    ..Default::default()
                });
            }
            "hash_put" => {
                let diag = number(&mut tokens);
                let level = number(&mut tokens);
                let len = number(&mut tokens);
                let flag = number(&mut tokens);
                let off = number(&mut tokens);
                let window = number(&mut tokens);
                let result = s_blast_diag_hash_insert(
                    hash.as_mut().unwrap(),
                    diag,
                    level,
                    len,
                    flag,
                    off,
                    window,
                );
                assert_eq!(result, number(&mut tokens), "line {line}");
                let h = hash.as_ref().unwrap();
                assert_eq!(h.capacity, number(&mut tokens) as u32, "line {line}");
                assert_eq!(h.occupancy, number(&mut tokens) as u32, "line {line}");
                assert_eq!(h.chain.len(), h.capacity as usize);
            }
            "hash_get" => {
                let diag = number(&mut tokens);
                let expected = number(&mut tokens);
                let mut level = -707;
                let mut len = -708;
                let mut flag = -709;
                assert_eq!(
                    s_blast_diag_hash_retrieve(
                        hash.as_ref().unwrap(),
                        diag,
                        &mut level,
                        &mut len,
                        &mut flag
                    ),
                    expected,
                    "line {line}"
                );
                assert_eq!(
                    (level, len, flag),
                    (
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens)
                    ),
                    "line {line}"
                );
            }
            "hash_snapshot" => {
                let h = hash.as_ref().unwrap();
                assert_eq!(h.capacity, number(&mut tokens) as u32);
                assert_eq!(h.occupancy, number(&mut tokens) as u32);
                let used = number(&mut tokens);
                let mut backbone = vec![0; h.num_buckets as usize];
                for _ in 0..used {
                    let i = number(&mut tokens) as usize;
                    backbone[i] = number(&mut tokens) as u32;
                }
                assert_eq!(h.backbone, backbone, "line {line}");
                for i in 0..h.occupancy as usize {
                    let diag = number(&mut tokens);
                    let level = number(&mut tokens);
                    let len = number(&mut tokens);
                    let flag = number(&mut tokens);
                    let next = number(&mut tokens);
                    let expected = DiagHashCell {
                        diag,
                        level_and_saved: ((level as u32) & 0x7fffffff) | ((flag as u32 & 1) << 31),
                        hit_len: len,
                        next: next as u32,
                    };
                    assert_eq!(h.chain[i], expected, "line {line}");
                }
            }
            _ => panic!("Unknown record on line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}
#[test]
fn original_nucleotide_ungapped_kernels() {
    verify(include_str!("fixtures/ungapped-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_UNGAPPED_FIXTURE"]
fn freshly_compiled_nucleotide_ungapped_kernels() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_UNGAPPED_FIXTURE").unwrap()).unwrap());
}

#[test]
fn original_diagonal_initial_hit_extension() {
    verify(include_str!("fixtures/initial-hit-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_INITIAL_HIT_FIXTURE"]
fn freshly_compiled_diagonal_initial_hit_extension() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_INITIAL_HIT_FIXTURE").unwrap()).unwrap());
}

#[test]
fn original_nucleotide_mini_extension() {
    verify(include_str!("fixtures/mini-extension-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_MINI_EXTENSION_FIXTURE"]
fn freshly_compiled_nucleotide_mini_extension() {
    verify(
        &std::fs::read_to_string(std::env::var("BLAST_MINI_EXTENSION_FIXTURE").unwrap()).unwrap(),
    );
}

#[test]
fn original_nucleotide_word_finder() {
    verify(include_str!("fixtures/word-finder-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_WORD_FINDER_FIXTURE"]
fn freshly_compiled_nucleotide_word_finder() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_WORD_FINDER_FIXTURE").unwrap()).unwrap());
}
