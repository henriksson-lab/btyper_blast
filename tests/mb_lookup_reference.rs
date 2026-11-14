use btyper_blast::core::blast_lookup::SSeqRange;
use btyper_blast::core::blast_nalookup::{
    blast_mb_lookup_table_new, BlastMBLookupTable, LookupTableOptions,
};
use btyper_blast::core::blast_nascan::*;
use btyper_blast::core::blast_options::{QuerySetUpOptions, SBlastFilterOptions};
use btyper_blast::core::blast_program::EBlastProgramType;
use btyper_blast::core::lookup_util::ilog2;

fn number<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> i32 {
    tokens.next().unwrap().parse().unwrap()
}

fn bytes(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

type Scanner = fn(
    &BlastMBLookupTable,
    &PackedNucleotideSubject<'_>,
    &mut [BlastOffsetPair],
    i32,
    &mut [i32; 2],
) -> i32;

fn verify(fixture: &str) {
    let mut current = None;
    let mut selected: Scanner = s_mb_scan_subject_any;
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "mb" => {
                let approx = number(&mut tokens);
                let mask_at_hash = number(&mut tokens) != 0;
                let query = bytes(tokens.next().unwrap());
                let count = number(&mut tokens);
                let locations: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let word_length = number(&mut tokens);
                let lut_word_length = number(&mut tokens);
                let scan_step = number(&mut tokens);
                let length = number(&mut tokens);
                let longest_chain = number(&mut tokens);
                let pv_array_bts = number(&mut tokens);
                selected = match number(&mut tokens) {
                    0 => s_mb_scan_subject_any,
                    1 => s_mb_scan_subject_10_2,
                    2 => s_mb_scan_subject_11_1_mod4,
                    3 => s_mb_scan_subject_11_2_mod4,
                    _ => unreachable!(),
                };
                let hashsize = 1_i64 << (2 * lut_word_length);
                let mut hashtable = vec![0; hashsize as usize];
                let count = number(&mut tokens);
                for _ in 0..count {
                    let index = number(&mut tokens) as usize;
                    hashtable[index] = number(&mut tokens);
                }
                let next_pos = (0..=length).map(|_| number(&mut tokens)).collect();
                let size = number(&mut tokens);
                let count = number(&mut tokens);
                let mut pv_array = vec![0; size as usize];
                for _ in 0..count {
                    let index = number(&mut tokens) as usize;
                    pv_array[index] = tokens.next().unwrap().parse().unwrap();
                }
                let expected = BlastMBLookupTable {
                    word_length,
                    lut_word_length,
                    scan_step,
                    hashsize,
                    longest_chain,
                    hashtable,
                    next_pos,
                    pv_array,
                    pv_array_bts,
                    ..Default::default()
                };
                let count = number(&mut tokens);
                let masks: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let options = LookupTableOptions {
                    word_size: word_length,
                    mb_template_length: 0,
                    program_number: EBlastProgramType::Blastn,
                    db_filter: false,
                    stride: 0,
                    max_db_word_count: 0,
                };
                let qopts = QuerySetUpOptions {
                    filtering_options: Some(SBlastFilterOptions { mask_at_hash }),
                    filter_string: None,
                };
                let lookup = blast_mb_lookup_table_new(
                    &query,
                    &locations,
                    &options,
                    Some(&qopts),
                    approx,
                    lut_word_length,
                )
                .unwrap();
                assert_eq!(lookup.word_length, expected.word_length);
                assert_eq!(lookup.lut_word_length, expected.lut_word_length);
                assert_eq!(lookup.scan_step, expected.scan_step);
                assert_eq!(lookup.hashsize, expected.hashsize);
                assert_eq!(lookup.longest_chain, expected.longest_chain, "line {line}");
                assert_eq!(lookup.pv_array_bts, expected.pv_array_bts);
                assert_eq!(lookup.hashtable, expected.hashtable, "line {line}");
                assert_eq!(lookup.next_pos, expected.next_pos, "line {line}");
                assert_eq!(lookup.pv_array, expected.pv_array, "line {line}");
                assert_eq!(lookup.masked_locations, masks, "line {line}");
                assert!(!lookup.stride);
                current = Some(lookup);
            }
            "ilog2" => {
                let value = tokens.next().unwrap().parse().unwrap();
                assert_eq!(ilog2(value), number(&mut tokens));
            }
            "mb_empty" => {
                let options = LookupTableOptions {
                    word_size: 28,
                    mb_template_length: 0,
                    program_number: EBlastProgramType::Blastn,
                    db_filter: false,
                    stride: 0,
                    max_db_word_count: 0,
                };
                assert_eq!(
                    blast_mb_lookup_table_new(&[], &[], &options, None, 0, 11).unwrap_err(),
                    number(&mut tokens) as i16
                );
            }
            "mb_probe" => {
                let index = i64::from(number(&mut tokens));
                let expected_pv = number(&mut tokens);
                let expected_hits = number(&mut tokens);
                let expected: Vec<_> = (0..expected_hits)
                    .map(|_| BlastOffsetPair {
                        q_off: number(&mut tokens),
                        s_off: number(&mut tokens),
                    })
                    .collect();
                let lookup = current.as_ref().unwrap();
                assert_eq!(
                    s_blast_mb_lookup_has_hits(lookup, index),
                    expected_pv,
                    "line {line}"
                );
                let mut pairs = vec![BlastOffsetPair::default(); lookup.next_pos.len()];
                let hits = s_blast_mb_lookup_retrieve(lookup, index, &mut pairs, 79);
                assert_eq!(hits, expected_hits, "line {line}");
                assert_eq!(&pairs[..hits as usize], expected, "line {line}");
            }
            "mb_scan" => {
                let generic = number(&mut tokens) != 0;
                let mask_type = match number(&mut tokens) {
                    0 => ESubjectMaskingType::NoSubjMasking,
                    1 => ESubjectMaskingType::SoftSubjMasking,
                    2 => ESubjectMaskingType::HardSubjMasking,
                    _ => unreachable!(),
                };
                let mut range = [number(&mut tokens), number(&mut tokens)];
                let capacity = number(&mut tokens);
                let sequence = bytes(tokens.next().unwrap());
                let expected_stop = number(&mut tokens);
                let expected_hits = number(&mut tokens);
                let expected: Vec<_> = (0..expected_hits)
                    .map(|_| BlastOffsetPair {
                        q_off: number(&mut tokens),
                        s_off: number(&mut tokens),
                    })
                    .collect();
                let mut pairs = vec![BlastOffsetPair::default(); capacity as usize];
                let subject = PackedNucleotideSubject {
                    sequence: &sequence,
                    mask_type,
                };
                let scanner = if generic {
                    s_mb_scan_subject_any
                } else {
                    selected
                };
                let hits = scanner(
                    current.as_ref().unwrap(),
                    &subject,
                    &mut pairs,
                    capacity,
                    &mut range,
                );
                assert_eq!(range[0], expected_stop, "line {line}");
                assert_eq!(hits, expected_hits, "line {line}");
                assert_eq!(&pairs[..hits as usize], expected, "line {line}");
            }
            _ => panic!("Unknown record at line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}

#[test]
fn original_c_mb_lookup() {
    verify(include_str!("fixtures/mb-lookup-2.17.0.txt"));
}

#[test]
#[ignore = "requires original C constructor/scanner output via BLAST_MB_LOOKUP_FIXTURE"]
fn freshly_compiled_c_mb_lookup() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_MB_LOOKUP_FIXTURE").unwrap()).unwrap());
}
