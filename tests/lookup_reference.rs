use btyper_blast::core::blast_lookup::*;
use btyper_blast::core::blast_nalookup::*;
use btyper_blast::core::blast_nascan::*;
use btyper_blast::core::blast_program::EBlastProgramType;
use btyper_blast::core::blast_util::blast_compress_blastna_sequence;
use btyper_blast::core::lookup_util::estimate_num_table_entries;

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

fn verify(fixture: &str) {
    let mut small = None;
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "choose" => {
                let word_size = number(&mut tokens);
                let program_number = match number(&mut tokens) {
                    12 => EBlastProgramType::Blastn,
                    524 => EBlastProgramType::Mapping,
                    _ => unreachable!(),
                };
                let options = LookupTableOptions {
                    stride: 0,
                    max_db_word_count: 0,
                    word_size,
                    program_number,
                    db_filter: number(&mut tokens) != 0,
                    mb_template_length: number(&mut tokens),
                };
                let entries = number(&mut tokens);
                let max_offset = number(&mut tokens);
                let expected_type = number(&mut tokens);
                let expected_width = number(&mut tokens);
                let mut width = 0;
                let kind = blast_choose_na_lookup_table(&options, entries, max_offset, &mut width);
                assert_eq!(
                    (kind as i32, width),
                    (expected_type, expected_width),
                    "line {line}"
                );
            }
            "index" => {
                let word_length = number(&mut tokens);
                let width = number(&mut tokens);
                let query = bytes(tokens.next().unwrap());
                let count = number(&mut tokens);
                let locations: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let expected_estimate = number(&mut tokens);
                let expected_max = number(&mut tokens);
                let mut max_offset = -1;
                assert_eq!(
                    estimate_num_table_entries(&locations, &mut max_offset),
                    expected_estimate
                );
                assert_eq!(max_offset, expected_max);
                let mut backbone = vec![None; 1 << (2 * width)];
                blast_lookup_index_query_exact_matches(
                    &mut backbone,
                    word_length,
                    2,
                    width,
                    &query,
                    &locations,
                );
                let occupied = number(&mut tokens) as usize;
                assert_eq!(
                    backbone.iter().filter(|chain| chain.is_some()).count(),
                    occupied,
                    "line {line}"
                );
                for _ in 0..occupied {
                    let index = number(&mut tokens) as usize;
                    let capacity = number(&mut tokens);
                    let hits = number(&mut tokens);
                    let expected: Vec<_> = (0..hits).map(|_| number(&mut tokens)).collect();
                    let chain = backbone[index].as_ref().unwrap();
                    assert_eq!(chain.len(), capacity as usize, "line {line}");
                    assert_eq!(&chain[..2], &[capacity, hits], "line {line}");
                    assert_eq!(&chain[2..2 + hits as usize], expected, "line {line}");
                }
            }
            "rolling" => {
                let query = bytes(tokens.next().unwrap());
                let expected = bytes(tokens.next().unwrap());
                assert_eq!(
                    blast_compress_blastna_sequence(&query),
                    expected,
                    "line {line}"
                );
            }
            "small" => {
                let word_size = number(&mut tokens);
                let width = number(&mut tokens);
                let query = bytes(tokens.next().unwrap());
                let count = number(&mut tokens);
                let locations: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let status = number(&mut tokens);
                let options = LookupTableOptions {
                    stride: 0,
                    max_db_word_count: 0,
                    word_size,
                    program_number: EBlastProgramType::Blastn,
                    db_filter: false,
                    mb_template_length: 0,
                };
                let result = blast_small_na_lookup_table_new(&query, &locations, &options, width);
                if status == 0 {
                    let (lookup, compressed) = result.unwrap();
                    assert_eq!(compressed, blast_compress_blastna_sequence(&query));
                    assert_eq!(lookup.scan_step, word_size - width + 1);
                    assert_eq!(lookup.mask, (1 << (2 * width)) - 1);
                    assert_eq!(lookup.backbone_size, 1 << (2 * width));
                    assert_eq!(lookup.longest_chain, number(&mut tokens), "line {line}");
                    assert_eq!(lookup.overflow_size, number(&mut tokens), "line {line}");
                    let occupied = number(&mut tokens);
                    assert_eq!(
                        lookup.final_backbone.iter().filter(|v| **v != -1).count(),
                        occupied as usize
                    );
                    for _ in 0..occupied {
                        let index = number(&mut tokens) as usize;
                        assert_eq!(
                            lookup.final_backbone[index],
                            number(&mut tokens) as i16,
                            "line {line}"
                        );
                    }
                    for value in &lookup.overflow[2..] {
                        assert_eq!(*value, number(&mut tokens) as i16, "line {line}");
                    }
                    let count = number(&mut tokens);
                    let masks: Vec<_> = (0..count)
                        .map(|_| SSeqRange {
                            left: number(&mut tokens),
                            right: number(&mut tokens),
                        })
                        .collect();
                    assert_eq!(lookup.masked_locations, masks, "line {line}");
                    small = Some(lookup);
                } else {
                    assert_eq!(result.unwrap_err(), status);
                    small = None;
                }
            }
            "scan" => {
                let mut range = [number(&mut tokens), number(&mut tokens)];
                let capacity = number(&mut tokens);
                let subject = bytes(tokens.next().unwrap());
                let expected_stop = number(&mut tokens);
                let expected_hits = number(&mut tokens);
                let expected: Vec<_> = (0..expected_hits)
                    .map(|_| BlastOffsetPair {
                        q_off: number(&mut tokens),
                        s_off: number(&mut tokens),
                    })
                    .collect();
                let lookup = small.as_ref().unwrap();
                let mut pairs = vec![BlastOffsetPair::default(); capacity as usize];
                let hits = if lookup.word_length == 11 {
                    s_blast_small_na_scan_subject_8_4(
                        lookup, &subject, &mut pairs, capacity, &mut range,
                    )
                } else {
                    s_blast_small_na_scan_subject_8_1_mod4(
                        lookup, &subject, &mut pairs, capacity, &mut range,
                    )
                };
                assert_eq!(range[0], expected_stop, "line {line}");
                assert_eq!(hits, expected_hits, "line {line}");
                assert_eq!(&pairs[..hits as usize], expected, "line {line}");
            }
            "finalize" => {
                let counts = [number(&mut tokens), number(&mut tokens)];
                let mut thin = vec![None; 4];
                for (k, count) in counts.into_iter().enumerate() {
                    if count > 0 {
                        let mut chain = vec![count + 2, count];
                        chain.extend(0..count);
                        thin[k] = Some(chain);
                    }
                }
                let mut lookup = BlastSmallNaLookupTable {
                    backbone_size: 4,
                    ..Default::default()
                };
                let mut compressed = None;
                let status = s_blast_small_na_lookup_finalize(
                    &mut thin,
                    &mut lookup,
                    &[0; 40],
                    &mut compressed,
                );
                assert_eq!(status, number(&mut tokens), "line {line}");
                assert!(thin.iter().all(|chain| chain.is_none()));
                assert_eq!(lookup.longest_chain, number(&mut tokens));
                assert_eq!(lookup.overflow_size, number(&mut tokens));
                for i in 0..4 {
                    let value = lookup.final_backbone.get(i).copied().unwrap_or(-99);
                    assert_eq!(value, number(&mut tokens) as i16);
                }
                if status == 0 {
                    for value in &lookup.overflow[2..] {
                        assert_eq!(*value, number(&mut tokens) as i16);
                    }
                }
                assert_eq!(compressed.is_some(), number(&mut tokens) != 0);
            }
            _ => panic!("Unknown record on line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}

#[test]
fn original_c_lookup_fixtures() {
    verify(include_str!("fixtures/lookup-2.17.0.txt"));
}

#[test]
#[ignore = "requires freshly compiled original C; see generate_lookup_reference.py"]
fn freshly_compiled_c_lookup_fixtures() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_LOOKUP_FIXTURE").unwrap()).unwrap());
}
