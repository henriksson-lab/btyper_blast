use btyper_blast::core::blast_lookup::SSeqRange;
use btyper_blast::core::blast_nalookup::*;
use btyper_blast::core::blast_nascan::*;
use btyper_blast::core::blast_options::*;
use btyper_blast::core::blast_program::EBlastProgramType;
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

fn query_options(mode: i32) -> Option<QuerySetUpOptions> {
    if mode == 0 {
        return None;
    }
    let filtering_options = matches!(mode, 2 | 3 | 9 | 10).then_some(SBlastFilterOptions {
        mask_at_hash: mode == 3 || mode == 9,
    });
    let filter_string: Option<&[u8]> = match mode {
        4 | 10 => Some(b"m"),
        5 | 9 => Some(b"F"),
        6 => Some(b"F; m "),
        7 => Some(b"am\0ignored"),
        8 => Some(b"F\0m"),
        _ => None,
    };
    Some(QuerySetUpOptions {
        filtering_options,
        filter_string: filter_string.map(Vec::from),
    })
}

fn verify(fixture: &str) {
    let mut current = None;
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "mask_at_hash" => {
                let options = query_options(number(&mut tokens));
                assert_eq!(
                    s_has_mask_at_hash_enabled(options.as_ref()),
                    number(&mut tokens) != 0
                );
            }
            "na" => {
                let word_size = number(&mut tokens);
                let width = number(&mut tokens);
                let query_opts = query_options(number(&mut tokens));
                let query = bytes(tokens.next().unwrap());
                let count = number(&mut tokens);
                let locations: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let options = LookupTableOptions {
                    stride: 0,
                    max_db_word_count: 0,
                    word_size,
                    program_number: EBlastProgramType::Blastn,
                    db_filter: false,
                    mb_template_length: 0,
                };
                if query.len() > 32768 {
                    let mut max_offset = 0;
                    let estimate = estimate_num_table_entries(&locations, &mut max_offset);
                    let mut selected_width = 0;
                    assert_eq!(
                        blast_choose_na_lookup_table(
                            &options,
                            estimate,
                            max_offset,
                            &mut selected_width
                        ),
                        ELookupTableType::NaLookupTable
                    );
                    assert_eq!(selected_width, width);
                }
                let mut lookup = blast_na_lookup_table_new(
                    &query,
                    &locations,
                    &options,
                    query_opts.as_ref(),
                    width,
                );
                assert_eq!(lookup.word_length, word_size);
                assert_eq!(lookup.lut_word_length, width);
                assert_eq!(lookup.scan_step, word_size - width + 1);
                assert_eq!(lookup.backbone_size, 1 << (2 * width));
                assert_eq!(lookup.mask, lookup.backbone_size - 1);
                assert_eq!(lookup.longest_chain, number(&mut tokens), "line {line}");
                assert_eq!(lookup.overflow_size, number(&mut tokens), "line {line}");
                let occupied = number(&mut tokens);
                assert_eq!(
                    lookup
                        .thick_backbone
                        .iter()
                        .filter(|cell| cell.num_used > 0)
                        .count(),
                    occupied as usize
                );
                for _ in 0..occupied {
                    let index = number(&mut tokens);
                    let expected_hits = number(&mut tokens);
                    let expected_payload = [
                        number(&mut tokens),
                        number(&mut tokens),
                        number(&mut tokens),
                    ];
                    let cell = &lookup.thick_backbone[index as usize];
                    assert_eq!(cell.num_used, expected_hits, "line {line}");
                    assert_eq!(cell.payload, expected_payload, "line {line}");
                    assert_eq!(s_blast_lookup_get_num_hits(&lookup, index), expected_hits);
                }
                for value in &lookup.overflow {
                    assert_eq!(*value, number(&mut tokens), "line {line}");
                }
                let pv_len = number(&mut tokens) as usize;
                let occupied = number(&mut tokens);
                let mut expected_pv = vec![0; pv_len];
                for _ in 0..occupied {
                    let index = number(&mut tokens) as usize;
                    expected_pv[index] = tokens.next().unwrap().parse::<u32>().unwrap();
                }
                assert_eq!(lookup.pv, expected_pv, "line {line}");
                let count = number(&mut tokens);
                let masks: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                assert_eq!(lookup.masked_locations, masks, "line {line}");
                let aligned = number(&mut tokens) != 0;
                assert_eq!(aligned, width == 8 && lookup.scan_step == 4);
                s_na_choose_scan_subject(&mut lookup);
                assert!(lookup.scansub_callback.is_some());
                current = Some(lookup);
            }
            "na_scan" => {
                let kind = number(&mut tokens);
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
                let lookup = current.as_ref().unwrap();
                let subject = PackedNucleotideSubject {
                    sequence: &sequence,
                    mask_type,
                };
                let scanner = if kind == 0 {
                    lookup.scansub_callback.unwrap()
                } else {
                    s_blast_na_scan_subject_any
                };
                let mut pairs = vec![BlastOffsetPair::default(); capacity as usize];
                let hits = scanner(lookup, &subject, &mut pairs, capacity, &mut range);
                assert_eq!(range[0], expected_stop, "line {line}");
                assert_eq!(hits, expected_hits, "line {line}");
                assert_eq!(&pairs[..hits as usize], expected, "line {line}");
            }
            _ => panic!("Unknown row on line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}

#[test]
fn original_c_ordinary_lookup() {
    verify(include_str!("fixtures/ordinary-lookup-2.17.0.txt"));
}

#[test]
#[ignore = "requires original C output via BLAST_ORDINARY_LOOKUP_FIXTURE"]
fn freshly_compiled_c_ordinary_lookup() {
    verify(
        &std::fs::read_to_string(std::env::var("BLAST_ORDINARY_LOOKUP_FIXTURE").unwrap()).unwrap(),
    );
}
