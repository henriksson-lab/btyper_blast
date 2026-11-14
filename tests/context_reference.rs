use btyper_blast::core::blast_encoding::EBlastEncoding;
use btyper_blast::core::blast_filter::*;
use btyper_blast::core::blast_lookup::SSeqRange;
use btyper_blast::core::blast_program::EBlastProgramType;
use btyper_blast::core::blast_query_info::*;
use btyper_blast::core::blast_setup::*;
use btyper_blast::{
    NucleotideTask, PreparedNucleotideBatch, PreparedNucleotideQuery, QueryBatchError,
};

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
fn program(value: i32) -> EBlastProgramType {
    match value {
        12 => EBlastProgramType::Blastn,
        41 => EBlastProgramType::Tblastn,
        22 => EBlastProgramType::Blastx,
        60 => EBlastProgramType::Tblastx,
        524 => EBlastProgramType::Mapping,
        _ => unreachable!(),
    }
}
fn encoding(value: i32) -> EBlastEncoding {
    match value {
        0 => EBlastEncoding::Protein,
        1 => EBlastEncoding::Nucleotide,
        2 => EBlastEncoding::Ncbi4na,
        3 => EBlastEncoding::Ncbi2na,
        _ => unreachable!(),
    }
}
fn strand(value: i32) -> ENaStrand {
    match value {
        0 => ENaStrand::Unknown,
        1 => ENaStrand::Plus,
        2 => ENaStrand::Minus,
        3 => ENaStrand::Both,
        4 => ENaStrand::BothRev,
        255 => ENaStrand::Other,
        _ => unreachable!(),
    }
}
fn sentinel(value: i32) -> ESentinelType {
    match value {
        0 => ESentinelType::Sentinels,
        1 => ESentinelType::NoSentinels,
        _ => unreachable!(),
    }
}
fn ranges<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> Vec<SSeqRange> {
    let count = number(tokens);
    (0..count)
        .map(|_| SSeqRange {
            left: number(tokens),
            right: number(tokens),
        })
        .collect()
}
fn info<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> BlastQueryInfo {
    let first_context = number(tokens);
    let last_context = number(tokens);
    let num_queries = number(tokens);
    let max_length = number(tokens) as u32;
    let min_length = number(tokens) as u32;
    let contexts = (0..=last_context)
        .map(|_| BlastContextInfo {
            query_offset: number(tokens),
            query_length: number(tokens),
            eff_searchsp: tokens.next().unwrap().parse().unwrap(),
            length_adjustment: number(tokens),
            query_index: number(tokens),
            frame: number(tokens) as i8,
            is_valid: number(tokens) != 0,
            segment_flags: number(tokens),
        })
        .collect();
    BlastQueryInfo {
        first_context,
        last_context,
        num_queries,
        max_length,
        min_length,
        contexts,
    }
}
fn verify(fixture: &str) {
    let mut current = None;
    let mut current_program = EBlastProgramType::Undefined;
    let mut batches = Vec::new();
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "new_none" => {
                let p = program(number(&mut tokens));
                let count = number(&mut tokens);
                assert_eq!(
                    blast_query_info_new(p, count).is_none(),
                    number(&mut tokens) != 0
                );
            }
            "new" => {
                let p = program(number(&mut tokens));
                let count = number(&mut tokens);
                assert_eq!(
                    blast_query_info_new(p, count).unwrap(),
                    info(&mut tokens),
                    "line {line}"
                );
            }
            "info" => {
                current_program = program(number(&mut tokens));
                let count = number(&mut tokens);
                let mut actual = blast_query_info_new(current_program, count).unwrap();
                let contexts = number(&mut tokens);
                for index in 0..contexts {
                    s_query_info_set_context(&mut actual, index as u32, number(&mut tokens) as u32);
                }
                let reassign = number(&mut tokens);
                if reassign >= 0 {
                    s_query_info_set_context(&mut actual, reassign as u32, 0);
                    s_query_info_set_context(&mut actual, reassign as u32, 3);
                }
                actual.first_context = number(&mut tokens);
                actual.max_length = actual
                    .contexts
                    .iter()
                    .map(|ctx| ctx.query_length as u32)
                    .max()
                    .unwrap();
                actual.min_length = actual
                    .contexts
                    .iter()
                    .map(|ctx| ctx.query_length as u32)
                    .min()
                    .unwrap();
                assert_eq!(actual, info(&mut tokens), "line {line}");
                assert_eq!(
                    query_info_get_seq_buf_len(&actual),
                    number(&mut tokens) as u32,
                    "line {line}"
                );
                current = Some(actual);
            }
            "offset" => {
                let offset = number(&mut tokens);
                let actual = current.as_ref().unwrap();
                assert_eq!(
                    b_search_context_info(offset, actual),
                    number(&mut tokens),
                    "line {line}"
                );
                assert_eq!(
                    blast_get_query_index_from_query_offset(offset, current_program, actual),
                    number(&mut tokens),
                    "line {line}"
                );
            }
            "complement" => {
                let p = program(number(&mut tokens));
                let present = number(&mut tokens) != 0;
                let output = number(&mut tokens) != 0;
                let qinfo = current.as_ref().unwrap();
                let mut masks = BlastMaskLoc {
                    total_size: qinfo.last_context + 1,
                    seqloc_array: (0..=qinfo.last_context)
                        .map(|_| ranges(&mut tokens))
                        .collect(),
                };
                let expected_status = number(&mut tokens);
                let expected = ranges(&mut tokens);
                let expected_masks: Vec<_> = (0..=qinfo.last_context)
                    .map(|_| ranges(&mut tokens))
                    .collect();
                let initial = vec![SSeqRange {
                    left: -100,
                    right: -100,
                }];
                let mut actual = initial.clone();
                let status = blast_complement_mask_locations(
                    p,
                    qinfo,
                    if present { Some(&mut masks) } else { None },
                    if output { Some(&mut actual) } else { None },
                );
                assert_eq!(status, expected_status as i16, "line {line}");
                assert_eq!(
                    actual,
                    if output { expected } else { initial },
                    "line {line}"
                );
                assert_eq!(masks.seqloc_array, expected_masks, "line {line}");
            }
            "buffer" => {
                let length = tokens.next().unwrap().parse().unwrap();
                let enc = encoding(number(&mut tokens));
                let dir = strand(number(&mut tokens));
                let sent = sentinel(number(&mut tokens));
                let expected = tokens.next().unwrap().parse::<u32>().unwrap();
                assert_eq!(
                    calculate_seq_buffer_length(length, enc, dir, sent),
                    expected,
                    "line {line}"
                );
            }
            "sentinel" => {
                let enc = encoding(number(&mut tokens));
                assert_eq!(get_sentinel_byte(enc), number(&mut tokens) as u8);
            }
            "strand" => {
                let enc = encoding(number(&mut tokens));
                let sent = sentinel(number(&mut tokens));
                let dir = strand(number(&mut tokens));
                let plus = bytes(tokens.next().unwrap());
                let minus = bytes(tokens.next().unwrap());
                let expected = bytes(tokens.next().unwrap());
                let actual = if dir == ENaStrand::Both {
                    get_sequence_nucleotide_both_strands(&plus, &minus, enc, sent)
                } else {
                    get_sequence_single_nucleotide_strand(&plus, &minus, enc, dir, sent)
                };
                assert_eq!(actual, expected, "line {line}");
            }
            "batch" => {
                let count = number(&mut tokens);
                let input: Vec<_> = (0..count).map(|_| bytes(tokens.next().unwrap())).collect();
                let expected_info = info(&mut tokens);
                let expected_buffer = bytes(tokens.next().unwrap());
                let expected_ranges = ranges(&mut tokens);
                batches.clear();
                for task in [NucleotideTask::Blastn, NucleotideTask::Megablast] {
                    let queries: Vec<_> = input
                        .iter()
                        .enumerate()
                        .map(|(i, input)| {
                            PreparedNucleotideQuery::new(format!("query_{i}"), input, task).unwrap()
                        })
                        .collect();
                    let actual = PreparedNucleotideBatch::new(&queries).unwrap();
                    assert_eq!(actual.query_info(), &expected_info, "line {line}");
                    assert_eq!(
                        actual.sequence_with_sentinels(),
                        expected_buffer,
                        "line {line}"
                    );
                    assert_eq!(
                        actual.blastna_contexts(),
                        &expected_buffer[1..expected_buffer.len() - 1],
                        "line {line}"
                    );
                    assert_eq!(actual.lookup_regions(), expected_ranges, "line {line}");
                    assert_eq!(
                        actual.identifiers(),
                        (0..count).map(|i| format!("query_{i}")).collect::<Vec<_>>()
                    );
                    assert_eq!(actual.task(), task);
                    batches.push(actual);
                }
            }
            "batch_lookup" => {
                let word = number(&mut tokens);
                let kind = number(&mut tokens);
                let width = number(&mut tokens);
                let batch = &batches[usize::from(word == 28)];
                let lookup = batch.seed_index().unwrap();
                assert_eq!(lookup.table_type() as i32, kind, "line {line}");
                assert_eq!(lookup.lookup_width(), width, "line {line}");
                assert_eq!(lookup.query_blastna(), batch.blastna_contexts());
            }
            _ => panic!("Unknown row on line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}

#[test]
fn original_query_context_setup() {
    verify(include_str!("fixtures/context-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C/C++ output via BLAST_CONTEXT_FIXTURE"]
fn freshly_compiled_query_context_setup() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_CONTEXT_FIXTURE").unwrap()).unwrap());
}
#[test]
fn native_batch_rejects_empty_and_mixed_tasks() {
    assert_eq!(
        PreparedNucleotideBatch::new(&[]).unwrap_err(),
        QueryBatchError::EmptyBatch
    );
    let plus = PreparedNucleotideQuery::new("same_id", b"ACGT", NucleotideTask::Blastn).unwrap();
    let other =
        PreparedNucleotideQuery::new("same_id", b"ACGT", NucleotideTask::Megablast).unwrap();
    assert_eq!(
        PreparedNucleotideBatch::new(&[plus, other]).unwrap_err(),
        QueryBatchError::MixedTasks
    );
}
