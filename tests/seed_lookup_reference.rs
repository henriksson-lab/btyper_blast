use btyper_blast::core::blast_lookup::SSeqRange;
use btyper_blast::core::blast_nalookup::ELookupTableType;
use btyper_blast::core::blast_nascan::BlastOffsetPair;
use btyper_blast::seed::{NucleotideSeedIndex, SeedError};
use btyper_blast::{DatabaseBuilder, NucleotideTask};

fn number<'a>(tokens: &mut impl Iterator<Item = &'a str>) -> i32 {
    tokens.next().unwrap().parse().unwrap()
}
fn bytes(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}
fn verify(fixture: &str) {
    let mut lookup = None;
    let mut database = None;
    let mut pairs = Vec::new();
    let mut next = 0;
    for (line, row) in fixture.lines().enumerate() {
        let mut tokens = row.split_whitespace();
        match tokens.next().unwrap() {
            "seed_query" => {
                let task = match number(&mut tokens) {
                    11 => NucleotideTask::Blastn,
                    28 => NucleotideTask::Megablast,
                    _ => unreachable!(),
                };
                let query = bytes(tokens.next().unwrap());
                let count = number(&mut tokens);
                let regions: Vec<_> = (0..count)
                    .map(|_| SSeqRange {
                        left: number(&mut tokens),
                        right: number(&mut tokens),
                    })
                    .collect();
                let kind = number(&mut tokens);
                let width = number(&mut tokens);
                let step = number(&mut tokens);
                let longest = number(&mut tokens);
                let index = NucleotideSeedIndex::from_segments(&query, &regions, task).unwrap();
                assert_eq!(index.query_blastna(), query);
                assert_eq!(index.table_type() as i32, kind, "line {line}");
                assert_eq!(index.lookup_width(), width, "line {line}");
                assert_eq!(index.scan_step(), step, "line {line}");
                assert_eq!(index.longest_chain(), longest, "line {line}");
                assert_eq!(index.recommended_capacity(), 4096 + longest as usize);
                assert_eq!(
                    index.compressed_query().is_some(),
                    index.table_type() == ELookupTableType::SmallNaLookupTable
                );
                lookup = Some(index);
            }
            "seed_subject" => {
                let raw = bytes(tokens.next().unwrap());
                let capacity = number(&mut tokens) as usize;
                let ascii: Vec<_> = raw.iter().map(|letter| b"ACGT"[*letter as usize]).collect();
                let mut builder = DatabaseBuilder::new();
                builder.add_sequence("subject", &ascii).unwrap();
                database = Some(builder.build().unwrap());
                pairs = vec![BlastOffsetPair::default(); capacity];
                next = 0;
            }
            "seed_batch" => {
                let previous = number(&mut tokens);
                let stop = number(&mut tokens);
                let count = number(&mut tokens);
                let expected: Vec<_> = (0..count)
                    .map(|_| BlastOffsetPair {
                        q_off: number(&mut tokens),
                        s_off: number(&mut tokens),
                    })
                    .collect();
                assert_eq!(next, previous, "line {line}");
                let hits = lookup
                    .as_ref()
                    .unwrap()
                    .scan_into(
                        database.as_ref().unwrap().sequence(0).unwrap(),
                        &mut next,
                        &mut pairs,
                    )
                    .unwrap();
                assert_eq!(hits, count as usize, "line {line}");
                assert_eq!(next, stop, "line {line}");
                assert_eq!(&pairs[..hits], expected, "line {line}");
            }
            "seed_done" => {
                assert_eq!(next, number(&mut tokens), "line {line}");
                let previous = next;
                assert_eq!(
                    lookup
                        .as_ref()
                        .unwrap()
                        .scan_into(
                            database.as_ref().unwrap().sequence(0).unwrap(),
                            &mut next,
                            &mut pairs
                        )
                        .unwrap(),
                    0
                );
                assert_eq!(next, previous);
            }
            _ => panic!("Unknown line {line}"),
        }
        assert_eq!(tokens.next(), None, "line {line}");
    }
}

#[test]
fn original_nucleotide_seed_orchestration() {
    verify(include_str!("fixtures/seed-lookup-2.17.0.txt"));
}

#[test]
#[ignore = "requires original C output via BLAST_SEED_LOOKUP_FIXTURE"]
fn freshly_compiled_nucleotide_seed_orchestration() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_SEED_LOOKUP_FIXTURE").unwrap()).unwrap());
}

#[test]
fn native_seed_validation_preserves_scan_state() {
    let query = vec![0; 64];
    let good = [SSeqRange { left: 0, right: 63 }];
    for regions in [
        vec![SSeqRange { left: -1, right: 3 }],
        vec![SSeqRange { left: 0, right: 64 }],
        vec![SSeqRange { left: 10, right: 5 }],
        vec![
            SSeqRange { left: 0, right: 31 },
            SSeqRange {
                left: 30,
                right: 60,
            },
        ],
    ] {
        assert_eq!(
            NucleotideSeedIndex::from_segments(&query, &regions, NucleotideTask::Blastn)
                .unwrap_err(),
            SeedError::InvalidRegions
        );
    }
    assert_eq!(
        NucleotideSeedIndex::from_segments(&[], &[], NucleotideTask::Blastn).unwrap_err(),
        SeedError::InvalidQuery
    );
    assert_eq!(
        NucleotideSeedIndex::from_segments(&[16], &[], NucleotideTask::Blastn).unwrap_err(),
        SeedError::InvalidQuery
    );
    let index = NucleotideSeedIndex::from_segments(&query, &good, NucleotideTask::Blastn).unwrap();
    let mut builder = DatabaseBuilder::new();
    builder.add_sequence("subject", &[b'A'; 64]).unwrap();
    let db = builder.build().unwrap();
    let subject = db.sequence(0).unwrap();
    let mut offset = 0;
    let mut small = vec![BlastOffsetPair::default(); index.minimum_capacity() - 1];
    assert_eq!(
        index.scan_into(subject, &mut offset, &mut small),
        Err(SeedError::InsufficientCapacity {
            minimum: index.minimum_capacity()
        })
    );
    assert_eq!(offset, 0);
    let mut buffer = vec![BlastOffsetPair::default(); index.recommended_capacity()];
    offset = 1;
    assert_eq!(
        index.scan_into(subject, &mut offset, &mut buffer),
        Err(SeedError::InvalidOffset)
    );
    assert_eq!(offset, 1);
}
