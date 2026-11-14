use btyper_blast::{NucleotideTask, PreparedNucleotideQuery, QueryError};

fn bytes(text: &str) -> Vec<u8> {
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn intervals(text: &str) -> Vec<(u32, u32)> {
    if text == "-" {
        return Vec::new();
    }
    text.split(',')
        .map(|pair| {
            let (start, stop) = pair.split_once(':').unwrap();
            (start.parse().unwrap(), stop.parse().unwrap())
        })
        .collect()
}

fn verify(fixture: &str) -> usize {
    let mut count = 0;
    for line in fixture.lines() {
        let fields: Vec<_> = line.split_whitespace().collect();
        assert_eq!(fields[0], "query");
        let task = match fields[1] {
            "megablast" => NucleotideTask::Megablast,
            "blastn" => NucleotideTask::Blastn,
            _ => panic!("Unknown task"),
        };
        let input = bytes(fields[3]);
        let query = PreparedNucleotideQuery::new(fields[2], &input, task).unwrap();
        assert_eq!(query.id(), fields[2]);
        assert_eq!(query.task(), task);
        assert_eq!(
            query.ncbi4na(),
            bytes(fields[4]),
            "encoding {} {}",
            fields[1],
            fields[2]
        );
        assert_eq!(
            query.dust_masks(),
            intervals(fields[5]),
            "mask {} {}",
            fields[1],
            fields[2]
        );
        let last = query.len() as u32 - 1;
        let expected_reverse: Vec<_> = intervals(fields[5])
            .into_iter()
            .rev()
            .map(|(start, stop)| (last - stop, last - start))
            .collect();
        assert_eq!(query.reverse_dust_masks(), expected_reverse);
        assert!(query.soft_masking());
        // Soft masking must leave both original residue contexts intact.
        let table = [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];
        let forward: Vec<_> = bytes(fields[4])
            .iter()
            .map(|base| table[*base as usize])
            .collect();
        let reverse: Vec<_> = bytes(fields[4])
            .iter()
            .rev()
            .map(|base| {
                let complement =
                    ((base & 1) << 3) | ((base & 2) << 1) | ((base & 4) >> 1) | ((base & 8) >> 3);
                table[complement as usize]
            })
            .collect();
        assert_eq!(query.forward_blastna(), forward);
        assert_eq!(query.reverse_blastna(), reverse);
        let expected_hyphens: Vec<_> = input
            .iter()
            .enumerate()
            .filter_map(|(i, base)| (*base == b'-').then_some(i))
            .collect();
        assert_eq!(query.ignored_hyphens(), expected_hyphens);
        count += 1;
    }
    count
}

#[test]
fn query_setup_matches_original_blast_archives() {
    assert_eq!(verify(include_str!("fixtures/query-2.17.0.txt")), 414);
}

#[test]
#[ignore = "requires fresh original BLAST archives via BLAST_QUERY_FIXTURE"]
fn freshly_generated_original_blast_query_setup() {
    let fixture = std::fs::read_to_string(
        std::env::var_os("BLAST_QUERY_FIXTURE").expect("Set BLAST_QUERY_FIXTURE"),
    )
    .unwrap();
    assert_eq!(verify(&fixture), 414);
}

#[test]
fn query_input_errors_and_task_distinction() {
    let default_task = NucleotideTask::default();
    assert_eq!(default_task, NucleotideTask::Megablast);
    assert_ne!(default_task, NucleotideTask::Blastn);
    assert!(matches!(
        PreparedNucleotideQuery::new("", b"ACGT", default_task),
        Err(QueryError::EmptyIdentifier)
    ));
    assert!(matches!(
        PreparedNucleotideQuery::new("empty", b"--", default_task),
        Err(QueryError::EmptySequence)
    ));
    assert!(matches!(
        PreparedNucleotideQuery::new("bad", b"AC?", default_task),
        Err(QueryError::InvalidResidue {
            offset: 2,
            byte: b'?'
        })
    ));
    assert!(matches!(
        PreparedNucleotideQuery::from_ncbi4na("bad", &[1, 16], default_task),
        Err(QueryError::InvalidNcbi4na {
            offset: 1,
            value: 16
        })
    ));
    let sequence = [1, 0, 2, 8];
    let query = PreparedNucleotideQuery::from_ncbi4na("raw", &sequence, default_task).unwrap();
    assert_eq!(query.ncbi4na(), &[1, 15, 2, 8]);
    assert_eq!(query.forward_blastna(), &[0, 14, 1, 3]);
    assert_eq!(query.reverse_blastna(), &[0, 2, 14, 3]);
    assert!(query.ignored_hyphens().is_empty());
}
