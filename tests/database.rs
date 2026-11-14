use btyper_blast::{DatabaseBuilder, DatabaseError, NucleotideDatabase};
use std::io::Cursor;

fn bytes(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn verify_makeblastdb(fixture: &str) {
    let mut fasta = Vec::new();
    let lines: Vec<_> = fixture
        .lines()
        .map(|line| line.split_whitespace().collect::<Vec<_>>())
        .collect();
    assert_eq!(lines[0][0], "stats");
    for fields in &lines[1..] {
        assert_eq!(fields[0], "sequence");
        fasta.extend_from_slice(format!(">{} description\r\n", fields[1]).as_bytes());
        fasta.extend_from_slice(&bytes(fields[2]));
        fasta.extend_from_slice(b"\r\n");
    }
    let database = NucleotideDatabase::from_fasta(Cursor::new(fasta)).unwrap();
    assert_eq!(database.len(), lines[0][1].parse::<usize>().unwrap());
    assert_eq!(database.total_length(), lines[0][2].parse::<u64>().unwrap());
    assert_eq!(database.max_length(), lines[0][3].parse::<usize>().unwrap());
    assert!(!database.warnings().is_empty()); // The original input includes invalid X, ?, and .
    for (ordinal, fields) in lines[1..].iter().enumerate() {
        let subject = database.sequence(ordinal).unwrap();
        assert_eq!(subject.id(), fields[1]);
        assert_eq!(
            subject.ncbi4na(),
            bytes(fields[3]),
            "canonical {}",
            fields[1]
        );
        assert_eq!(
            subject.packed_ncbi2na(),
            bytes(fields[4]),
            "packed {}",
            fields[1]
        );
        assert_eq!(
            subject.ambiguity_data(),
            bytes(fields[5]),
            "ambiguity {}",
            fields[1]
        );
        assert_eq!(subject.len(), bytes(fields[3]).len());
    }
    assert!(database.sequence(database.len()).is_none());
}

#[test]
fn memory_database_matches_original_makeblastdb_payloads() {
    verify_makeblastdb(include_str!("fixtures/makeblastdb-2.17.0.txt"));
}

#[test]
#[ignore = "requires pinned makeblastdb via BLAST_MAKEBLASTDB_FIXTURE"]
fn freshly_generated_makeblastdb_payloads() {
    let fixture = std::fs::read_to_string(
        std::env::var_os("BLAST_MAKEBLASTDB_FIXTURE").expect("Set BLAST_MAKEBLASTDB_FIXTURE"),
    )
    .unwrap();
    verify_makeblastdb(&fixture);
}

#[test]
fn builder_rejects_bad_input_without_adding_partial_records() {
    let mut builder = DatabaseBuilder::new();
    builder.add_sequence("first", b"ACGT").unwrap();
    assert!(matches!(
        builder.add_sequence("bad", b"AC?T"),
        Err(DatabaseError::InvalidResidue {
            offset: 2,
            byte: b'?'
        })
    ));
    assert!(matches!(
        builder.add_ncbi4na("bad", &[1, 201]),
        Err(DatabaseError::InvalidNcbi4na {
            offset: 1,
            value: 201
        })
    ));
    assert!(matches!(
        builder.add_sequence("", b"ACGT"),
        Err(DatabaseError::EmptyIdentifier)
    ));
    assert!(matches!(
        builder.add_sequence("empty", b""),
        Err(DatabaseError::EmptySequence)
    ));
    builder.add_ncbi4na("first", &[1, 0, 2]).unwrap();
    let db = builder.build().unwrap();
    assert_eq!(db.len(), 2);
    assert_eq!(db.total_length(), 7);
    assert_eq!(db.max_length(), 4);
    assert_eq!(db.sequence(1).unwrap().ncbi4na(), &[1, 15, 2]);
    assert_eq!(db.sequence(1).unwrap().blastna(), &[0, 14, 1]);
}

#[test]
fn fasta_reports_omitted_bytes_without_printing_or_losing_locations() {
    let fasta = b"; comment\n>one title\nAu?C\n>two\nA.XC\n";
    let db = NucleotideDatabase::from_fasta(Cursor::new(fasta)).unwrap();
    assert_eq!(db.sequence(0).unwrap().ncbi4na(), &[1, 8, 2]);
    assert_eq!(db.sequence(1).unwrap().ncbi4na(), &[1, 2]);
    assert_eq!(db.warnings().len(), 3);
    assert_eq!(
        (
            db.warnings()[0].sequence_ordinal,
            db.warnings()[0].line,
            db.warnings()[0].column,
            db.warnings()[0].byte
        ),
        (0, 3, 3, b'?')
    );
    assert!(matches!(
        NucleotideDatabase::from_fasta(Cursor::new(b"ACGT")),
        Err(DatabaseError::MissingFastaHeader { line: 1 })
    ));
    assert!(matches!(
        NucleotideDatabase::from_fasta(Cursor::new(b"")),
        Err(DatabaseError::EmptyDatabase)
    ));
    assert!(matches!(
        NucleotideDatabase::from_fasta(Cursor::new(b">empty\n>other\nA\n")),
        Err(DatabaseError::EmptySequence)
    ));
}
