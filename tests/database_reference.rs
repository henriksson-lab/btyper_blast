//! Differential checks against unchanged selected NCBI C++ writer/RNG bodies.
use btyper_blast::seqdb::writedb_convert::{write_db_ncbi4na_to_binary, CAmbigDataBuilder};
use btyper_blast::util::random_gen::CRandom;

fn bytes(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
    assert_eq!(text.len() % 2, 0);
    (0..text.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&text[i..i + 2], 16).unwrap())
        .collect()
}

fn verify(fixture: &str) -> usize {
    let mut generator = None;
    let mut current_seed = 0;
    let mut cases = 0;
    for line in fixture.lines() {
        let fields: Vec<_> = line.split_whitespace().collect();
        match fields[0] {
            "rng" => {
                let seed = fields[1].parse::<u32>().unwrap();
                if fields[2] == "0" {
                    current_seed = seed;
                    generator = Some(CRandom::new(seed));
                }
                assert_eq!(current_seed, seed);
                assert_eq!(
                    generator.as_mut().unwrap().get_rand(),
                    fields[3].parse::<u32>().unwrap(),
                    "{line}"
                );
            }
            "reseed" => {
                let seed = fields[1].parse::<u32>().unwrap();
                generator.as_mut().unwrap().set_seed(seed);
                assert_eq!(
                    generator.as_mut().unwrap().get_rand(),
                    fields[2].parse::<u32>().unwrap(),
                    "{line}"
                );
            }
            "rng32" => {
                assert_eq!(current_seed, fields[1].parse::<u32>().unwrap());
                assert_eq!(
                    generator.as_mut().unwrap().x_get_rand32_bits(),
                    fields[3].parse::<u32>().unwrap(),
                    "{line}"
                );
            }
            "record" => {
                let sequence = bytes(fields[1]);
                let mut input = vec![0; (sequence.len() + 1) / 2];
                for (i, base) in sequence.iter().enumerate() {
                    input[i / 2] |= *base << if i % 2 == 0 { 4 } else { 0 };
                }
                let (mut packed, mut ambiguity) = (Vec::new(), Vec::new());
                write_db_ncbi4na_to_binary(
                    &input,
                    input.len() as i32,
                    sequence.len() as i32,
                    &mut packed,
                    &mut ambiguity,
                );
                assert_eq!(
                    packed,
                    bytes(fields[2]),
                    "packed length {} case {cases}",
                    sequence.len()
                );
                assert_eq!(
                    ambiguity,
                    bytes(fields[3]),
                    "ambiguity length {} case {cases}",
                    sequence.len()
                );
            }
            "large_ambiguity" => {
                let mut builder = CAmbigDataBuilder::new(fields[1].parse().unwrap());
                for offset in [0, 17, 0xfffffe] {
                    builder.check(15, offset);
                }
                let mut ambiguity = Vec::new();
                builder.get_ambig(&mut ambiguity);
                assert_eq!(ambiguity, bytes(fields[2]), "{line}");
            }
            _ => panic!("Unknown reference record {line}"),
        }
        cases += 1;
    }
    cases
}

#[test]
fn pinned_cpp_database_fixtures() {
    assert_eq!(verify(include_str!("fixtures/database-2.17.0.txt")), 77_795);
}

#[test]
#[ignore = "requires a freshly compiled C++ oracle via BLAST_DATABASE_FIXTURE"]
fn freshly_compiled_cpp_database_fixtures() {
    let path = std::env::var_os("BLAST_DATABASE_FIXTURE").expect("Set BLAST_DATABASE_FIXTURE");
    let fixture = std::fs::read_to_string(path).expect("Read original C++ oracle records");
    assert_eq!(verify(&fixture), 77_795);
}
