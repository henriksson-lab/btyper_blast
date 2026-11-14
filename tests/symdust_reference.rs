//! Symmetric DUST against the original class and unchanged nested algorithms.
use btyper_blast::dustmask::symdust::CSymDustMasker;
use btyper_blast::util::random_gen::CRandom;

fn bytes(text: &str) -> Vec<u8> {
    if text == "-" {
        return Vec::new();
    }
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
    let mut rng = CRandom::new_default();
    let mut reused = CSymDustMasker::new(20, 64, 1);
    let mut count = 0;
    for line in fixture.lines() {
        let fields: Vec<_> = line.split_whitespace().collect();
        match fields[0] {
            "default_rng" => {
                assert_eq!(rng.get_rand(), fields[2].parse::<u32>().unwrap(), "{line}")
            }
            "reset_rng" => {
                if fields[1] == "0" {
                    rng.reset();
                }
                assert_eq!(rng.get_rand(), fields[2].parse::<u32>().unwrap(), "{line}");
            }
            "after_seed_reset" => {
                rng.set_seed(4123);
                rng.reset();
                assert_eq!(rng.get_rand(), fields[1].parse::<u32>().unwrap());
            }
            "mask" => {
                let (level, window, linker) = (
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                );
                let mut masker = CSymDustMasker::new(level, window, linker);
                assert_eq!(
                    masker.mask(&bytes(fields[4])),
                    intervals(fields[5]),
                    "{line}"
                );
            }
            "reuse" => assert_eq!(
                reused.mask(&bytes(fields[2])),
                intervals(fields[3]),
                "{line}"
            ),
            "range" => {
                let mut masker = CSymDustMasker::new(20, 64, 1);
                assert_eq!(
                    masker.mask_range(
                        &bytes(fields[3]),
                        fields[1].parse().unwrap(),
                        fields[2].parse().unwrap()
                    ),
                    intervals(fields[4]),
                    "{line}"
                );
            }
            _ => panic!("Unknown reference record: {line}"),
        }
        count += 1;
    }
    count
}

#[test]
fn pinned_cpp_symmetric_dust() {
    assert_eq!(verify(include_str!("fixtures/symdust-2.17.0.txt")), 22_473);
}

#[test]
#[ignore = "requires a freshly compiled original class via BLAST_SYMDUST_FIXTURE"]
fn freshly_compiled_cpp_symmetric_dust() {
    let path = std::env::var_os("BLAST_SYMDUST_FIXTURE").expect("Set BLAST_SYMDUST_FIXTURE");
    let fixture = std::fs::read_to_string(path).unwrap();
    assert_eq!(verify(&fixture), 22_473);
}
