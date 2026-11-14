//! Differential checks against fixtures emitted by the unmodified pinned C.
//! Set BLAST_KERNEL_FIXTURE to exercise a freshly compiled C reference as well.
use btyper_blast::core::blast_encoding::*;
use btyper_blast::core::blast_program::*;
use btyper_blast::core::blast_util::*;
use btyper_blast::core::gapinfo::*;
use btyper_blast::core::genetic_code::STANDARD_GENETIC_CODE;
use btyper_blast::core::ncbi_math::*;

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

fn int(text: &str) -> i32 {
    text.parse().unwrap()
}

fn float(text: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(text, 16).unwrap())
}

fn program(text: &str) -> EBlastProgramType {
    use EBlastProgramType::*;
    match int(text) {
        0 => Undefined,
        3 => Blastp,
        12 => Blastn,
        22 => Blastx,
        41 => Tblastn,
        60 => Tblastx,
        67 => PsiBlast,
        105 => PsiTblastn,
        131 => RpsBlast,
        150 => RpsTblastn,
        259 => PhiBlastp,
        268 => PhiBlastn,
        524 => Mapping,
        _ => panic!("Unknown reference program"),
    }
}

fn operation(text: &str) -> EGapAlignOpType {
    use EGapAlignOpType::*;
    match int(text) {
        0 => Del,
        1 => Del2,
        2 => Del1,
        3 => Sub,
        4 => Ins1,
        5 => Ins2,
        6 => Ins,
        7 => Decline,
        8 => Invalid,
        _ => panic!("Unknown reference edit operation"),
    }
}

fn check_float(actual: f64, expected: f64, context: &str) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "{context}: expected NaN, got {actual}");
    } else if expected == 0.0 || !expected.is_finite() {
        assert_eq!(actual.to_bits(), expected.to_bits(), "{context}");
    } else {
        // C exp/log implementations can differ slightly across OS libm versions.
        // Polynomial, integer, and sequence outputs use exact comparisons below.
        assert_eq!(
            actual.is_sign_negative(),
            expected.is_sign_negative(),
            "{context}"
        );
        assert!(
            actual.to_bits().abs_diff(expected.to_bits()) <= 8,
            "{context}: actual {actual:?}, C {expected:?}"
        );
    }
}

fn verify(fixture: &str) -> usize {
    let mut count = 0;
    let mut prelim = Some(gap_prelim_edit_block_new());
    let mut next_prelim_step = 0;
    for line in fixture.lines() {
        let fields: Vec<_> = line.split_whitespace().collect();
        match fields[0] {
            "table" => {
                let actual: &[u8] = match fields[1] {
                    "NCBI4NA_TO_BLASTNA" => &NCBI4NA_TO_BLASTNA,
                    "BLASTNA_TO_NCBI4NA" => &BLASTNA_TO_NCBI4NA,
                    "BLASTNA_TO_IUPACNA" => &BLASTNA_TO_IUPACNA,
                    "NCBI4NA_TO_IUPACNA" => &NCBI4NA_TO_IUPACNA,
                    "IUPACNA_TO_BLASTNA" => &IUPACNA_TO_BLASTNA,
                    "IUPACNA_TO_NCBI4NA" => &IUPACNA_TO_NCBI4NA,
                    "AMINOACID_TO_NCBISTDAA" => &AMINOACID_TO_NCBISTDAA,
                    "NCBISTDAA_TO_AMINOACID" => &NCBISTDAA_TO_AMINOACID,
                    "codes" => &STANDARD_GENETIC_CODE,
                    _ => panic!("Unexpected table"),
                };
                assert_eq!(actual, bytes(fields[2]), "{line}");
            }
            "sentinels" => {
                assert_eq!(i32::from(K_PROT_SENTINEL), int(fields[1]));
                assert_eq!(i32::from(K_NUCL_SENTINEL), int(fields[2]));
            }
            "packed_table" => {
                let table =
                    s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), fields[1] == "1")
                        .unwrap();
                assert_eq!(table.as_slice(), bytes(fields[2]));
            }
            "packed_table_null" => assert_eq!(
                s_blast_get_translation_table(None, false).is_none() as i32,
                int(fields[1])
            ),
            "protein_length" => {
                for (context, value) in fields[2..].iter().enumerate() {
                    assert_eq!(
                        blast_get_translated_protein_length(
                            int(fields[1]) as usize,
                            context as u32
                        ),
                        int(value) as usize,
                        "{line}"
                    );
                }
            }
            "compressed_null" => {
                let table =
                    s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), false).unwrap();
                let mut output = [0xee; 8];
                assert_eq!(
                    blast_translate_compressed_sequence(&table, 12, None, 1, Some(&mut output)),
                    int(fields[1])
                );
                assert_eq!(
                    blast_translate_compressed_sequence(&table, 12, Some(&[0; 4]), 1, None),
                    int(fields[2])
                );
                assert_eq!(i32::from(output[0]), int(fields[3]));
            }
            "compressed" => {
                let frame = int(fields[1]) as i16;
                let length = int(fields[2]);
                let sequence = bytes(fields[3]);
                let table =
                    s_blast_get_translation_table(Some(&STANDARD_GENETIC_CODE), frame < 0).unwrap();
                let mut output = vec![0xee; (length / 3 + 4) as usize];
                let protein_length = blast_translate_compressed_sequence(
                    &table,
                    length,
                    Some(&sequence),
                    frame,
                    Some(&mut output),
                );
                assert_eq!(
                    protein_length,
                    int(fields[4]),
                    "compressed frame {frame} length {length}"
                );
                assert_eq!(
                    &output[..protein_length as usize + 2],
                    bytes(fields[5]),
                    "compressed frame {frame} length {length}"
                );
                assert!(output[protein_length as usize + 2..]
                    .iter()
                    .all(|b| *b == 0xee));
            }
            "alltranslations" => {
                let encoding = if fields[1] == "2" {
                    EBlastEncoding::Ncbi4na
                } else {
                    EBlastEncoding::Ncbi2na
                };
                let length = int(fields[2]);
                let sequence = bytes(fields[3]);
                let mut all = Vec::new();
                let mut offsets = [0; 7];
                let mut mixed = Vec::new();
                assert_eq!(
                    blast_get_all_translations(
                        &sequence,
                        encoding,
                        length,
                        &STANDARD_GENETIC_CODE,
                        Some(&mut all),
                        Some(&mut offsets),
                        Some(&mut mixed)
                    ),
                    0
                );
                for (context, offset) in offsets.iter().enumerate() {
                    assert_eq!(
                        *offset,
                        int(fields[4 + context]) as u32,
                        "alltranslations enc {} length {length}",
                        fields[1]
                    );
                }
                assert_eq!(
                    &all[..offsets[6] as usize],
                    bytes(fields[11]),
                    "alltranslations enc {} length {length}",
                    fields[1]
                );
                assert_eq!(
                    mixed,
                    bytes(fields[12]),
                    "mixed translation enc {} length {length}",
                    fields[1]
                );
            }
            "partial" => {
                let frame = int(fields[1]) as i16;
                let use_mixed = fields[2] == "1";
                let sequence = bytes(fields[3]);
                let mut partial = Vec::new();
                let mut mixed = Vec::new();
                let mut protein_length = -1;
                assert_eq!(
                    blast_get_partial_translation(
                        &sequence,
                        sequence.len() as i32,
                        frame,
                        &STANDARD_GENETIC_CODE,
                        Some(&mut partial),
                        Some(&mut protein_length),
                        if use_mixed { Some(&mut mixed) } else { None },
                    ),
                    0
                );
                assert_eq!(
                    protein_length,
                    int(fields[4]),
                    "partial frame {frame} mixed {use_mixed} length {}",
                    sequence.len()
                );
                let defined_length = if use_mixed {
                    sequence.len() + 2
                } else {
                    protein_length as usize + 2
                };
                assert_eq!(
                    &partial[..defined_length],
                    bytes(fields[5]),
                    "partial frame {frame} mixed {use_mixed} length {}",
                    sequence.len()
                );
                if use_mixed {
                    assert_eq!(&mixed[..sequence.len() + 1], bytes(fields[6]));
                } else {
                    assert_eq!(mixed, bytes(fields[6]));
                }
            }
            "program" => {
                let p = program(fields[1]);
                let flags = blast_query_is_protein(p) as i32
                    | ((blast_query_is_nucleotide(p) as i32) << 1)
                    | ((blast_subject_is_protein(p) as i32) << 2)
                    | ((blast_subject_is_nucleotide(p) as i32) << 3)
                    | ((blast_query_is_translated(p) as i32) << 4)
                    | ((blast_subject_is_translated(p) as i32) << 5)
                    | ((blast_program_is_nucleotide(p) as i32) << 6)
                    | ((blast_program_is_valid(p) as i32) << 7);
                assert_eq!(flags, int(fields[2]), "{line}");
                assert_eq!(blast_get_number_of_contexts(p), int(fields[3]) as u32);
                for (j, value) in fields[4..].iter().enumerate() {
                    assert_eq!(
                        i32::from(blast_context_to_frame(p, j as u32)),
                        int(value),
                        "{line}"
                    );
                }
            }
            "frame" => assert_eq!(
                blast_frame_to_context(int(fields[2]) as i16, program(fields[1])),
                int(fields[3]),
                "{line}"
            ),
            "codon" => {
                let codon = [
                    int(fields[1]) as u8,
                    int(fields[2]) as u8,
                    int(fields[3]) as u8,
                ];
                assert_eq!(
                    i32::from(s_codon_to_aa(&codon, &STANDARD_GENETIC_CODE)),
                    int(fields[4]),
                    "{line}"
                );
            }
            "pack" => {
                let sequence = bytes(fields[2]);
                let encoding = if fields[1] == "1" {
                    EBlastEncoding::Nucleotide
                } else {
                    EBlastEncoding::Ncbi4na
                };
                assert_eq!(
                    blast_pack_dna(&sequence, sequence.len() as i32, encoding),
                    bytes(fields[3]),
                    "pack {} len {}",
                    fields[1],
                    sequence.len()
                );
            }
            "reverse" => {
                let sequence = bytes(fields[1]);
                assert_eq!(
                    get_reverse_nucl_sequence(&sequence, sequence.len() as i32),
                    bytes(fields[2]),
                    "reverse len {}",
                    sequence.len()
                );
            }
            "translation" => {
                let frame = int(fields[1]) as i16;
                let sequence = bytes(fields[2]);
                let reverse = get_reverse_nucl_sequence(&sequence, sequence.len() as i32);
                let mut protein = vec![255; sequence.len() / 3 + 2];
                let length = blast_get_translation(
                    &sequence,
                    &reverse,
                    sequence.len() as i32,
                    frame,
                    &mut protein,
                    &STANDARD_GENETIC_CODE,
                );
                assert_eq!(
                    length,
                    int(fields[3]),
                    "translation frame {frame} len {}",
                    sequence.len()
                );
                assert_eq!(
                    &protein[..length as usize + 2],
                    bytes(fields[4]),
                    "translation frame {frame} len {}",
                    sequence.len()
                );
                // The original writes only the returned sequence and sentinels.
                assert!(protein[length as usize + 2..].iter().all(|b| *b == 255));
            }
            "edit" => {
                let script = gap_edit_script_new(int(fields[1]));
                assert_eq!(script.is_some() as i32, int(fields[2]));
                if let Some(script) = &script {
                    assert_eq!(script.size, int(fields[3]));
                    for (i, pair) in fields[4..fields.len() - 1].iter().enumerate() {
                        let (op, number) = pair.split_once(':').unwrap();
                        assert_eq!(script.op_type[i], operation(op));
                        assert_eq!(script.num[i], int(number));
                    }
                }
                assert_eq!(
                    gap_edit_script_delete(script).is_none() as i32,
                    int(fields.last().unwrap())
                );
            }
            "prelim" => {
                let block = prelim.as_mut().unwrap();
                let final_step = int(fields[1]);
                while next_prelim_step <= final_step {
                    let i = next_prelim_step;
                    let op = match (i / 3) % 3 {
                        0 => EGapAlignOpType::Del,
                        1 => EGapAlignOpType::Sub,
                        _ => EGapAlignOpType::Ins,
                    };
                    let count = if i % 5 == 0 { 0 } else { i % 11 };
                    if i == 900 {
                        gap_prelim_edit_block_reset(Some(block));
                    } else {
                        gap_prelim_edit_block_add(block, op, count);
                    }
                    next_prelim_step += 1;
                }
                assert_eq!(block.num_ops_allocated, int(fields[2]));
                assert_eq!(block.num_ops, int(fields[3]));
                assert_eq!(block.last_op, operation(fields[4]));
                assert_eq!(block.num_ops as usize, fields.len() - 5);
                for (i, pair) in fields[5..].iter().enumerate() {
                    let (op, number) = pair.split_once(':').unwrap();
                    assert_eq!(block.edit_ops[i].op_type, operation(op));
                    assert_eq!(block.edit_ops[i].num, int(number));
                }
            }
            "free" => assert_eq!(
                gap_prelim_edit_block_free(prelim.take()).is_none() as i32,
                int(fields[1])
            ),
            "expm1" | "log1p" | "powi" => {
                let x = float(fields[1]);
                let (actual, expected) = match fields[0] {
                    "expm1" => (blast_expm1(x), float(fields[2])),
                    "log1p" => (blast_log1p(x), float(fields[2])),
                    _ => (blast_powi(x, int(fields[2])), float(fields[3])),
                };
                if (fields[0] == "expm1" && x.abs() <= 0.33)
                    || (fields[0] == "log1p" && x.abs() < 0.2)
                    || fields[0] == "powi"
                {
                    assert_eq!(actual.to_bits(), expected.to_bits(), "{line}");
                } else {
                    check_float(actual, expected, line);
                }
            }
            "nint" => assert_eq!(
                blast_nint(float(fields[1])) as i64,
                fields[2].parse::<i64>().unwrap(),
                "{line}"
            ),
            "gcd" => assert_eq!(
                blast_gcd(int(fields[1]), int(fields[2])),
                int(fields[3]),
                "{line}"
            ),
            "gdb3" => {
                let (mut a, mut b, mut c) = (int(fields[1]), int(fields[2]), int(fields[3]));
                assert_eq!(blast_gdb3(&mut a, &mut b, &mut c), int(fields[4]), "{line}");
                assert_eq!(
                    (a, b, c),
                    (int(fields[5]), int(fields[6]), int(fields[7])),
                    "{line}"
                );
            }
            _ => panic!("Unexpected C reference record: {line}"),
        }
        count += 1;
    }
    count
}

#[test]
fn pinned_c_kernel_fixtures() {
    let count = verify(include_str!("fixtures/kernels-2.17.0.txt"));
    assert_eq!(count, 19_322, "Incomplete reference fixture");
}

#[test]
#[ignore = "requires a freshly compiled C oracle via BLAST_KERNEL_FIXTURE"]
fn freshly_compiled_c_kernel_fixtures() {
    let path = std::env::var_os("BLAST_KERNEL_FIXTURE").expect("Set BLAST_KERNEL_FIXTURE");
    let fixture = std::fs::read_to_string(path).expect("Read fresh C oracle output");
    let count = verify(&fixture);
    assert_eq!(count, 19_322);
}
