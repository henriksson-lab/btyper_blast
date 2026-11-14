// NCBI BLAST+ 2.17.0+, blast_util.c/h and blast_hits.c target translation.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_hits::BlastHSP;
use super::blast_program::EBlastProgramType;
use super::blast_util::{
    blast_frame_to_context, blast_get_translation, get_reverse_nucl_sequence, FENCE_SENTRY,
};

/// Source SBlastTargetTranslation, scoped to ordinary non-OOF tblastn.
/// The first input byte is the first NCBI4na nucleotide, without a sentinel.
#[derive(Debug)]
pub struct BlastTargetTranslation<'a> {
    pub gen_code_string: &'a [u8; 64],
    pub program_number: EBlastProgramType,
    pub partial: bool,
    pub num_frames: i32,
    pub translations: [Option<Vec<u8>>; 6],
    pub range: [i32; 12],
    pub subject: &'a [u8],
}
/// Defined slice representation of the source's virtual global-AA pointer.
/// A global residue i is residues[i-start_shift]; end is global exclusive.
/// Both leading and trailing source sentinel/FENCE bytes are excluded.
#[derive(Clone, Copy, Debug)]
pub struct BlastTargetTranslationView<'a> {
    pub residues: &'a [u8],
    pub start_shift: i32,
    pub end: i32,
}
pub fn blast_target_translation_new<'a>(
    subject: &'a [u8],
    gen_code_string: &'a [u8; 64],
    program_number: EBlastProgramType,
    is_ooframe: bool,
) -> BlastTargetTranslation<'a> {
    assert_eq!(program_number, EBlastProgramType::Tblastn);
    assert!(
        !is_ooframe,
        "out-of-frame target translation is outside standard tblastn"
    );
    assert!(subject.len() <= i32::MAX as usize);
    BlastTargetTranslation {
        gen_code_string,
        program_number,
        partial: !is_ooframe,
        num_frames: 6,
        translations: std::array::from_fn(|_| None),
        range: [0; 12],
        subject,
    }
}
pub fn blast_target_translation_free(
    target: Option<BlastTargetTranslation<'_>>,
) -> Option<BlastTargetTranslation<'_>> {
    drop(target);
    None
}
pub fn blast_hsp_get_target_translation<'t>(
    target: &'t mut BlastTargetTranslation<'_>,
    hsp: Option<&BlastHSP>,
) -> Option<BlastTargetTranslationView<'t>> {
    let hsp = hsp?;
    let context = blast_frame_to_context(hsp.subject.frame, target.program_number);
    assert!((0..6).contains(&context));
    let context = context as usize;
    let start = target.range[2 * context];
    let stop = target.range[2 * context + 1];
    let subject_length = target.subject.len() as i32;
    if target.partial && (start != 0 || stop < subject_length / 3 - 3) {
        let (nucl_start, mut nucl_end) = if hsp.subject.offset < 0 {
            (0, subject_length)
        } else {
            (
                (3 * hsp.subject.offset - 99).max(0),
                subject_length.min(3 * hsp.subject.end + 99),
            )
        };
        if hsp.subject.offset >= 0 && subject_length - nucl_end <= 21 {
            nucl_end = subject_length;
        }
        let nucl_length = nucl_end - nucl_start;
        let translation_length = 1 + nucl_length / 3;
        let start_shift = nucl_start / 3;
        let nucl_shift = if hsp.subject.frame < 0 {
            subject_length - nucl_start - nucl_length
        } else {
            nucl_start
        };
        if start_shift < start || start_shift + translation_length > stop {
            target.range[2 * context] = start_shift;
            if translation_length > stop - start {
                target.translations[context] = Some(vec![0; (translation_length + 2) as usize]);
            }
            let nucl_seq = &target.subject[nucl_shift as usize..];
            let nucl_seq_rev = if hsp.subject.frame < 0 {
                get_reverse_nucl_sequence(nucl_seq, nucl_length)
            } else {
                Vec::new()
            };
            let translation = target.translations[context]
                .as_mut()
                .expect("source target buffer allocated before translation");
            let length = blast_get_translation(
                nucl_seq,
                &nucl_seq_rev,
                nucl_length,
                hsp.subject.frame,
                translation,
                target.gen_code_string,
            );
            target.range[2 * context + 1] = start_shift + length;
            if hsp.subject.offset >= 0 {
                translation[0] = FENCE_SENTRY;
                translation[length as usize + 1] = FENCE_SENTRY;
            }
        }
    }
    let start_shift = target.range[2 * context];
    let end = target.range[2 * context + 1];
    // Source returns an invalid pointer when an uninitialized tiny subject
    // never enters the translation branch. There is no ordinary HSP there;
    // an absent owned buffer is represented as None without pointer creation.
    let translation = target.translations[context].as_ref()?;
    Some(BlastTargetTranslationView {
        residues: &translation[1..(end - start_shift + 1) as usize],
        start_shift,
        end,
    })
}
