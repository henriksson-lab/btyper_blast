// NCBI BLAST+ 2.17.0+, selected ordinary protein preliminary orchestration.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_dynamic_align::{blast_semi_gapped_align, DynamicSequence};
use super::blast_extend::BlastInitHSP;
use super::blast_gapalign::BlastGapAlignStruct;
use super::blast_parameters::BlastScoringParameters;
use super::blast_program::EBlastProgramType;

pub const MAX_SUBJECT_OFFSET: i32 = 90_000;
pub const MAX_TOTAL_GAPS: i32 = 3000;

pub fn adjust_subject_range(
    subject_offset: &mut i32,
    subject_length: &mut i32,
    query_offset: i32,
    query_length: i32,
    start_shift: &mut i32,
) {
    let length = *subject_length;
    if length < MAX_SUBJECT_OFFSET {
        *start_shift = 0;
        return;
    }
    let offset = *subject_offset;
    let max_extension_left = query_offset + MAX_TOTAL_GAPS;
    let max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;
    if offset <= max_extension_left {
        *start_shift = 0;
    } else {
        *start_shift = offset - max_extension_left;
        *subject_offset = max_extension_left;
    }
    *subject_length = length.min(offset + max_extension_right) - *start_shift;
}

/// Standard in-frame, unrestricted protein preliminary extension. The source
/// sequence pointers are represented by signed origins into sentinel buffers.
/// Out-of-frame and restricted approximate DP are outside default tblastn.
pub fn s_blast_prot_gapped_alignment(
    _program: EBlastProgramType,
    query: &DynamicSequence<'_>,
    query_length: i32,
    subject: &DynamicSequence<'_>,
    mut subject_length: i32,
    gap_align: &mut BlastGapAlignStruct<'_>,
    score_params: &BlastScoringParameters,
    init_hsp: &BlastInitHSP,
    restricted_alignment: bool,
    mut fence_hit: Option<&mut bool>,
) -> i16 {
    assert!(
        !score_params.options.is_ooframe,
        "out-of-frame alignment is outside standard scope"
    );
    assert!(
        !restricted_alignment,
        "restricted approximate DP is outside default tblastn scope"
    );
    let q_length = init_hsp.offsets.q_off + 1;
    let mut s_length = init_hsp.offsets.s_off + 1;
    let mut subject_shift = 0;
    adjust_subject_range(
        &mut s_length,
        &mut subject_length,
        q_length,
        query_length,
        &mut subject_shift,
    );
    let mut found_start = false;
    let mut found_end = false;
    let mut score_left = 0;
    let mut private_q_start = 0;
    let mut private_s_start = 0;
    if q_length != 0 && s_length != 0 {
        found_start = true;
        score_left = blast_semi_gapped_align(
            query,
            &subject.shifted(subject_shift),
            q_length,
            s_length,
            &mut private_q_start,
            &mut private_s_start,
            true,
            None,
            &mut gap_align.dynamic,
            score_params,
            init_hsp.offsets.q_off,
            false,
            true,
            fence_hit.as_deref_mut(),
        );
        gap_align.dynamic.query_start = q_length - private_q_start;
        gap_align.dynamic.subject_start = s_length - private_s_start + subject_shift;
    }
    let mut score_right = 0;
    if q_length < query_length && s_length < subject_length {
        found_end = true;
        let mut q_stop = gap_align.dynamic.query_stop;
        let mut s_stop = gap_align.dynamic.subject_stop;
        score_right = blast_semi_gapped_align(
            &query.shifted(init_hsp.offsets.q_off),
            &subject.shifted(init_hsp.offsets.s_off),
            query_length - q_length,
            subject_length - s_length,
            &mut q_stop,
            &mut s_stop,
            true,
            None,
            &mut gap_align.dynamic,
            score_params,
            init_hsp.offsets.q_off,
            false,
            false,
            fence_hit.as_deref_mut(),
        );
        gap_align.dynamic.query_stop = q_stop + init_hsp.offsets.q_off + 1;
        gap_align.dynamic.subject_stop = s_stop + init_hsp.offsets.s_off + 1;
    }
    if !found_start {
        gap_align.dynamic.query_start = init_hsp.offsets.q_off;
        gap_align.dynamic.subject_start = init_hsp.offsets.s_off;
    }
    if !found_end {
        gap_align.dynamic.query_stop = init_hsp.offsets.q_off;
        gap_align.dynamic.subject_stop = init_hsp.offsets.s_off;
    }
    gap_align.dynamic.score = score_right + score_left;
    0
}
