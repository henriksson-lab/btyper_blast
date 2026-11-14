// NCBI BLAST+ 2.17.0+, selected nucleotide blast_traceback.c orchestration.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_dynamic_align::DynamicSequence;
use super::blast_gapalign::{
    blast_gapped_alignment_with_traceback, blast_greedy_gapped_alignment, BlastGapAlignStruct,
};
use super::blast_hits::*;
use super::blast_hsp_evaluate::*;
use super::blast_hsp_link::blast_link_hsps;
use super::blast_options::{BlastExtensionOptions, EBlastTbackExt};
use super::blast_parameters::{BlastHitSavingParameters, BlastScoringParameters};
use super::blast_prelim::{GapSequence, GapSubject, HSP_MAX_WINDOW};
use super::blast_program::{blast_program_is_rps_blast, EBlastProgramType};
use super::blast_protein_gapalign::adjust_subject_range;
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::BlastScoreBlk;
use super::interval_tree::*;
use std::ptr::NonNull;

pub fn blast_hsp_update_with_traceback(
    gap_align: Option<&mut BlastGapAlignStruct<'_>>,
    hsp: Option<&mut BlastHSP>,
) -> i16 {
    let (Some(gap), Some(hsp)) = (gap_align, hsp) else {
        return -1;
    };
    hsp.score = gap.dynamic.score;
    hsp.query.offset = gap.dynamic.query_start;
    hsp.subject.offset = gap.dynamic.subject_start;
    hsp.query.end = gap.dynamic.query_stop;
    hsp.subject.end = gap.dynamic.subject_stop;
    if let Some(script) = gap.edit_script.take() {
        hsp.gap_info = Some(Box::new(script));
    }
    0
}

pub fn s_hsp_list_rescale_scores(list: &mut BlastHSPList, scale_factor: f64) {
    for hsp in &mut list.hsp_array[..list.hspcnt as usize] {
        let hsp = hsp.as_mut().unwrap();
        hsp.score = ((hsp.score as f64 + 0.5 * scale_factor) / scale_factor) as i32;
    }
    blast_hsp_list_sort_by_score(Some(list));
}

pub fn s_hsp_list_post_traceback_update(
    program: EBlastProgramType,
    list: &mut BlastHSPList,
    query_info: &BlastQueryInfo,
    score_params: &BlastScoringParameters,
    hit_params: &BlastHitSavingParameters,
    sbp: &mut BlastScoreBlk,
    subject_length: i32,
) -> i16 {
    assert!(
        !blast_program_is_rps_blast(program),
        "RPS traceback is outside standard scope"
    );
    let gapped = score_params.options.gapped_calculation;
    if let Some(link) = &hit_params.link_hsp_params {
        blast_link_hsps(
            program,
            Some(list),
            query_info,
            subject_length,
            sbp,
            link,
            gapped,
        );
    } else {
        blast_hsp_list_get_evalues(
            program,
            query_info,
            subject_length,
            Some(list),
            gapped,
            false,
            sbp,
            0.0,
            1.0,
        );
    }
    blast_hsp_list_reap_by_evalue(
        Some(list),
        hit_params.options.as_ref().unwrap().expect_value,
    );
    s_hsp_list_rescale_scores(list, score_params.scale_factor);
    blast_hsp_list_get_bit_scores(Some(list), gapped, sbp);
    0
}

pub fn blast_get_offsets_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    matrix: &[Vec<i32>],
    pssm: Option<&[Vec<i32>]>,
    hsp: &BlastHSP,
    q_retval: &mut i32,
    s_retval: &mut i32,
) -> bool {
    let q_length = hsp.query.end - hsp.query.offset;
    let s_length = hsp.subject.end - hsp.subject.offset;
    let q_start = hsp.query.offset;
    let s_start = hsp.subject.offset;
    if q_length <= HSP_MAX_WINDOW {
        *q_retval = q_start + q_length / 2;
        *s_retval = s_start + q_length / 2;
        return true;
    }
    let mut hsp_end = q_start + HSP_MAX_WINDOW;
    let mut q = q_start;
    let mut s = s_start;
    let mut score = 0;
    for index in q_start..hsp_end {
        score += if let Some(pssm) = pssm {
            pssm[index as usize][subject[s as usize] as usize]
        } else {
            matrix[query[q as usize] as usize][subject[s as usize] as usize]
        };
        q += 1;
        s += 1;
    }
    let mut max_score = score;
    let mut max_offset = hsp_end - 1;
    hsp_end = q_start + q_length.min(s_length);
    for index in q_start + HSP_MAX_WINDOW..hsp_end {
        if let Some(pssm) = pssm {
            score -= pssm[(index - HSP_MAX_WINDOW) as usize]
                [subject[(s - HSP_MAX_WINDOW) as usize] as usize];
            score += pssm[index as usize][subject[s as usize] as usize];
        } else {
            score -= matrix[query[(q - HSP_MAX_WINDOW) as usize] as usize]
                [subject[(s - HSP_MAX_WINDOW) as usize] as usize];
            score += matrix[query[q as usize] as usize][subject[s as usize] as usize];
        }
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
        q += 1;
        s += 1;
    }
    if max_score > 0 {
        *q_retval = max_offset;
        *s_retval = max_offset - q_start + s_start;
        return true;
    }
    score = 0;
    q = q_start + q_length - HSP_MAX_WINDOW;
    s = s_start + s_length - HSP_MAX_WINDOW;
    for index in hsp.query.end - HSP_MAX_WINDOW..hsp.query.end {
        score += if let Some(pssm) = pssm {
            pssm[index as usize][subject[s as usize] as usize]
        } else {
            matrix[query[q as usize] as usize][subject[s as usize] as usize]
        };
        q += 1;
        s += 1;
    }
    if score > 0 {
        *q_retval = hsp.query.end - HSP_MAX_WINDOW / 2;
        *s_retval = hsp.subject.end - HSP_MAX_WINDOW / 2;
        return true;
    }
    false
}

fn nucleotide_at(sequence: &DynamicSequence<'_>, index: i32) -> u8 {
    sequence.data[(sequence.offset + index) as usize]
}
fn nucleotide_slice<'a>(sequence: &DynamicSequence<'a>) -> &'a [u8] {
    &sequence.data[sequence.offset as usize..]
}
fn nucleotide_shift<'a>(sequence: &DynamicSequence<'a>, offset: i32) -> DynamicSequence<'a> {
    DynamicSequence {
        data: sequence.data,
        offset: sequence.offset + offset,
    }
}

pub fn blast_get_start_for_gapped_alignment_nucl(
    query: &DynamicSequence<'_>,
    subject: &DynamicSequence<'_>,
    hsp: &mut BlastHSP,
) {
    let mut hsp_max_ident_run = 10;
    let offset = (hsp.subject.gapped_start - hsp.subject.offset)
        .min(hsp.query.gapped_start - hsp.query.offset);
    let mut q_start = hsp.query.gapped_start;
    let mut s_start = hsp.subject.gapped_start;
    let mut score = -1;
    let (mut q, mut s) = (q_start, s_start);
    let mut q_len = hsp.query.end;
    while q < q_len {
        let equal = nucleotide_at(query, q) == nucleotide_at(subject, s);
        q += 1;
        s += 1;
        if !equal {
            break;
        }
        score += 1;
        if score > hsp_max_ident_run {
            return;
        }
    }
    q = q_start;
    s = s_start;
    while q >= 0 {
        let equal = nucleotide_at(query, q) == nucleotide_at(subject, s);
        q -= 1;
        s -= 1;
        if !equal {
            break;
        }
        score += 1;
        if score > hsp_max_ident_run {
            return;
        }
    }
    hsp_max_ident_run = (hsp_max_ident_run as f64 * 1.5) as i32;
    q_start = hsp.query.gapped_start - offset;
    s_start = hsp.subject.gapped_start - offset;
    q_len = (hsp.subject.end - s_start).min(hsp.query.end - q_start);
    q = q_start;
    s = s_start;
    let (mut max_score, mut max_offset) = (0, q_start);
    score = 0;
    let (mut is_match, mut prev_match) = (false, false);
    let mut index = q_start;
    while index < q_start + q_len {
        is_match = nucleotide_at(query, q) == nucleotide_at(subject, s);
        q += 1;
        s += 1;
        if is_match != prev_match {
            prev_match = is_match;
            if is_match {
                score = 1;
            } else if score > max_score {
                max_score = score;
                max_offset = index - score / 2;
            }
        } else if is_match {
            score += 1;
            if score > hsp_max_ident_run {
                max_offset = index - hsp_max_ident_run / 2;
                hsp.query.gapped_start = max_offset;
                hsp.subject.gapped_start = max_offset + s_start - q_start;
                return;
            }
        }
        index += 1;
    }
    if is_match && score > max_score {
        max_score = score;
        max_offset = index - score / 2;
    }
    if max_score > 0 {
        hsp.query.gapped_start = max_offset;
        hsp.subject.gapped_start = max_offset + s_start - q_start;
    }
}

pub fn blast_hsp_adjust_subject_offset(hsp: &mut BlastHSP, start_shift: i32) {
    if start_shift > 0 {
        hsp.subject.offset += start_shift;
        hsp.subject.end += start_shift;
        hsp.subject.gapped_start += start_shift;
    }
}

/// Complete in-memory, ordinary blastn traceback. Both query views include
/// concatenated contexts and preserve their sentinel origin. Subject is a full
/// uncompressed sequence; partial fetching/fence rollback are not exposed.
pub fn blast_traceback_from_hsp_list(
    program: EBlastProgramType,
    list: &mut BlastHSPList,
    query_blk: &GapSequence<'_>,
    query_nomask: &DynamicSequence<'_>,
    subject_blk: &GapSubject<'_>,
    query_info: &BlastQueryInfo,
    gap_align: &mut BlastGapAlignStruct<'_>,
    sbp: &mut BlastScoreBlk,
    score_params: &BlastScoringParameters,
    ext_options: &BlastExtensionOptions,
    hit_params: &BlastHitSavingParameters,
) -> i16 {
    let num_initial = list.hspcnt as usize;
    if num_initial == 0 {
        return 0;
    }
    assert_eq!(program, EBlastProgramType::Blastn);
    assert!(!score_params.options.is_ooframe);
    assert!(
        !score_params.options.complexity_adjusted_scoring,
        "nondefault complexity adjustment is outside scope"
    );
    assert!(query_blk.oof_sequence.is_none());
    let greedy = ext_options.e_tback_ext == EBlastTbackExt::GreedyTbck;
    assert!(
        greedy || ext_options.e_tback_ext == EBlastTbackExt::DynProgTbck,
        "Smith-Waterman traceback is outside standard scope"
    );
    let hit_options = hit_params.options.as_ref().expect("hit options required");
    let subject = subject_blk.sequence;
    let subject_length = subject_blk.length;
    assert!(blast_hsp_list_is_sorted_by_score(Some(list)));
    let Some(mut tree) = blast_interval_tree_init(0, query_blk.length + 1, 0, subject_length + 1)
    else {
        return -1;
    };
    for index in 0..num_initial {
        let hsp = list.hsp_array[index].as_mut().expect("nonnull initial HSP");
        let ctx = &query_info.contexts[hsp.context as usize];
        let query = nucleotide_shift(&query_blk.sequence.unwrap(), ctx.query_offset);
        let nomask = nucleotide_shift(query_nomask, ctx.query_offset);
        let query_length = ctx.query_length;
        // SAFETY: only accepted, previously processed HSP boxes are indexed;
        // they stay live and have stable coordinates through this loop.
        if unsafe {
            blast_interval_tree_contains_hsp(
                &tree,
                hsp,
                query_info,
                hit_options.min_diag_separation,
            )
        } {
            list.hsp_array[index] = None;
            continue;
        }
        let (mut q_start, mut s_start) = (0, 0);
        if hsp.query.gapped_start == 0 && hsp.subject.gapped_start == 0 {
            if !blast_get_offsets_for_gapped_alignment(
                nucleotide_slice(&query),
                nucleotide_slice(&subject),
                &sbp.matrix,
                None,
                hsp,
                &mut q_start,
                &mut s_start,
            ) {
                list.hsp_array[index] = None;
                continue;
            }
            hsp.query.gapped_start = q_start;
            hsp.subject.gapped_start = s_start;
        } else {
            blast_get_start_for_gapped_alignment_nucl(&query, &subject, hsp);
            q_start = hsp.query.gapped_start;
            s_start = hsp.subject.gapped_start;
        }
        let mut adjusted_s_length = subject_length;
        let mut start_shift = 0;
        adjust_subject_range(
            &mut s_start,
            &mut adjusted_s_length,
            q_start,
            query_length,
            &mut start_shift,
        );
        let adjusted_subject = nucleotide_shift(&subject, start_shift);
        hsp.subject.gapped_start = s_start;
        if greedy {
            blast_greedy_gapped_alignment(
                nucleotide_slice(&query),
                nucleotide_slice(&adjusted_subject),
                query_length,
                adjusted_s_length,
                gap_align,
                score_params,
                q_start,
                s_start,
                false,
                true,
                None,
            );
        } else {
            blast_gapped_alignment_with_traceback(
                program,
                &query,
                &adjusted_subject,
                Some(gap_align),
                score_params,
                q_start,
                s_start,
                query_length,
                adjusted_s_length,
                None,
            );
        }
        blast_hsp_update_with_traceback(Some(gap_align), Some(hsp));
        let mut delete_hsp = false;
        if !greedy {
            let mut align_length = 0;
            blast_hsp_get_num_identities_and_positives(
                Some(nucleotide_slice(&nomask)),
                Some(nucleotide_slice(&adjusted_subject)),
                hsp,
                &score_params.options,
                Some(&mut align_length),
                Some(sbp),
            );
            delete_hsp = blast_hsp_test(hsp, hit_options, align_length);
        }
        if !delete_hsp {
            blast_hsp_adjust_subject_offset(hsp, start_shift);
            // SAFETY: boxed HSP stays allocated until the explicit tree reset
            // below; later loop iterations mutate different HSP owners only.
            let status = unsafe {
                blast_interval_tree_add_hsp(
                    NonNull::from(&mut **hsp),
                    &mut tree,
                    query_info,
                    EITreeIndexMethod::QueryAndSubject,
                )
            };
            if status != 0 {
                return status;
            }
        } else {
            list.hsp_array[index] = None;
        }
    }
    // End the raw-pointer borrow before clipping/freeing any indexed HSP.
    // Source does not consult this tree again until its later reset.
    blast_interval_tree_reset(&mut tree);
    blast_hsp_list_purge_null_hs_ps(Some(list));
    let mut extra_start =
        blast_hsp_list_purge_hs_ps_with_common_endpoints(program, Some(list), false);
    if greedy {
        extra_start = 0;
    }
    for index in extra_start as usize..list.hspcnt as usize {
        let Some(hsp) = list.hsp_array[index].as_mut() else {
            continue;
        };
        let ctx = &query_info.contexts[hsp.context as usize];
        let query = nucleotide_shift(&query_blk.sequence.unwrap(), ctx.query_offset);
        let nomask = nucleotide_shift(query_nomask, ctx.query_offset);
        let mut delete_hsp = blast_hsp_reevaluate_with_ambiguities_gapped(
            hsp,
            nucleotide_slice(&query),
            ctx.query_length,
            nucleotide_slice(&subject),
            subject_length,
            hit_params,
            score_params,
            sbp,
        );
        if !delete_hsp {
            delete_hsp = blast_hsp_test_identity_and_length(
                program,
                hsp,
                nucleotide_slice(&nomask),
                nucleotide_slice(&subject),
                &score_params.options,
                hit_options,
            );
        }
        if delete_hsp {
            list.hsp_array[index] = None;
        }
    }
    blast_hsp_list_purge_null_hs_ps(Some(list));
    blast_hsp_list_purge_hs_ps_with_common_endpoints(program, Some(list), true);
    blast_hsp_list_sort_by_score(Some(list));
    blast_interval_tree_reset(&mut tree);
    for index in 0..list.hspcnt as usize {
        let hsp = list.hsp_array[index].as_mut().unwrap();
        // SAFETY: accepted HSPs precede this index and remain live/unchanged;
        // contained HSPs were never inserted into this tree.
        if unsafe {
            blast_interval_tree_contains_hsp(
                &tree,
                hsp,
                query_info,
                hit_options.min_diag_separation,
            )
        } {
            list.hsp_array[index] = None;
        } else {
            let status = unsafe {
                blast_interval_tree_add_hsp(
                    NonNull::from(&mut **hsp),
                    &mut tree,
                    query_info,
                    EITreeIndexMethod::QueryAndSubject,
                )
            };
            if status != 0 {
                return status;
            }
        }
    }
    drop(tree);
    blast_hsp_list_purge_null_hs_ps(Some(list));
    s_hsp_list_post_traceback_update(
        program,
        list,
        query_info,
        score_params,
        hit_params,
        sbp,
        subject_length,
    );
    0
}
