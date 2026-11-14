// NCBI BLAST+ 2.17.0+, selected blast_gapalign.c preliminary orchestration.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_diagnostics::BlastGappedStats;
use super::blast_dynamic_align::{s_blast_dyn_prog_nt_gapped_alignment, DynamicSequence};
use super::blast_extend::{blast_init_hit_list_is_sorted_by_score, BlastInitHSP, BlastInitHitList};
use super::blast_gapalign::{blast_greedy_gapped_alignment, BlastGapAlignStruct};
use super::blast_hits::*;
use super::blast_options::EBlastPrelimGapExt;
use super::blast_parameters::{
    BlastExtensionParameters, BlastHitSavingParameters, BlastInitialWordParameters,
    BlastScoringParameters,
};
use super::blast_program::EBlastProgramType;
use super::blast_protein_gapalign::s_blast_prot_gapped_alignment;
use super::blast_query_info::{
    b_search_context_info, blast_get_query_index_from_query_offset, BlastQueryInfo,
};
use super::interval_tree::*;
use std::ptr::NonNull;

pub const HSP_MAX_WINDOW: i32 = 11;
/// Only the source sequence-block fields consumed by query localization.
#[derive(Clone, Copy, Debug)]
pub struct GapSequence<'a> {
    pub sequence: Option<DynamicSequence<'a>>,
    pub oof_sequence: Option<DynamicSequence<'a>>,
    pub length: i32,
}
#[derive(Clone, Copy, Debug)]
pub struct GapSubject<'a> {
    pub sequence: DynamicSequence<'a>,
    pub length: i32,
    pub frame: i16,
}

pub fn s_get_ungapped_hsp_context(query_info: &BlastQueryInfo, init_hsp: &BlastInitHSP) -> i32 {
    b_search_context_info(init_hsp.offsets.q_off, query_info)
}
pub fn s_adjust_initial_hsp_offsets(init_hsp: &mut BlastInitHSP, query_start: i32) {
    init_hsp.offsets.q_off -= query_start;
    if let Some(data) = &mut init_hsp.ungapped_data {
        data.q_start -= query_start;
    }
    assert!(init_hsp.ungapped_data.is_none_or(|data| data.q_start >= 0));
}
pub fn s_set_up_local_blast_sequence_blk<'a>(
    concatenated_query: &GapSequence<'a>,
    query_info: &BlastQueryInfo,
    context: i32,
    single_query: &mut GapSequence<'a>,
    query_start: &mut i32,
) {
    let query_length;
    if let Some(oof) = concatenated_query.oof_sequence {
        let mixed_frame_context = context - context % 3;
        *query_start = query_info.contexts[mixed_frame_context as usize].query_offset;
        let end = &query_info.contexts[(mixed_frame_context + 2) as usize];
        query_length = end.query_offset + end.query_length - *query_start;
        single_query.sequence = None;
        single_query.oof_sequence = Some(oof.shifted(*query_start));
    } else {
        let ctx = &query_info.contexts[context as usize];
        *query_start = ctx.query_offset;
        query_length = ctx.query_length;
        single_query.sequence = Some(concatenated_query.sequence.unwrap().shifted(*query_start));
        single_query.oof_sequence = None;
    }
    single_query.length = query_length;
}
pub fn s_adjust_hsp_offsets_and_get_query_data<'a>(
    query: &GapSequence<'a>,
    query_info: &BlastQueryInfo,
    init_hsp: &mut BlastInitHSP,
    query_out: &mut GapSequence<'a>,
    context: &mut i32,
) {
    assert!(init_hsp.ungapped_data.is_none_or(|data| data.q_start >= 0));
    *context = s_get_ungapped_hsp_context(query_info, init_hsp);
    let mut query_start = 0;
    s_set_up_local_blast_sequence_blk(query, query_info, *context, query_out, &mut query_start);
    assert!(init_hsp
        .ungapped_data
        .is_none_or(|data| data.q_start - query_start >= 0));
    s_adjust_initial_hsp_offsets(init_hsp, query_start);
}
pub fn blast_get_start_for_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    matrix: &[Vec<i32>],
    pssm: Option<&[Vec<i32>]>,
    q_start: u32,
    q_length: u32,
    s_start: u32,
    s_length: u32,
) -> i32 {
    if q_length <= HSP_MAX_WINDOW as u32 {
        return q_start.wrapping_add(q_length / 2) as i32;
    }
    let mut hsp_end = q_start.wrapping_add(HSP_MAX_WINDOW as u32) as i32;
    let mut q = q_start as usize;
    let mut s = s_start as usize;
    let mut score = 0i32;
    for index in q_start as i32..hsp_end {
        score = score.wrapping_add(if let Some(pssm) = pssm {
            pssm[index as usize][subject[s] as usize]
        } else {
            matrix[query[q] as usize][subject[s] as usize]
        });
        q += 1;
        s += 1;
    }
    let mut max_score = score;
    let mut max_offset = hsp_end - 1;
    hsp_end = q_start.wrapping_add(q_length.min(s_length)) as i32;
    for index in q_start.wrapping_add(HSP_MAX_WINDOW as u32) as i32..hsp_end {
        let w = HSP_MAX_WINDOW as usize;
        if let Some(pssm) = pssm {
            score = score
                .wrapping_sub(pssm[(index - HSP_MAX_WINDOW) as usize][subject[s - w] as usize]);
            score = score.wrapping_add(pssm[index as usize][subject[s] as usize]);
        } else {
            score = score.wrapping_sub(matrix[query[q - w] as usize][subject[s - w] as usize]);
            score = score.wrapping_add(matrix[query[q] as usize][subject[s] as usize]);
        }
        if score > max_score {
            max_score = score;
            max_offset = index;
        }
        q += 1;
        s += 1;
    }
    if max_score > 0 {
        max_offset -= HSP_MAX_WINDOW / 2;
    } else {
        max_offset = q_start as i32;
    }
    max_offset
}

/// Standard blastn/tblastn paths, excluding RPS, OOF, chaining and restricted
/// approximate protein alignment. Preliminary lists use the default unbounded
/// per-subject HSP count; boxed HSPs stay live until the interval tree is freed.
pub fn blast_get_gapped_score(
    program: EBlastProgramType,
    query: &GapSequence<'_>,
    query_info: &BlastQueryInfo,
    subject: &GapSubject<'_>,
    gap_align: &mut BlastGapAlignStruct<'_>,
    score_params: &BlastScoringParameters,
    ext_params: &BlastExtensionParameters,
    hit_params: &BlastHitSavingParameters,
    _word_params: &BlastInitialWordParameters,
    init_hitlist: &BlastInitHitList,
    hsp_list_ptr: Option<&mut Option<Box<BlastHSPList>>>,
    mut gapped_stats: Option<&mut BlastGappedStats>,
    mut fence_hit: Option<&mut bool>,
) -> i16 {
    let Some(hsp_list_ptr) = hsp_list_ptr else {
        return 1;
    };
    if init_hitlist.total == 0 {
        return 0;
    }
    assert!(matches!(
        program,
        EBlastProgramType::Blastn | EBlastProgramType::Tblastn
    ));
    assert!(
        !score_params.options.is_ooframe
            && !ext_params.options.chaining
            && !hit_params.restricted_align
    );
    let hit_options = hit_params.options.as_ref().unwrap();
    let max_hsps = blast_hsp_num_max(true, hit_options.hsp_num_max);
    assert_eq!(
        max_hsps,
        i32::MAX,
        "bounded preliminary HSP storage is outside the default path"
    );
    let is_prot = program == EBlastProgramType::Tblastn;
    let is_greedy = ext_params.options.e_prelim_gap_ext == EBlastPrelimGapExt::GreedyScoreOnly;
    assert!(blast_init_hit_list_is_sorted_by_score(init_hitlist));
    if hsp_list_ptr.is_none() {
        *hsp_list_ptr = Some(blast_hsp_list_new(max_hsps));
    }
    let list = hsp_list_ptr.as_deref_mut().unwrap();
    assert_eq!(list.hsp_max, i32::MAX);
    assert!(!list.do_not_reallocate);
    let Some(mut tree) = blast_interval_tree_init(0, query.length + 1, 0, subject.length + 1)
    else {
        return -1;
    };
    let mut found_high_score = vec![0; query_info.num_queries as usize];
    if let Some(low) = &hit_params.low_score {
        for hsp in &init_hitlist.init_hsp_array[..init_hitlist.total as usize] {
            let index =
                blast_get_query_index_from_query_offset(hsp.offsets.q_off, program, query_info)
                    as usize;
            if hsp.ungapped_data.as_ref().unwrap().score > low[index] {
                found_high_score[index] = 1;
            }
        }
    }
    let mut status = 0;
    for hsp in &init_hitlist.init_hsp_array[..init_hitlist.total as usize] {
        let mut init_hsp = BlastInitHSP {
            offsets: hsp.offsets,
            ungapped_data: hsp.ungapped_data,
        };
        let mut local_query = *query;
        let mut context = 0;
        s_adjust_hsp_offsets_and_get_query_data(
            query,
            query_info,
            &mut init_hsp,
            &mut local_query,
            &mut context,
        );
        if hit_params.low_score.is_some() {
            let index =
                super::blast_query_info::blast_get_query_index_from_context(context, program)
                    as usize;
            if found_high_score[index] == 0 {
                continue;
            }
        }
        let (q_start, q_end, s_start, s_end, score) = if let Some(data) = init_hsp.ungapped_data {
            (
                data.q_start,
                data.q_start + data.length,
                data.s_start,
                data.s_start + data.length,
                data.score,
            )
        } else {
            (
                init_hsp.offsets.q_off,
                init_hsp.offsets.q_off,
                init_hsp.offsets.s_off,
                init_hsp.offsets.s_off,
                i32::MIN,
            )
        };
        let mut tmp = BlastHSP {
            score,
            context,
            ..Default::default()
        };
        tmp.query = BlastSeg {
            offset: q_start,
            end: q_end,
            frame: query_info.contexts[context as usize].frame as i16,
            ..Default::default()
        };
        tmp.subject = BlastSeg {
            offset: s_start,
            end: s_end,
            frame: subject.frame,
            ..Default::default()
        };
        // Indexed boxes are retained by the unbounded list until tree drop.
        if unsafe {
            blast_interval_tree_contains_hsp(
                &tree,
                &tmp,
                query_info,
                hit_options.min_diag_separation,
            )
        } {
            continue;
        }
        let cutoff = hit_params.cutoffs[context as usize].cutoff_score;
        if let Some(stats) = gapped_stats.as_deref_mut() {
            stats.extensions = stats.extensions.wrapping_add(1);
        }
        let q = local_query.sequence.unwrap();
        if is_prot {
            let data = init_hsp.ungapped_data.as_ref().unwrap();
            let start = blast_get_start_for_gapped_alignment(
                &q.data[q.offset as usize..],
                &subject.sequence.data[subject.sequence.offset as usize..],
                gap_align.dynamic.matrix,
                gap_align.dynamic.pssm,
                data.q_start as u32,
                data.length as u32,
                data.s_start as u32,
                data.length as u32,
            );
            init_hsp.offsets.s_off += start - init_hsp.offsets.q_off;
            init_hsp.offsets.q_off = start;
            status = s_blast_prot_gapped_alignment(
                program,
                &q,
                local_query.length,
                &subject.sequence,
                subject.length,
                gap_align,
                score_params,
                &init_hsp,
                false,
                fence_hit.as_deref_mut(),
            );
        } else if is_greedy {
            if let Some(data) = init_hsp.ungapped_data {
                init_hsp.offsets.q_off = data.q_start + data.length / 2;
                init_hsp.offsets.s_off = data.s_start + data.length / 2;
            }
            status = blast_greedy_gapped_alignment(
                &q.data[q.offset as usize..],
                &subject.sequence.data[subject.sequence.offset as usize..],
                local_query.length,
                subject.length,
                gap_align,
                score_params,
                init_hsp.offsets.q_off,
                init_hsp.offsets.s_off,
                true,
                false,
                fence_hit.as_deref_mut(),
            );
            init_hsp.offsets.q_off = gap_align.greedy_query_seed_start;
            init_hsp.offsets.s_off = gap_align.greedy_subject_seed_start;
        } else {
            if s_end >= init_hsp.offsets.s_off + 8 {
                init_hsp.offsets.s_off += 3;
                init_hsp.offsets.q_off += 3;
            }
            status = s_blast_dyn_prog_nt_gapped_alignment(
                &q,
                local_query.length,
                &subject.sequence,
                subject.length,
                &mut gap_align.dynamic,
                score_params,
                &init_hsp,
            );
        }
        if status != 0 {
            return status;
        }
        if gap_align.dynamic.score >= cutoff {
            let mut new_hsp = None;
            status = blast_hsp_init(
                gap_align.dynamic.query_start,
                gap_align.dynamic.query_stop,
                gap_align.dynamic.subject_start,
                gap_align.dynamic.subject_stop,
                init_hsp.offsets.q_off,
                init_hsp.offsets.s_off,
                context,
                query_info.contexts[context as usize].frame as i16,
                subject.frame,
                gap_align.dynamic.score,
                Some(&mut gap_align.edit_script),
                Some(&mut new_hsp),
            );
            if status != 0 {
                return status;
            }
            let previous = list.hspcnt;
            status = blast_hsp_list_save_hsp(list, new_hsp);
            if status != 0 {
                break;
            }
            assert_eq!(
                list.hspcnt,
                previous + 1,
                "default preliminary list must retain each accepted box"
            );
            let pointer = NonNull::from(list.hsp_array[previous as usize].as_deref_mut().unwrap());
            // The list owns this box, and no indexed HSP is freed or moved.
            status = unsafe {
                blast_interval_tree_add_hsp(
                    pointer,
                    &mut tree,
                    query_info,
                    EITreeIndexMethod::QueryAndSubject,
                )
            };
            if status != 0 {
                break;
            }
        } else if is_greedy {
            gap_align.edit_script = None;
        }
    }
    status
}
