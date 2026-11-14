// NCBI BLAST+ 2.17.0+, selected standard gapped blast_engine.c/hits.c closure.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_diagnostics::BlastGappedStats;
use super::blast_hits::*;
use super::blast_hsp_evaluate::{
    blast_hsp_list_get_evalues, blast_hsp_list_purge_hs_ps_with_common_endpoints,
};
use super::blast_hsp_link::blast_link_hsps;
use super::blast_parameters::{BlastHitSavingParameters, BlastScoringParameters};
use super::blast_program::{blast_subject_is_translated, EBlastProgramType};
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::BlastScoreBlk;

pub fn s_blast_hsp_list_reap_by_prelim_evalue(
    hsp_list: Option<&mut BlastHSPList>,
    hit_params: &BlastHitSavingParameters,
) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    let mut count = 0;
    for index in 0..list.hspcnt as usize {
        let hsp = list.hsp_array[index].as_ref().expect("nonnull active HSP");
        if hsp.evalue > hit_params.prelim_evalue {
            list.hsp_array[index] = None;
        } else {
            if index > count {
                list.hsp_array[count] = list.hsp_array[index].take();
            }
            count += 1;
        }
    }
    list.hspcnt = count as i32;
    0
}

pub fn s_blast_hsp_lists_combine_by_score(
    hsp_list: &mut BlastHSPList,
    combined: &mut BlastHSPList,
    new_hspcnt: i32,
) {
    assert!(new_hspcnt <= combined.allocated);
    if new_hspcnt >= hsp_list.hspcnt + combined.hspcnt {
        let mut index = combined.hspcnt as usize;
        for index1 in 0..hsp_list.hspcnt as usize {
            if hsp_list.hsp_array[index1].is_some() {
                combined.hsp_array[index] = hsp_list.hsp_array[index1].take();
                index += 1;
            }
        }
        combined.hspcnt = new_hspcnt;
        blast_hsp_list_sort_by_score(Some(combined));
    } else {
        let mut new_array: Vec<Option<Box<BlastHSP>>> =
            (0..combined.allocated).map(|_| None).collect();
        blast_hsp_list_sort_by_score(Some(combined));
        blast_hsp_list_sort_by_score(Some(hsp_list));
        let (mut index1, mut index2) = (0, 0);
        for index in 0..new_hspcnt as usize {
            if index1 < combined.hspcnt as usize
                && (index2 >= hsp_list.hspcnt as usize
                    || score_compare_hs_ps(
                        combined.hsp_array[index1].as_deref(),
                        hsp_list.hsp_array[index2].as_deref(),
                    )
                    .is_le())
            {
                new_array[index] = combined.hsp_array[index1].take();
                index1 += 1;
            } else {
                new_array[index] = hsp_list.hsp_array[index2].take();
                index2 += 1;
            }
        }
        for index in index1..combined.hspcnt as usize {
            combined.hsp_array[index] = None;
        }
        for index in index2..hsp_list.hspcnt as usize {
            hsp_list.hsp_array[index] = None;
        }
        combined.hsp_array = new_array;
        combined.hspcnt = new_hspcnt;
    }
    hsp_list.hspcnt = 0;
}

pub fn blast_hsp_list_append(
    old_hsp_list: &mut Option<Box<BlastHSPList>>,
    combined_hsp_list: &mut Option<Box<BlastHSPList>>,
    hsp_num_max: i32,
) -> i16 {
    let Some(old) = old_hsp_list.as_ref() else {
        return 0;
    };
    if old.hspcnt == 0 {
        return 0;
    }
    if combined_hsp_list.is_none() {
        *combined_hsp_list = old_hsp_list.take();
        return 0;
    }
    let combined = combined_hsp_list.as_mut().unwrap();
    let mut new_hspcnt = (combined.hspcnt + old.hspcnt).min(hsp_num_max);
    if new_hspcnt > combined.allocated && !combined.do_not_reallocate {
        let new_allocated = (2 * new_hspcnt).min(hsp_num_max);
        let additional = new_allocated as usize - combined.hsp_array.len();
        if combined.hsp_array.try_reserve_exact(additional).is_ok() {
            combined
                .hsp_array
                .resize_with(new_allocated as usize, || None);
            combined.allocated = new_allocated;
        } else {
            combined.do_not_reallocate = true;
            new_hspcnt = combined.allocated;
        }
    }
    if combined.allocated == hsp_num_max {
        combined.do_not_reallocate = true;
    }
    s_blast_hsp_lists_combine_by_score(old_hsp_list.as_mut().unwrap(), combined, new_hspcnt);
    *old_hsp_list = None;
    0
}

/// Unsplit complete subject context. The callback binds the source word finder
/// and BLAST_GetGappedScore calls; supplied result HSPs already use local query
/// contexts. Only the ordinary gapped branch is translated here.
pub fn s_blast_search_engine_one_context<F>(
    program: EBlastProgramType,
    subject_oid: i32,
    search_context: F,
    hsp_list_out: &mut Option<Box<BlastHSPList>>,
) -> i16
where
    F: FnOnce() -> Result<Option<Box<BlastHSPList>>, i16>,
{
    *hsp_list_out = None;
    let mut list = match search_context() {
        Ok(list) => list,
        Err(status) => return status,
    };
    blast_hsp_list_purge_hs_ps_with_common_endpoints(program, list.as_deref_mut(), true);
    // The source odd-blastn-score adjustment at this stage is #if 0.
    blast_hsp_list_sort_by_score(list.as_deref_mut());
    if list.as_ref().is_none_or(|list| list.hspcnt == 0) {
        return 0;
    }
    list.as_mut().unwrap().oid = subject_oid;
    // A full unsplit subject has backup.offset=0 and no combined list;
    // the source first-chunk merge is an ownership transfer.
    *hsp_list_out = list;
    0
}

/// Binding stage of the callback result: either the gapped-score output for
/// a complete unsplit context, or an already completed source OneContext call
/// (including per-chunk endpoint purges, coordinate shifts and overlap merges).
pub enum PreliminaryContextResult {
    GappedScore(Option<Box<BlastHSPList>>),
    Finished(Option<Box<BlastHSPList>>),
}
impl From<Option<Box<BlastHSPList>>> for PreliminaryContextResult {
    fn from(list: Option<Box<BlastHSPList>>) -> Self {
        Self::GappedScore(list)
    }
}

/// Source core frame loop and statistical/threshold tail. Complete translated
/// frame buffers (or the single blastn context) are prepared by the caller;
/// the callback binds each existing word-finder/gapped-score search.
/// For tblastn `original_subject_length` is nucleotide length and frame_lengths
/// lists amino-acid lengths in the original +1,+2,+3,-1,-2,-3 traversal order.
pub fn s_blast_search_engine_core<F, R>(
    program: EBlastProgramType,
    frame_lengths: &[i32],
    original_subject_length: i32,
    subject_oid: i32,
    query_info: &BlastQueryInfo,
    score_params: &BlastScoringParameters,
    hit_params: &BlastHitSavingParameters,
    sbp: &mut BlastScoreBlk,
    mut diagnostics: Option<&mut BlastGappedStats>,
    mut search_context: F,
    hsp_list_out: &mut Option<Box<BlastHSPList>>,
) -> i16
where
    F: FnMut(usize) -> Result<R, i16>,
    R: Into<PreliminaryContextResult>,
{
    *hsp_list_out = None;
    assert!(matches!(
        program,
        EBlastProgramType::Blastn | EBlastProgramType::Tblastn
    ));
    assert!(score_params.options.gapped_calculation);
    assert!(!score_params.options.is_ooframe);
    let translated = blast_subject_is_translated(program);
    if !translated {
        assert_eq!(frame_lengths.len(), 1);
    }
    let hit_options = hit_params.options.as_ref().expect("hit options required");
    let hsp_num_max = blast_hsp_num_max(true, hit_options.hsp_num_max);
    let mut combined = None;
    let mut stat_length = original_subject_length;
    for (context, &frame_length) in frame_lengths.iter().enumerate() {
        assert!(frame_length >= 0);
        if frame_length == 0 {
            continue;
        }
        if translated && frame_length > 0 {
            stat_length = frame_length;
        }
        let result = match search_context(context) {
            Ok(result) => result.into(),
            Err(status) => return status,
        };
        let mut context_list = match result {
            PreliminaryContextResult::GappedScore(list) => {
                assert!(
                    frame_length <= 5_000_000,
                    "split contexts require a completed OneContext result"
                );
                let mut output = None;
                let status = s_blast_search_engine_one_context(
                    program,
                    subject_oid,
                    || Ok(list),
                    &mut output,
                );
                if status != 0 {
                    return status;
                }
                output
            }
            PreliminaryContextResult::Finished(list) => list,
        };
        if blast_hsp_list_append(&mut context_list, &mut combined, hsp_num_max) != 0 {
            return 1;
        }
    }
    if let Some(link) = &hit_params.link_hsp_params {
        blast_link_hsps(
            program,
            combined.as_deref_mut(),
            query_info,
            original_subject_length,
            sbp,
            link,
            true,
        );
    } else {
        blast_hsp_list_get_evalues(
            program,
            query_info,
            stat_length,
            combined.as_deref_mut(),
            true,
            false,
            sbp,
            0.0,
            1.0,
        );
    }
    let status = if sbp.matrix_only_scoring {
        blast_hsp_list_reap_by_raw_score(combined.as_deref_mut(), hit_options.cutoff_score)
    } else {
        s_blast_hsp_list_reap_by_prelim_evalue(combined.as_deref_mut(), hit_params)
    };
    if combined.as_ref().is_some_and(|list| list.hspcnt == 0) {
        combined = None;
    }
    if let (Some(stats), Some(list)) = (diagnostics.as_mut(), combined.as_ref()) {
        if list.hspcnt > 0 {
            stats.num_seqs_passed += 1;
            stats.good_extensions += list.hspcnt;
        }
    }
    *hsp_list_out = combined;
    status
}
