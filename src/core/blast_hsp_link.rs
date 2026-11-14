// NCBI BLAST+ 2.17.0+, uneven-gap link_hsps.c/h.
// Author: Ilya Dondoshansky. See LICENSE-NCBI.
use super::blast_hits::*;
use super::blast_hsp_evaluate::blast_hsp_list_get_evalues;
use super::blast_parameters::BlastLinkHSPParameters;
use super::blast_program::{blast_subject_is_translated, EBlastProgramType};
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::{blast_gap_decay_divisor, BlastKarlinBlk, BlastScoreBlk};
use super::hsp_sumstat::blast_uneven_gap_sum_e;
use std::cmp::Ordering;

/// Stable arena indexes replace wrapper/list pointers. `hsp` indexes the
/// original HSP array, whose owners remain stationary throughout linking.
#[derive(Clone, Debug)]
pub struct BlastLinkedHSPSet {
    pub hsp: usize,
    pub query_id: u32,
    pub next: Option<usize>,
    pub prev: Option<usize>,
    pub sum_score: f64,
}

pub fn s_sum_hsp_evalue(
    program: EBlastProgramType,
    query_info: &BlastQueryInfo,
    subject_length: i32,
    params: &BlastLinkHSPParameters,
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    head: usize,
    new: usize,
    sum_score: &mut f64,
) -> f64 {
    assert_ne!(program, EBlastProgramType::Tblastx);
    let h = list.hsp_array[sets[head].hsp].as_ref().unwrap();
    let new_h = list.hsp_array[sets[new].hsp].as_ref().unwrap();
    let mut subject_eff_length = if blast_subject_is_translated(program) {
        subject_length / 3
    } else {
        subject_length
    };
    let num = (h.num + new_h.num) as i16;
    let context = &query_info.contexts[h.context as usize];
    let len_adj = context.length_adjustment;
    let query_eff_length = (context.query_length - len_adj).max(1);
    subject_eff_length = (subject_eff_length - len_adj).max(1);
    *sum_score = sets[new].sum_score + sets[head].sum_score;
    let query_window = params.overlap_size + params.gap_size + 1;
    let subject_window = params.overlap_size + params.longest_intron + 1;
    blast_uneven_gap_sum_e(
        query_window,
        subject_window,
        num,
        *sum_score,
        query_eff_length,
        subject_eff_length,
        context.eff_searchsp,
        blast_gap_decay_divisor(params.gap_decay_rate, num as u32),
    )
}

pub fn s_fwd_compare_linked_hsp_sets(
    a: &BlastLinkedHSPSet,
    b: &BlastLinkedHSPSet,
    list: &BlastHSPList,
) -> Ordering {
    if a.query_id != b.query_id {
        return a.query_id.cmp(&b.query_id);
    }
    let a = list.hsp_array[a.hsp].as_ref().unwrap();
    let b = list.hsp_array[b.hsp].as_ref().unwrap();
    a.query
        .offset
        .cmp(&b.query.offset)
        .then_with(|| a.subject.offset.cmp(&b.subject.offset))
}

pub fn s_sum_score_compare_linked_hsp_sets(
    a: Option<&BlastLinkedHSPSet>,
    b: Option<&BlastLinkedHSPSet>,
    list: &BlastHSPList,
) -> Ordering {
    match (a, b) {
        (None, None) => Ordering::Equal,
        (None, _) => Ordering::Greater,
        (_, None) => Ordering::Less,
        (Some(a), Some(b)) => {
            if a.sum_score < b.sum_score {
                return Ordering::Greater;
            }
            if a.sum_score > b.sum_score {
                return Ordering::Less;
            }
            score_compare_hs_ps(
                list.hsp_array[a.hsp].as_deref(),
                list.hsp_array[b.hsp].as_deref(),
            )
        }
    }
}

pub fn s_hsp_offset_binary_search(
    array: &[usize],
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    query_id: u32,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, array.len());
    while begin < end {
        let index = (begin + end) / 2;
        let set = &sets[array[index]];
        if set.query_id < query_id {
            begin = index + 1;
        } else if set.query_id > query_id {
            end = index;
        } else if list.hsp_array[set.hsp].as_ref().unwrap().query.offset >= offset {
            end = index;
        } else {
            begin = index + 1;
        }
    }
    end
}

pub fn s_hsp_offset_end_binary_search(
    array: &[usize],
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    qend_index_array: &[usize],
    query_id: u32,
    offset: i32,
) -> usize {
    let (mut begin, mut end) = (0, array.len());
    while begin < end {
        let right = (begin + end) / 2;
        let left = qend_index_array[right];
        if sets[array[right]].query_id < query_id {
            begin = right + 1;
        } else if sets[array[right]].query_id > query_id {
            end = left;
        } else if list.hsp_array[sets[array[left]].hsp]
            .as_ref()
            .unwrap()
            .query
            .end
            >= offset
        {
            end = left;
        } else {
            begin = right + 1;
        }
    }
    end
}

pub fn s_merge_linked_hsp_sets(
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    mut set1: usize,
    mut set2: usize,
    merged_size: &mut i32,
) -> Vec<usize> {
    while let Some(prev) = sets[set1].prev {
        set1 = prev;
    }
    while let Some(prev) = sets[set2].prev {
        set2 = prev;
    }
    *merged_size = list.hsp_array[sets[set1].hsp].as_ref().unwrap().num
        + list.hsp_array[sets[set2].hsp].as_ref().unwrap().num;
    let mut merged = Vec::with_capacity(*merged_size as usize);
    let (mut first, mut second) = (Some(set1), Some(set2));
    while first.is_some() || second.is_some() {
        if second.is_none()
            || first.is_some_and(|a| {
                list.hsp_array[sets[a].hsp].as_ref().unwrap().query.offset
                    < list.hsp_array[sets[second.unwrap()].hsp]
                        .as_ref()
                        .unwrap()
                        .query
                        .offset
            })
        {
            let a = first.unwrap();
            merged.push(a);
            first = sets[a].next;
        } else {
            let b = second.unwrap();
            merged.push(b);
            second = sets[b].next;
        }
    }
    merged
}

pub fn s_combine_linked_hsp_sets(
    sets: &mut [BlastLinkedHSPSet],
    list: &mut BlastHSPList,
    set1: Option<usize>,
    set2: Option<usize>,
    sum_score: f64,
    evalue: f64,
) -> Option<usize> {
    let (first, second) = match (set1, set2) {
        (None, _) => return set2,
        (_, None) => return set1,
        (Some(a), Some(b)) => (a, b),
    };
    let mut new_num = 0;
    let merged = s_merge_linked_hsp_sets(sets, list, first, second, &mut new_num);
    let head = merged[0];
    sets[head].prev = None;
    for index in 0..new_num as usize {
        let link = merged[index];
        if index < new_num as usize - 1 {
            let next = merged[index + 1];
            sets[link].next = Some(next);
            sets[next].prev = Some(link);
        } else {
            sets[link].next = None;
        }
        sets[link].sum_score = sum_score;
        let hsp = list.hsp_array[sets[link].hsp].as_mut().unwrap();
        hsp.evalue = evalue;
        hsp.num = new_num;
    }
    Some(head)
}

pub fn s_linked_hsp_sets_admissible(
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    set1: Option<usize>,
    set2: Option<usize>,
    params: Option<&BlastLinkHSPParameters>,
    program: EBlastProgramType,
) -> bool {
    let (Some(first), Some(mut second), Some(params)) = (set1, set2, params) else {
        return false;
    };
    if sets[first].prev.is_some() {
        return false;
    }
    while let Some(prev) = sets[second].prev {
        second = prev;
    }
    if first == second || sets[first].query_id != sets[second].query_id {
        return false;
    }
    if list.hsp_array[sets[first].hsp]
        .as_ref()
        .unwrap()
        .subject
        .frame
        .signum()
        != list.hsp_array[sets[second].hsp]
            .as_ref()
            .unwrap()
            .subject
            .frame
            .signum()
    {
        return false;
    }
    let mut combined_size = 0;
    let merged = s_merge_linked_hsp_sets(sets, list, first, second, &mut combined_size);
    let (mut gap_s, mut gap_q) = (params.longest_intron, params.gap_size);
    if program == EBlastProgramType::Blastx {
        (gap_s, gap_q) = (params.gap_size, params.longest_intron);
    }
    let overlap = params.overlap_size;
    for index in 0..combined_size as usize - 1 {
        let left = list.hsp_array[sets[merged[index]].hsp].as_ref().unwrap();
        let right = list.hsp_array[sets[merged[index + 1]].hsp]
            .as_ref()
            .unwrap();
        if left.query.end < right.query.offset - gap_q {
            return false;
        }
        if left.query.offset >= right.query.offset {
            return false;
        }
        if left.query.end > right.query.offset + overlap || left.query.end >= right.query.end {
            return false;
        }
        if left.subject.end > right.subject.offset + overlap
            || left.subject.end < right.subject.offset - gap_s
            || left.subject.offset >= right.subject.offset
            || left.subject.end >= right.subject.end
        {
            return false;
        }
    }
    true
}

pub fn s_linked_hsp_set_array_set_up(
    list: &mut BlastHSPList,
    kbp_array: &[Option<BlastKarlinBlk>],
    program: EBlastProgramType,
) -> Vec<BlastLinkedHSPSet> {
    let mut array = Vec::with_capacity(list.hspcnt as usize);
    for index in 0..list.hspcnt as usize {
        let hsp = list.hsp_array[index].as_mut().unwrap();
        let kbp = kbp_array[hsp.context as usize].as_ref().unwrap();
        array.push(BlastLinkedHSPSet {
            hsp: index,
            sum_score: kbp.lambda * hsp.score as f64 - kbp.log_k,
            query_id: if program == EBlastProgramType::Blastx {
                (hsp.context / 3) as u32
            } else {
                hsp.context as u32
            },
            next: None,
            prev: None,
        });
        hsp.num = 1;
    }
    array
}

pub fn s_linked_hsp_set_array_clean_up(array: Vec<BlastLinkedHSPSet>) {
    drop(array);
}

pub fn s_linked_hsp_set_array_index_query_ends(
    array: &[usize],
    sets: &[BlastLinkedHSPSet],
    list: &BlastHSPList,
    qend_index_ptr: &mut Vec<usize>,
) -> i16 {
    let mut indices = Vec::new();
    if indices.try_reserve_exact(array.len()).is_err() {
        return -1;
    }
    indices.resize(array.len(), 0);
    let mut current_end = list.hsp_array[sets[array[0]].hsp]
        .as_ref()
        .unwrap()
        .query
        .end;
    let mut current_index = 0;
    for index in 1..array.len() {
        let link = &sets[array[index]];
        let hsp = list.hsp_array[link.hsp].as_ref().unwrap();
        if link.query_id > sets[array[current_index]].query_id || hsp.query.end > current_end {
            current_index = index;
            current_end = hsp.query.end;
        }
        indices[index] = current_index;
    }
    *qend_index_ptr = indices;
    0
}

pub fn s_blast_uneven_gap_link_hs_ps(
    program: EBlastProgramType,
    hsp_list: Option<&mut BlastHSPList>,
    query_info: Option<&BlastQueryInfo>,
    subject_length: i32,
    sbp: Option<&BlastScoreBlk>,
    params: Option<&BlastLinkHSPParameters>,
    gapped_calculation: bool,
) -> i16 {
    let (Some(params), Some(sbp), Some(query_info)) = (params, sbp, query_info) else {
        return -1;
    };
    let Some(list) = hsp_list else {
        return 0;
    };
    if list.hspcnt <= 1 {
        return 0;
    }
    let kbp_array = if gapped_calculation {
        sbp.kbp_gap.as_deref().unwrap()
    } else {
        if sbp.kbp_is_psi {
            &sbp.kbp_psi
        } else {
            &sbp.kbp_std
        }
    };
    let gap_size = if program == EBlastProgramType::Blastx {
        params.longest_intron
    } else {
        params.gap_size
    };
    let hspcnt = list.hspcnt as usize;
    let mut sets = s_linked_hsp_set_array_set_up(list, kbp_array, program);
    let mut score_array: Vec<usize> = (0..hspcnt).collect();
    score_array.sort_by(|a, b| {
        s_sum_score_compare_linked_hsp_sets(Some(&sets[*a]), Some(&sets[*b]), list)
    });
    let mut offset_array: Vec<usize> = (0..hspcnt).collect();
    offset_array.sort_by(|a, b| s_fwd_compare_linked_hsp_sets(&sets[*a], &sets[*b], list));
    let mut qends = Vec::new();
    let status = s_linked_hsp_set_array_index_query_ends(&offset_array, &sets, list, &mut qends);
    if status != 0 {
        return status;
    }
    let mut head = None;
    let mut index = 0;
    while index < hspcnt {
        let mut best_hsp = None;
        if head.is_none() {
            while index < hspcnt
                && (sets[score_array[index]].next.is_some()
                    || sets[score_array[index]].prev.is_some())
            {
                index += 1;
            }
            if index == hspcnt {
                break;
            }
            head = Some(score_array[index]);
        }
        let current = head.unwrap();
        let mut tail = current;
        while let Some(next) = sets[tail].next {
            tail = next;
        }
        let hsp = list.hsp_array[sets[current].hsp].as_ref().unwrap();
        let mut best_evalue = hsp.evalue;
        let mut best_sum_score = sets[current].sum_score;
        let left_offset = hsp.query.offset - gap_size;
        let left = s_hsp_offset_end_binary_search(
            &offset_array,
            &sets,
            list,
            &qends,
            sets[current].query_id,
            left_offset,
        );
        let right = s_hsp_offset_binary_search(
            &offset_array,
            &sets,
            list,
            sets[tail].query_id,
            list.hsp_array[sets[tail].hsp].as_ref().unwrap().query.end + gap_size,
        );
        for index1 in left..right {
            let candidate = offset_array[index1];
            if sets[candidate].prev.is_some_and(|prev| {
                list.hsp_array[sets[prev].hsp].as_ref().unwrap().query.end >= left_offset
            }) {
                continue;
            }
            if s_linked_hsp_sets_admissible(
                &sets,
                list,
                head,
                Some(candidate),
                Some(params),
                program,
            ) {
                let mut sum_score = 0.0;
                let evalue = s_sum_hsp_evalue(
                    program,
                    query_info,
                    subject_length,
                    params,
                    &sets,
                    list,
                    current,
                    candidate,
                    &mut sum_score,
                );
                let candidate_evalue = list.hsp_array[sets[candidate].hsp].as_ref().unwrap().evalue;
                // C MIN comparison, including its asymmetric NaN choice.
                let minimum = if best_evalue < candidate_evalue {
                    best_evalue
                } else {
                    candidate_evalue
                };
                if evalue < minimum {
                    best_hsp = Some(candidate);
                    best_evalue = evalue;
                    best_sum_score = sum_score;
                }
            }
        }
        if best_hsp.is_some() {
            head = s_combine_linked_hsp_sets(
                &mut sets,
                list,
                head,
                best_hsp,
                best_sum_score,
                best_evalue,
            );
        } else {
            head = None;
            index += 1;
        }
    }
    s_linked_hsp_set_array_clean_up(sets);
    0
}

/// Default translated gapped searches use positive longest_intron and the
/// uneven-gap method. The historical even-gap engine is outside this closure.
pub fn blast_link_hsps(
    program: EBlastProgramType,
    hsp_list: Option<&mut BlastHSPList>,
    query_info: &BlastQueryInfo,
    subject_length: i32,
    sbp: &mut BlastScoreBlk,
    params: &BlastLinkHSPParameters,
    gapped_calculation: bool,
) -> i16 {
    let Some(list) = hsp_list else {
        return 0;
    };
    if list.hspcnt == 0 {
        return 0;
    }
    assert!(
        params.longest_intron > 0,
        "even-gap linking is outside standard tblastn scope"
    );
    assert!(
        matches!(
            program,
            EBlastProgramType::Tblastn | EBlastProgramType::Blastx
        ),
        "ordinary translated linking required"
    );
    for hsp in &mut list.hsp_array[..list.hspcnt as usize] {
        hsp.as_mut().unwrap().num = 1;
    }
    blast_hsp_list_adjust_odd_blastn_scores(Some(list), gapped_calculation, sbp.round_down);
    blast_hsp_list_get_evalues(
        program,
        query_info,
        if blast_subject_is_translated(program) {
            subject_length / 3
        } else {
            subject_length
        },
        Some(list),
        gapped_calculation,
        false,
        sbp,
        params.gap_decay_rate,
        1.0,
    );
    s_blast_uneven_gap_link_hs_ps(
        program,
        Some(list),
        Some(query_info),
        subject_length,
        Some(sbp),
        Some(params),
        gapped_calculation,
    );
    blast_hsp_list_sort_by_score(Some(list));
    list.best_evalue = list.hsp_array[0].as_ref().unwrap().evalue;
    for index in 1..list.hspcnt as usize {
        let value = list.hsp_array[index].as_ref().unwrap().evalue;
        if value < list.best_evalue {
            list.best_evalue = value;
        }
    }
    0
}
