// NCBI BLAST+ 2.17.0+, selected blast_hits.c/gapinfo.c HSP evaluation.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_hits::*;
use super::blast_options::{BlastHitSavingOptions, BlastScoringOptions};
use super::blast_parameters::{
    BlastHitSavingParameters, BlastInitialWordParameters, BlastScoringParameters,
};
use super::blast_program::{
    blast_program_is_phi_blast, blast_program_is_rps_blast, EBlastProgramType,
};
use super::blast_query_info::BlastQueryInfo;
use super::blast_stat::{blast_gap_decay_divisor, blast_karlin_sto_e_simple, BlastScoreBlk};
use super::blast_stat_spouge::blast_spouge_sto_e;
use super::gapinfo::{gap_edit_script_new, EGapAlignOpType, GapEditScript};
use std::cmp::Ordering;

pub fn gap_edit_script_partial_copy(
    new: Option<&mut GapEditScript>,
    offset: i32,
    old: Option<&GapEditScript>,
    start: i32,
    stop: i32,
) -> i32 {
    let (Some(new), Some(old)) = (new, old) else {
        return -1;
    };
    let size = 1 + stop - start;
    if size > new.size {
        return -1;
    }
    let mut old_index = start;
    for new_index in offset..size + offset {
        new.num[new_index as usize] = old.num[old_index as usize];
        new.op_type[new_index as usize] = old.op_type[old_index as usize];
        old_index += 1;
    }
    0
}
/// Sequence positions are offsets from the original query/subject base pointers.
pub fn s_update_reevaluated_hsp(
    hsp: &mut BlastHSP,
    gapped: bool,
    cutoff_score: i32,
    score: i32,
    best_q_start: i32,
    best_q_end: i32,
    best_s_start: i32,
    best_s_end: i32,
    best_start_esp_index: i32,
    best_end_esp_index: i32,
    best_end_esp_num: i32,
) -> bool {
    hsp.score = score;
    if hsp.score >= cutoff_score {
        hsp.query.offset = best_q_start;
        hsp.query.end = best_q_end;
        hsp.subject.offset = best_s_start;
        hsp.subject.end = best_s_end;
        if gapped {
            let esp = hsp.gap_info.as_ref().expect("gapped edit script required");
            let last = esp.size - 1;
            if best_end_esp_index != last || best_start_esp_index > 0 {
                let mut new = gap_edit_script_new(best_end_esp_index - best_start_esp_index + 1)
                    .expect("nonempty selected script required");
                gap_edit_script_partial_copy(
                    Some(&mut new),
                    0,
                    Some(esp),
                    best_start_esp_index,
                    best_end_esp_index,
                );
                hsp.gap_info = Some(Box::new(new));
            }
            let esp = hsp.gap_info.as_mut().unwrap();
            let last = esp.size - 1;
            esp.num[last as usize] = best_end_esp_num;
            assert!(best_end_esp_num >= 0);
        }
        false
    } else {
        true
    }
}
pub fn blast_hsp_reevaluate_with_ambiguities_gapped(
    hsp: &mut BlastHSP,
    q: &[u8],
    qlen: i32,
    s: &[u8],
    slen: i32,
    hit_params: &BlastHitSavingParameters,
    score_params: &BlastScoringParameters,
    sbp: &BlastScoreBlk,
) -> bool {
    let cutoff_score = hit_params.cutoffs[hsp.context as usize].cutoff_score;
    let (mut factor, mut gap_open, mut gap_extend) =
        (1, score_params.gap_open, score_params.gap_extend);
    if gap_open == 0 && gap_extend == 0 {
        if score_params.reward % 2 == 1 {
            factor = 2;
        }
        gap_open = 0;
        gap_extend = (score_params.reward - 2 * score_params.penalty) * factor / 2;
    }
    let mut query = hsp.query.offset;
    let mut subject = hsp.subject.offset;
    let mut score = 0;
    let mut sum = 0;
    let mut best_q_start = query;
    let mut best_q_end = query;
    let mut current_q_start = query;
    let mut best_s_start = subject;
    let mut best_s_end = subject;
    let mut current_s_start = subject;
    let mut best_start_esp_index = 0;
    let mut best_end_esp_index = 0;
    let mut current_start_esp_index = 0;
    let mut best_end_esp_num = -1;
    let Some(esp) = hsp.gap_info.as_mut() else {
        return true;
    };
    for index in 0..esp.size as usize {
        let mut op_index = 0;
        while op_index < esp.num[index] {
            match esp.op_type[index] {
                EGapAlignOpType::Sub => {
                    sum += factor
                        * sbp.matrix[(q[query as usize] & 0x0f) as usize]
                            [s[subject as usize] as usize];
                    query += 1;
                    subject += 1;
                    op_index += 1;
                }
                EGapAlignOpType::Del => {
                    sum -= gap_open + gap_extend * esp.num[index];
                    subject += esp.num[index];
                    op_index += esp.num[index];
                }
                EGapAlignOpType::Ins => {
                    sum -= gap_open + gap_extend * esp.num[index];
                    query += esp.num[index];
                    op_index += esp.num[index];
                }
                _ => panic!("only ordinary substitution/indel scripts are supported"),
            }
            if sum < 0 {
                if op_index < esp.num[index] {
                    esp.num[index] -= op_index;
                    current_start_esp_index = index as i32;
                    op_index = 0;
                } else {
                    current_start_esp_index = index as i32 + 1;
                }
                sum = 0;
                current_q_start = query;
                current_s_start = subject;
                if score < cutoff_score {
                    best_q_start = query;
                    best_s_start = subject;
                    score = 0;
                    best_start_esp_index = current_start_esp_index;
                    best_end_esp_index = current_start_esp_index;
                }
            } else if sum > score {
                score = sum;
                best_q_start = current_q_start;
                best_s_start = current_s_start;
                best_q_end = query;
                best_s_end = subject;
                best_start_esp_index = current_start_esp_index;
                best_end_esp_index = index as i32;
                best_end_esp_num = op_index;
            }
        }
    }
    score /= factor;
    if best_start_esp_index < esp.size && best_end_esp_index < esp.size {
        assert_eq!(
            esp.op_type[best_start_esp_index as usize],
            EGapAlignOpType::Sub
        );
        assert_eq!(
            esp.op_type[best_end_esp_index as usize],
            EGapAlignOpType::Sub
        );
        let mut qp = best_q_start;
        let mut sp = best_s_start;
        let mut ext = 0;
        while qp > 0 && sp > 0 {
            qp -= 1;
            sp -= 1;
            if q[qp as usize] != s[sp as usize] || q[qp as usize] >= 4 {
                break;
            }
            ext += 1;
        }
        best_q_start -= ext;
        best_s_start -= ext;
        esp.num[best_start_esp_index as usize] += ext;
        if best_end_esp_index == best_start_esp_index {
            best_end_esp_num += ext;
        }
        score += ext * score_params.reward;
        qp = best_q_end;
        sp = best_s_end;
        ext = 0;
        while qp < qlen && sp < slen && q[qp as usize] < 4 {
            let equal = q[qp as usize] == s[sp as usize];
            qp += 1;
            sp += 1;
            if !equal {
                break;
            }
            ext += 1;
        }
        best_q_end += ext;
        best_s_end += ext;
        esp.num[best_end_esp_index as usize] += ext;
        best_end_esp_num += ext;
        score += ext * score_params.reward;
    }
    s_update_reevaluated_hsp(
        hsp,
        true,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
        best_start_esp_index,
        best_end_esp_index,
        best_end_esp_num,
    )
}
pub fn s_update_reevaluated_hsp_ungapped(
    hsp: &mut BlastHSP,
    cutoff_score: i32,
    score: i32,
    best_q_start: i32,
    best_q_end: i32,
    best_s_start: i32,
    best_s_end: i32,
) -> bool {
    s_update_reevaluated_hsp(
        hsp,
        false,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
        0,
        0,
        0,
    )
}
pub fn blast_hsp_reevaluate_with_ambiguities_ungapped(
    hsp: &mut BlastHSP,
    query_start: &[u8],
    subject_start: &[u8],
    word_params: &BlastInitialWordParameters,
    sbp: &BlastScoreBlk,
    translated: bool,
) -> bool {
    let mask = if translated { 0xff } else { 0x0f };
    let hsp_length = hsp.query.end - hsp.query.offset;
    let cutoff_score = word_params.cutoffs[hsp.context as usize].cutoff_score;
    let mut query = hsp.query.offset;
    let mut subject = hsp.subject.offset;
    let mut score = 0;
    let mut sum = 0;
    let mut best_q_start = query;
    let mut best_q_end = query;
    let mut current_q_start = query;
    let mut best_s_start = subject;
    let mut best_s_end = subject;
    let mut current_s_start = subject;
    for _ in 0..hsp_length {
        sum += sbp.matrix[(query_start[query as usize] & mask) as usize]
            [subject_start[subject as usize] as usize];
        query += 1;
        subject += 1;
        if sum < 0 {
            sum = 0;
            current_q_start = query;
            current_s_start = subject;
            if score < cutoff_score {
                best_q_start = query;
                best_q_end = query;
                best_s_start = subject;
                best_s_end = subject;
                score = 0;
            }
        } else if sum > score {
            score = sum;
            best_q_end = query;
            best_s_end = subject;
            best_q_start = current_q_start;
            best_s_start = current_s_start;
        }
    }
    s_update_reevaluated_hsp_ungapped(
        hsp,
        cutoff_score,
        score,
        best_q_start,
        best_q_end,
        best_s_start,
        best_s_end,
    )
}

pub fn s_blast_hsp_get_num_identities_and_positives(
    query: Option<&[u8]>,
    subject: Option<&[u8]>,
    hsp: &BlastHSP,
    num_ident_ptr: &mut i32,
    align_length_ptr: Option<&mut i32>,
    sbp: Option<&BlastScoreBlk>,
    num_pos_ptr: Option<&mut i32>,
) -> i16 {
    let q_length = hsp.query.end - hsp.query.offset;
    let s_length = hsp.subject.end - hsp.subject.offset;
    let (Some(query), Some(subject)) = (query, subject) else {
        return -1;
    };
    let mut q = hsp.query.offset;
    let mut s = hsp.subject.offset;
    let mut num_ident = 0;
    let mut align_length = 0;
    let mut num_pos = 0;
    let matrix = sbp
        .filter(|s| s.protein_alphabet)
        .map(|s| s.matrix.as_slice());
    if let Some(esp) = &hsp.gap_info {
        for index in 0..esp.size as usize {
            align_length += esp.num[index];
            match esp.op_type[index] {
                EGapAlignOpType::Sub => {
                    for _ in 0..esp.num[index] {
                        if query[q as usize] == subject[s as usize] {
                            num_ident += 1;
                        } else if let Some(m) = matrix {
                            if m[query[q as usize] as usize][subject[s as usize] as usize] > 0 {
                                num_pos += 1;
                            }
                        }
                        q += 1;
                        s += 1;
                    }
                }
                EGapAlignOpType::Del => s += esp.num[index],
                EGapAlignOpType::Ins => q += esp.num[index],
                _ => {
                    q += esp.num[index];
                    s += esp.num[index];
                }
            }
        }
    } else {
        if q_length != s_length {
            return -1;
        }
        align_length = q_length;
        for _ in 0..align_length {
            if query[q as usize] == subject[s as usize] {
                num_ident += 1;
            } else if let Some(m) = matrix {
                if m[query[q as usize] as usize][subject[s as usize] as usize] > 0 {
                    num_pos += 1;
                }
            }
            q += 1;
            s += 1;
        }
    }
    if let Some(out) = align_length_ptr {
        *out = align_length;
    }
    *num_ident_ptr = num_ident;
    if matrix.is_some() {
        *num_pos_ptr.expect("protein positives output required") = num_pos + num_ident;
    }
    0
}
pub fn blast_hsp_get_num_identities(
    query: Option<&[u8]>,
    subject: Option<&[u8]>,
    hsp: &mut BlastHSP,
    score_options: &BlastScoringOptions,
    align_length_ptr: Option<&mut i32>,
) -> i16 {
    assert!(
        !score_options.is_ooframe,
        "out-of-frame identity counting is outside standard scope"
    );
    let mut num_ident = hsp.num_ident;
    let status = s_blast_hsp_get_num_identities_and_positives(
        query,
        subject,
        hsp,
        &mut num_ident,
        align_length_ptr,
        None,
        None,
    );
    hsp.num_ident = num_ident;
    status
}
pub fn blast_hsp_get_num_identities_and_positives(
    query: Option<&[u8]>,
    subject: Option<&[u8]>,
    hsp: &mut BlastHSP,
    score_options: &BlastScoringOptions,
    align_length_ptr: Option<&mut i32>,
    sbp: Option<&BlastScoreBlk>,
) -> i16 {
    assert!(
        !score_options.is_ooframe,
        "out-of-frame identity counting is outside standard scope"
    );
    let mut num_ident = hsp.num_ident;
    let mut num_pos = hsp.num_positives;
    let status = s_blast_hsp_get_num_identities_and_positives(
        query,
        subject,
        hsp,
        &mut num_ident,
        align_length_ptr,
        sbp,
        Some(&mut num_pos),
    );
    hsp.num_ident = num_ident;
    hsp.num_positives = num_pos;
    status
}
pub fn s_hsp_test(hsp: &BlastHSP, hit_options: &BlastHitSavingOptions, align_length: i32) -> bool {
    hsp.num_ident as f64 * 100.0 < (align_length as f64 * hit_options.percent_identity)
        || align_length < hit_options.min_hit_length
}
pub fn blast_hsp_test_identity_and_length(
    _program: EBlastProgramType,
    hsp: &mut BlastHSP,
    query: &[u8],
    subject: &[u8],
    score_options: &BlastScoringOptions,
    hit_options: &BlastHitSavingOptions,
) -> bool {
    let mut align_length = 0;
    let status = blast_hsp_get_num_identities(
        Some(query),
        Some(subject),
        hsp,
        score_options,
        Some(&mut align_length),
    );
    assert_eq!(status, 0);
    s_hsp_test(hsp, hit_options, align_length)
}
pub fn blast_hsp_test(
    hsp: &BlastHSP,
    hit_options: &BlastHitSavingOptions,
    align_length: i32,
) -> bool {
    s_hsp_test(hsp, hit_options, align_length)
}
pub fn blast_hsp_get_query_coverage(hsp: &BlastHSP, query_length: i32) -> f64 {
    let mut pct = 0.0;
    if query_length > 0 {
        pct = 100.0 * (hsp.query.end - hsp.query.offset) as f64 / query_length as f64;
        if pct < 99.0 {
            pct += 0.5;
        }
    }
    pct
}
pub fn blast_hsp_query_coverage_test(
    hsp: &BlastHSP,
    min_query_coverage_pct: f64,
    query_length: i32,
) -> bool {
    blast_hsp_get_query_coverage(hsp, query_length) < min_query_coverage_pct
}

pub fn blast_hsp_list_get_evalues(
    program_number: EBlastProgramType,
    query_info: &BlastQueryInfo,
    subject_length: i32,
    hsp_list: Option<&mut BlastHSPList>,
    gapped_calculation: bool,
    rps_prelim: bool,
    sbp: &mut BlastScoreBlk,
    gap_decay_rate: f64,
    scaling_factor: f64,
) -> i16 {
    let Some(list) = hsp_list else { return 0 };
    if list.hspcnt == 0 {
        return 0;
    }
    let divisor = if gap_decay_rate != 0.0 {
        blast_gap_decay_divisor(gap_decay_rate, 1)
    } else {
        1.0
    };
    let is_rps = blast_program_is_rps_blast(program_number);
    let round_down = sbp.round_down;
    let gbp = sbp.gbp.as_ref();
    let kbp = if gapped_calculation {
        sbp.kbp_gap.as_mut().expect("gapped blocks required")
    } else if sbp.kbp_is_psi {
        &mut sbp.kbp_psi
    } else {
        &mut sbp.kbp_std
    };
    for slot in &mut list.hsp_array[..list.hspcnt as usize] {
        let hsp = slot.as_mut().expect("nonnull active HSP");
        assert!(scaling_factor != 0.0);
        let ctx = if rps_prelim {
            kbp.iter()
                .position(Option::is_some)
                .expect("valid Karlin block required")
        } else {
            hsp.context as usize
        };
        let karlin = kbp[ctx].as_mut().expect("valid Karlin block required");
        karlin.lambda /= scaling_factor;
        let score = if gapped_calculation && round_down {
            hsp.score & !1
        } else {
            hsp.score
        };
        let context = &query_info.contexts[hsp.context as usize];
        hsp.evalue = if let Some(gbp) = gbp {
            if !is_rps {
                blast_spouge_sto_e(score, karlin, gbp, context.query_length, subject_length)
            } else {
                blast_spouge_sto_e(score, karlin, gbp, subject_length, context.query_length)
            }
        } else {
            blast_karlin_sto_e_simple(score, karlin, context.eff_searchsp)
        };
        hsp.evalue /= divisor;
        karlin.lambda *= scaling_factor;
    }
    assert!(blast_hsp_list_is_sorted_by_score(Some(list)));
    list.best_evalue = s_blast_get_best_evalue(list);
    0
}

pub fn blast_hsp_list_reap_by_query_coverage(
    hsp_list: Option<&mut BlastHSPList>,
    hit_options: &BlastHitSavingOptions,
    query_info: &BlastQueryInfo,
    _program_number: EBlastProgramType,
) -> i16 {
    let Some(list) = hsp_list else { return 0 };
    if list.hspcnt == 0 || hit_options.query_cov_hsp_perc == 0.0 {
        return 0;
    }
    let mut count = 0;
    let mut removed = false;
    for index in 0..list.hspcnt as usize {
        let hsp = list.hsp_array[index].as_ref().expect("nonnull active HSP");
        if blast_hsp_query_coverage_test(
            hsp,
            hit_options.query_cov_hsp_perc,
            query_info.contexts[hsp.context as usize].query_length,
        ) {
            list.hsp_array[index] = None;
            removed = true;
        } else {
            if index > count {
                list.hsp_array[count] = list.hsp_array[index].take();
            }
            count += 1;
        }
    }
    list.hspcnt = count as i32;
    if removed {
        list.best_evalue = s_blast_get_best_evalue(list);
    }
    0
}

pub fn s_query_offset_compare_hs_ps(a: Option<&BlastHSP>, b: Option<&BlastHSP>) -> Ordering {
    match (a, b) {
        (None, None) => Ordering::Equal,
        (None, _) => Ordering::Greater,
        (_, None) => Ordering::Less,
        (Some(a), Some(b)) => a
            .context
            .cmp(&b.context)
            .then_with(|| a.query.offset.cmp(&b.query.offset))
            .then_with(|| a.subject.offset.cmp(&b.subject.offset))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.query.end.cmp(&a.query.end))
            .then_with(|| b.subject.end.cmp(&a.subject.end)),
    }
}
pub fn s_query_end_compare_hs_ps(a: Option<&BlastHSP>, b: Option<&BlastHSP>) -> Ordering {
    match (a, b) {
        (None, None) => Ordering::Equal,
        (None, _) => Ordering::Greater,
        (_, None) => Ordering::Less,
        (Some(a), Some(b)) => a
            .context
            .cmp(&b.context)
            .then_with(|| a.query.end.cmp(&b.query.end))
            .then_with(|| a.subject.end.cmp(&b.subject.end))
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| b.query.offset.cmp(&a.query.offset))
            .then_with(|| b.subject.offset.cmp(&a.subject.offset)),
    }
}

pub fn s_cut_off_gap_edit_script(
    hsp: &mut BlastHSP,
    mut q_cut: i32,
    mut s_cut: i32,
    cut_begin: bool,
) {
    let esp = hsp.gap_info.as_mut().expect("edit script required");
    let (mut qid, mut sid, mut opid, mut index, mut found) = (0, 0, 0, 0usize, false);
    q_cut -= hsp.query.offset;
    s_cut -= hsp.subject.offset;
    while index < esp.size as usize {
        opid = 0;
        while opid < esp.num[index] {
            match esp.op_type[index] {
                EGapAlignOpType::Sub => {
                    qid += 1;
                    sid += 1;
                    opid += 1;
                }
                EGapAlignOpType::Del => {
                    sid += esp.num[index];
                    opid += esp.num[index];
                }
                EGapAlignOpType::Ins => {
                    qid += esp.num[index];
                    opid += esp.num[index];
                }
                _ => panic!("ordinary edit script required"),
            }
            if qid >= q_cut && sid >= s_cut {
                found = true;
                break;
            }
        }
        if found {
            break;
        }
        index += 1;
    }
    if found {
        if cut_begin {
            let mut new_index = 0;
            if opid < esp.num[index] {
                assert_eq!(esp.op_type[index], EGapAlignOpType::Sub);
                esp.op_type[0] = esp.op_type[index];
                esp.num[0] = esp.num[index] - opid;
                new_index += 1;
            }
            index += 1;
            while index < esp.size as usize {
                esp.op_type[new_index] = esp.op_type[index];
                esp.num[new_index] = esp.num[index];
                index += 1;
                new_index += 1;
            }
            esp.size = new_index as i32;
            hsp.query.offset += qid;
            hsp.subject.offset += sid;
        } else {
            if opid < esp.num[index] {
                assert_eq!(esp.op_type[index], EGapAlignOpType::Sub);
                esp.num[index] = opid;
            }
            esp.size = index as i32 + 1;
            hsp.query.end = hsp.query.offset + qid;
            hsp.subject.end = hsp.subject.offset + sid;
        }
    }
}

pub fn blast_hsp_list_purge_hs_ps_with_common_endpoints(
    program: EBlastProgramType,
    hsp_list: Option<&mut BlastHSPList>,
    mut purge: bool,
) -> i32 {
    purge |= program != EBlastProgramType::Blastn;
    let Some(list) = hsp_list else { return 0 };
    if list.hspcnt == 0 {
        return 0;
    }
    if blast_program_is_phi_blast(program) {
        return list.hspcnt;
    }
    let mut count = list.hspcnt as usize;
    list.hsp_array[..count]
        .sort_by(|a, b| s_query_offset_compare_hs_ps(a.as_deref(), b.as_deref()));
    let mut i = 0;
    while i < count {
        while i + 1 < count
            && match (
                list.hsp_array[i].as_deref(),
                list.hsp_array[i + 1].as_deref(),
            ) {
                (Some(a), Some(b)) => {
                    a.context == b.context
                        && a.query.offset == b.query.offset
                        && a.subject.offset == b.subject.offset
                        && a.subject.frame == b.subject.frame
                }
                _ => false,
            }
        {
            count -= 1;
            let mut hsp = list.hsp_array[i + 1].take();
            let best = list.hsp_array[i].as_ref().unwrap();
            if !purge && hsp.as_ref().unwrap().query.end > best.query.end {
                s_cut_off_gap_edit_script(
                    hsp.as_mut().unwrap(),
                    best.query.end,
                    best.subject.end,
                    true,
                );
            } else {
                hsp = None;
            }
            for k in i + 1..count {
                list.hsp_array[k] = list.hsp_array[k + 1].take();
            }
            list.hsp_array[count] = hsp;
        }
        i += 1;
    }
    list.hsp_array[..count].sort_by(|a, b| s_query_end_compare_hs_ps(a.as_deref(), b.as_deref()));
    i = 0;
    while i < count {
        while i + 1 < count
            && match (
                list.hsp_array[i].as_deref(),
                list.hsp_array[i + 1].as_deref(),
            ) {
                (Some(a), Some(b)) => {
                    a.context == b.context
                        && a.query.end == b.query.end
                        && a.subject.end == b.subject.end
                        && a.subject.frame == b.subject.frame
                }
                _ => false,
            }
        {
            count -= 1;
            let mut hsp = list.hsp_array[i + 1].take();
            let best = list.hsp_array[i].as_ref().unwrap();
            if !purge && hsp.as_ref().unwrap().query.offset < best.query.offset {
                s_cut_off_gap_edit_script(
                    hsp.as_mut().unwrap(),
                    best.query.offset,
                    best.subject.offset,
                    false,
                );
            } else {
                hsp = None;
            }
            for k in i + 1..count {
                list.hsp_array[k] = list.hsp_array[k + 1].take();
            }
            list.hsp_array[count] = hsp;
        }
        i += 1;
    }
    if purge {
        blast_hsp_list_purge_null_hs_ps(Some(list));
    }
    count as i32
}
