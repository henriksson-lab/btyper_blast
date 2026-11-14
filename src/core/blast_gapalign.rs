// NCBI BLAST+ 2.17.0+, selected blast_gapalign.c/h orchestration.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_dynamic_align::{align_ex, BlastDynamicAlign, BlastGapDP, DynamicSequence};
use super::blast_options::EBlastPrelimGapExt;
use super::blast_parameters::{BlastExtensionParameters, BlastScoringParameters};
use super::blast_program::EBlastProgramType;
use super::blast_stat::BlastScoreBlk;
use super::gapinfo::{
    gap_edit_script_new, gap_prelim_edit_block_new, gap_prelim_edit_block_reset, EGapAlignOpType,
    GapEditScript, GapPrelimEditBlock,
};
use super::greedy_align::{
    blast_affine_greedy_align, s_blast_greedy_align_mem_alloc, GreedyExtensionParameters,
    GreedyScoringParameters, SGreedyAlignMem, SGreedySeed, GREEDY_MAX_COST,
    GREEDY_MAX_COST_FRACTION,
};

pub const MAX_DBSEQ_LEN: u32 = 5_000_000;
/// Scoped owner for ordinary nucleotide/protein alignment. Dynamic scratch
/// stores the shared matrix, X-drop, score and output coordinates.
#[derive(Debug)]
pub struct BlastGapAlignStruct<'a> {
    pub dynamic: BlastDynamicAlign<'a>,
    pub greedy_align_mem: Option<SGreedyAlignMem>,
    pub edit_script: Option<GapEditScript>,
    pub fwd_prelim_tback: GapPrelimEditBlock,
    pub rev_prelim_tback: GapPrelimEditBlock,
    pub greedy_query_seed_start: i32,
    pub greedy_subject_seed_start: i32,
}
pub fn blast_gap_align_struct_new<'a>(
    score_params: Option<&BlastScoringParameters>,
    ext_params: Option<&BlastExtensionParameters>,
    mut max_subject_length: u32,
    sbp: Option<&'a BlastScoreBlk>,
    gap_align_ptr: Option<&mut Option<BlastGapAlignStruct<'a>>>,
) -> i16 {
    let (Some(score), Some(ext), Some(sbp), Some(output)) =
        (score_params, ext_params, sbp, gap_align_ptr)
    else {
        return -1;
    };
    assert_ne!(
        ext.options.e_prelim_gap_ext,
        EBlastPrelimGapExt::JumperWithTraceback,
        "jumper alignment is outside standard search scope"
    );
    assert!(
        !ext.options.chaining,
        "chaining is outside standard search scope"
    );
    let mut dynamic = BlastDynamicAlign {
        matrix: &sbp.matrix,
        gap_x_dropoff: ext.gap_x_dropoff,
        ..Default::default()
    };
    let greedy_align_mem = if ext.options.e_prelim_gap_ext == EBlastPrelimGapExt::DynProgScoreOnly {
        dynamic.dp_mem_alloc = 1000;
        dynamic.dp_mem = vec![BlastGapDP::default(); 1000];
        None
    } else {
        max_subject_length = max_subject_length.min(MAX_DBSEQ_LEN);
        max_subject_length =
            (max_subject_length / GREEDY_MAX_COST_FRACTION as u32 + 1).min(GREEDY_MAX_COST as u32);
        // The source applies the GREEDY_MAX_COST clamp after division.
        let cfg = GreedyScoringParameters {
            reward: score.reward,
            penalty: score.penalty,
            gap_open: score.gap_open,
            gap_extend: score.gap_extend,
        };
        let extcfg = GreedyExtensionParameters {
            gap_x_dropoff: ext.gap_x_dropoff,
            gap_x_dropoff_final: ext.gap_x_dropoff_final,
        };
        let Some(mem) =
            s_blast_greedy_align_mem_alloc(Some(&cfg), Some(&extcfg), max_subject_length as i32, 0)
        else {
            *output = None;
            return -1;
        };
        Some(mem)
    };
    *output = Some(BlastGapAlignStruct {
        dynamic,
        greedy_align_mem,
        edit_script: None,
        fwd_prelim_tback: gap_prelim_edit_block_new(),
        rev_prelim_tback: gap_prelim_edit_block_new(),
        greedy_query_seed_start: 0,
        greedy_subject_seed_start: 0,
    });
    0
}
pub fn blast_gap_align_struct_free(
    gap_align: Option<BlastGapAlignStruct<'_>>,
) -> Option<BlastGapAlignStruct<'_>> {
    drop(gap_align);
    None
}

pub fn blast_prelim_edit_block_to_gap_edit_script(
    rev: Option<&GapPrelimEditBlock>,
    fwd: Option<&GapPrelimEditBlock>,
) -> Option<GapEditScript> {
    let (Some(rev), Some(fwd)) = (rev, fwd) else {
        return None;
    };
    let merge = fwd.num_ops > 0
        && rev.num_ops > 0
        && fwd.edit_ops[fwd.num_ops as usize - 1].op_type
            == rev.edit_ops[rev.num_ops as usize - 1].op_type;
    let size = fwd.num_ops + rev.num_ops - if merge { 1 } else { 0 };
    let mut esp = gap_edit_script_new(size)?;
    let mut index = 0;
    for i in 0..rev.num_ops as usize {
        esp.op_type[index] = rev.edit_ops[i].op_type;
        esp.num[index] = rev.edit_ops[i].num;
        index += 1;
    }
    if fwd.num_ops == 0 {
        return Some(esp);
    }
    if merge {
        esp.num[index - 1] += fwd.edit_ops[fwd.num_ops as usize - 1].num;
    }
    let mut i = fwd.num_ops - if merge { 2 } else { 1 };
    while i >= 0 {
        esp.op_type[index] = fwd.edit_ops[i as usize].op_type;
        esp.num[index] = fwd.edit_ops[i as usize].num;
        index += 1;
        i -= 1;
    }
    Some(esp)
}
pub fn s_blast_greedy_gap_align_struct_fill(
    gap_align: &mut BlastGapAlignStruct<'_>,
    q_start: i32,
    s_start: i32,
    q_end: i32,
    s_end: i32,
    q_seed_start: i32,
    s_seed_start: i32,
    score: i32,
    esp: Option<GapEditScript>,
) -> i16 {
    gap_align.dynamic.query_start = q_start;
    gap_align.dynamic.query_stop = q_end;
    gap_align.dynamic.subject_start = s_start;
    gap_align.dynamic.subject_stop = s_end;
    gap_align.greedy_query_seed_start = q_seed_start;
    gap_align.greedy_subject_seed_start = s_seed_start;
    gap_align.dynamic.score = score;
    if esp.is_some() {
        gap_align.edit_script = esp;
    }
    0
}

pub fn s_update_edit_script(esp: &mut GapEditScript, pos: i32, bf: i32, af: i32) {
    if bf > 0 {
        let mut op = pos;
        let mut qd = bf;
        let mut sd = bf;
        loop {
            op -= 1;
            if op < 0 {
                return;
            }
            match esp.op_type[op as usize] {
                EGapAlignOpType::Sub => {
                    qd -= esp.num[op as usize];
                    sd -= esp.num[op as usize];
                }
                EGapAlignOpType::Ins => qd -= esp.num[op as usize],
                EGapAlignOpType::Del => sd -= esp.num[op as usize],
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.num[op as usize] = -qd.max(sd);
        esp.op_type[op as usize] = EGapAlignOpType::Sub;
        op += 1;
        while op < pos - 1 {
            esp.num[op as usize] = 0;
            op += 1;
        }
        esp.num[pos as usize] += bf;
        qd -= sd;
        esp.op_type[(pos - 1) as usize] = if qd > 0 {
            EGapAlignOpType::Del
        } else {
            EGapAlignOpType::Ins
        };
        esp.num[(pos - 1) as usize] = qd.abs();
    }
    if af > 0 {
        let mut op = pos;
        let mut qd = af;
        let mut sd = af;
        loop {
            op += 1;
            if op >= esp.size {
                return;
            }
            match esp.op_type[op as usize] {
                EGapAlignOpType::Sub => {
                    qd -= esp.num[op as usize];
                    sd -= esp.num[op as usize];
                }
                EGapAlignOpType::Ins => qd -= esp.num[op as usize],
                EGapAlignOpType::Del => sd -= esp.num[op as usize],
                _ => {}
            }
            if qd <= 0 && sd <= 0 {
                break;
            }
        }
        esp.num[op as usize] = -qd.max(sd);
        esp.op_type[op as usize] = EGapAlignOpType::Sub;
        op -= 1;
        while op > pos + 1 {
            esp.num[op as usize] = 0;
            op -= 1;
        }
        esp.num[pos as usize] += af;
        qd -= sd;
        esp.op_type[(pos + 1) as usize] = if qd > 0 {
            EGapAlignOpType::Del
        } else {
            EGapAlignOpType::Ins
        };
        esp.num[(pos + 1) as usize] = qd.abs();
    }
}
pub fn s_rebuild_edit_script(esp: &mut GapEditScript) {
    let mut j = -1i32;
    for i in 0..esp.size as usize {
        if esp.num[i] == 0 {
            continue;
        }
        if j >= 0 && esp.op_type[i] == esp.op_type[j as usize] {
            esp.num[j as usize] += esp.num[i];
        } else if j == -1
            || esp.op_type[i] == EGapAlignOpType::Sub
            || esp.op_type[j as usize] == EGapAlignOpType::Sub
        {
            j += 1;
            esp.op_type[j as usize] = esp.op_type[i];
            esp.num[j as usize] = esp.num[i];
        } else {
            let d = esp.num[j as usize] - esp.num[i];
            if d > 0 {
                esp.num[(j - 1) as usize] += esp.num[i];
                esp.num[j as usize] = d;
            } else if d < 0 {
                if j == 0 && i as i32 - j > 0 {
                    esp.op_type[j as usize] = EGapAlignOpType::Sub;
                    j += 1;
                } else {
                    esp.num[(j - 1) as usize] += esp.num[j as usize];
                }
                esp.num[j as usize] = -d;
                esp.op_type[j as usize] = esp.op_type[i];
            } else {
                esp.num[(j - 1) as usize] += esp.num[j as usize];
                j -= 1;
            }
        }
    }
    esp.size = j + 1;
}
pub fn s_reduce_gaps(
    esp: &mut GapEditScript,
    query: &DynamicSequence<'_>,
    subject: &DynamicSequence<'_>,
    q_end: i32,
    s_end: i32,
) {
    let mut q1 = 0;
    let mut s1 = 0;
    for i in 0..esp.size as usize {
        if esp.num[i] == 0 {
            continue;
        }
        if esp.op_type[i] == EGapAlignOpType::Sub {
            if esp.num[i] >= 12 {
                let mut nm1 = 1;
                if i > 0 {
                    while q1 - nm1 >= 0 && query.at(q1 - nm1) == subject.at(s1 - nm1) {
                        nm1 += 1;
                    }
                }
                q1 += esp.num[i];
                s1 += esp.num[i];
                let mut nm2 = 0;
                if i < esp.size as usize - 1 {
                    while q1 + 1 < q_end && s1 + 1 < s_end {
                        let equal = query.at(q1) == subject.at(s1);
                        q1 += 1;
                        s1 += 1;
                        if !equal {
                            break;
                        }
                        nm2 += 1;
                    }
                }
                if nm1 > 1 || nm2 > 0 {
                    s_update_edit_script(esp, i as i32, nm1 - 1, nm2);
                }
                q1 -= 1;
                s1 -= 1;
            } else {
                q1 += esp.num[i];
                s1 += esp.num[i];
            }
        } else if esp.op_type[i] == EGapAlignOpType::Ins {
            q1 += esp.num[i];
        } else {
            s1 += esp.num[i];
        }
    }
    s_rebuild_edit_script(esp);
    let mut q = 0;
    let mut s = 0;
    for i in 0..esp.size as usize {
        if esp.op_type[i] == EGapAlignOpType::Sub {
            q += esp.num[i];
            s += esp.num[i];
            continue;
        }
        if i > 1 && esp.op_type[i] != esp.op_type[i - 2] && esp.num[i - 2] > 0 {
            let mut d = esp.num[i] + esp.num[i - 1] + esp.num[i - 2];
            if d == 3 {
                esp.num[i - 2] = 0;
                esp.num[i - 1] = 2;
                esp.num[i] = 0;
                if esp.op_type[i] == EGapAlignOpType::Ins {
                    q += 1;
                } else {
                    s += 1;
                }
            } else if d < 12 {
                let mut nm1 = 0;
                let mut nm2 = 0;
                d = esp.num[i].min(esp.num[i - 2]);
                q -= esp.num[i - 1];
                s -= esp.num[i - 1];
                q1 = q;
                s1 = s;
                if esp.op_type[i] == EGapAlignOpType::Ins {
                    s -= d;
                } else {
                    q -= d;
                }
                for _ in 0..esp.num[i - 1] {
                    if query.at(q1) == subject.at(s1) {
                        nm1 += 1;
                    }
                    if query.at(q) == subject.at(s) {
                        nm2 += 1;
                    }
                    q1 += 1;
                    s1 += 1;
                    q += 1;
                    s += 1;
                }
                for _ in 0..d {
                    if query.at(q) == subject.at(s) {
                        nm2 += 1;
                    }
                    q += 1;
                    s += 1;
                }
                if nm2 >= nm1 - d {
                    esp.num[i - 2] -= d;
                    esp.num[i - 1] += d;
                    esp.num[i] -= d;
                } else {
                    q = q1;
                    s = s1;
                }
            }
        }
        if esp.op_type[i] == EGapAlignOpType::Ins {
            q += esp.num[i];
        } else {
            s += esp.num[i];
        }
    }
    s_rebuild_edit_script(esp);
}

pub fn blast_greedy_gapped_alignment(
    query: &[u8],
    subject: &[u8],
    query_length: i32,
    subject_length: i32,
    gap_align: &mut BlastGapAlignStruct<'_>,
    score_params: &BlastScoringParameters,
    q_off: i32,
    s_off: i32,
    compressed_subject: bool,
    do_traceback: bool,
    mut fence_hit: Option<&mut bool>,
) -> i16 {
    assert!(
        query_length >= 0
            && subject_length >= 0
            && q_off >= 0
            && q_off <= query_length
            && s_off >= 0
            && s_off <= subject_length
    );
    let q_avail = query_length - q_off;
    let s_avail = subject_length - s_off;
    let q = &query[q_off as usize..];
    let (s, mut rem) = if compressed_subject {
        (&subject[(s_off / 4) as usize..], (s_off % 4) as u8)
    } else {
        (&subject[s_off as usize..], 4)
    };
    let x = gap_align.dynamic.gap_x_dropoff;
    if do_traceback {
        gap_prelim_edit_block_reset(Some(&mut gap_align.fwd_prelim_tback));
        gap_prelim_edit_block_reset(Some(&mut gap_align.rev_prelim_tback));
    }
    let cfg = GreedyScoringParameters {
        reward: score_params.reward,
        penalty: score_params.penalty,
        gap_open: score_params.gap_open,
        gap_extend: score_params.gap_extend,
    };
    let mut q_ext_r = 0;
    let mut s_ext_r = 0;
    let mut q_ext_l = 0;
    let mut s_ext_l = 0;
    let mut fwd_start_point = SGreedySeed::default();
    let mut rev_start_point = SGreedySeed::default();
    let mut score;
    loop {
        score = blast_affine_greedy_align(
            q,
            q_avail,
            s,
            s_avail,
            false,
            x,
            score_params.reward,
            -score_params.penalty,
            score_params.gap_open,
            score_params.gap_extend,
            &mut q_ext_r,
            &mut s_ext_r,
            gap_align
                .greedy_align_mem
                .as_mut()
                .expect("greedy scratch required"),
            if do_traceback {
                Some(&mut gap_align.fwd_prelim_tback)
            } else {
                None
            },
            rem,
            fence_hit.as_deref_mut(),
            &mut fwd_start_point,
        );
        if fence_hit.as_deref().copied() == Some(true) {
            return -1;
        }
        if score >= 0 {
            break;
        }
        let mem = gap_align.greedy_align_mem.as_ref().unwrap();
        let new_dist = mem.max_dist * 2;
        let xdrop = mem.xdrop;
        gap_align.greedy_align_mem = None;
        gap_align.greedy_align_mem =
            s_blast_greedy_align_mem_alloc(Some(&cfg), None, new_dist, xdrop);
        if gap_align.greedy_align_mem.is_none() {
            return -1;
        }
    }
    if compressed_subject {
        rem = 0;
    }
    loop {
        let score1 = blast_affine_greedy_align(
            query,
            q_off,
            subject,
            s_off,
            true,
            x,
            score_params.reward,
            -score_params.penalty,
            score_params.gap_open,
            score_params.gap_extend,
            &mut q_ext_l,
            &mut s_ext_l,
            gap_align.greedy_align_mem.as_mut().unwrap(),
            if do_traceback {
                Some(&mut gap_align.rev_prelim_tback)
            } else {
                None
            },
            rem,
            fence_hit.as_deref_mut(),
            &mut rev_start_point,
        );
        if fence_hit.as_deref().copied() == Some(true) {
            return -1;
        }
        if score1 >= 0 {
            score += score1;
            break;
        }
        let mem = gap_align.greedy_align_mem.as_ref().unwrap();
        let new_dist = mem.max_dist * 2;
        let xdrop = mem.xdrop;
        gap_align.greedy_align_mem = None;
        gap_align.greedy_align_mem =
            s_blast_greedy_align_mem_alloc(Some(&cfg), None, new_dist, xdrop);
        if gap_align.greedy_align_mem.is_none() {
            return -1;
        }
    }
    if score_params.gap_open == 0 && score_params.gap_extend == 0 {
        score = (q_ext_r + s_ext_r + q_ext_l + s_ext_l) * score_params.reward / 2
            - score * (score_params.reward - score_params.penalty);
    } else if score_params.reward % 2 == 1 {
        score /= 2;
    }
    let mut q_seed_start = q_off;
    let mut s_seed_start = s_off;
    let esp;
    if do_traceback {
        let mut script = blast_prelim_edit_block_to_gap_edit_script(
            Some(&gap_align.rev_prelim_tback),
            Some(&gap_align.fwd_prelim_tback),
        );
        assert!(
            !compressed_subject,
            "greedy final traceback requires the unpacked subject"
        );
        if let Some(e) = &mut script {
            s_reduce_gaps(
                e,
                &DynamicSequence {
                    data: query,
                    offset: q_off - q_ext_l,
                },
                &DynamicSequence {
                    data: subject,
                    offset: s_off - s_ext_l,
                },
                q_ext_l + q_ext_r,
                s_ext_l + s_ext_r,
            );
        }
        esp = script;
    } else {
        let q_box_l = q_off - q_ext_l;
        let s_box_l = s_off - s_ext_l;
        let q_box_r = q_off + q_ext_r;
        let s_box_r = s_off + s_ext_r;
        let mut q_seed_start_l = q_off - rev_start_point.start_q;
        let mut s_seed_start_l = s_off - rev_start_point.start_s;
        let mut q_seed_start_r = q_off + fwd_start_point.start_q;
        let mut s_seed_start_r = s_off + fwd_start_point.start_s;
        let mut valid_seed_len_l = 0;
        let mut valid_seed_len_r = 0;
        if q_seed_start_r < q_box_r && s_seed_start_r < s_box_r {
            valid_seed_len_r = (q_box_r - q_seed_start_r)
                .min(s_box_r - s_seed_start_r)
                .min(fwd_start_point.match_length)
                / 2;
        } else {
            q_seed_start_r = q_off;
            s_seed_start_r = s_off;
        }
        if q_seed_start_l > q_box_l && s_seed_start_l > s_box_l {
            valid_seed_len_l = (q_seed_start_l - q_box_l)
                .min(s_seed_start_l - s_box_l)
                .min(rev_start_point.match_length)
                / 2;
        } else {
            q_seed_start_l = q_off;
            s_seed_start_l = s_off;
        }
        if valid_seed_len_r > valid_seed_len_l {
            q_seed_start = q_seed_start_r + valid_seed_len_r;
            s_seed_start = s_seed_start_r + valid_seed_len_r;
        } else {
            q_seed_start = q_seed_start_l - valid_seed_len_l;
            s_seed_start = s_seed_start_l - valid_seed_len_l;
        }
        esp = None;
    }
    s_blast_greedy_gap_align_struct_fill(
        gap_align,
        q_off - q_ext_l,
        s_off - s_ext_l,
        q_off + q_ext_r,
        s_off + s_ext_r,
        q_seed_start,
        s_seed_start,
        score,
        esp,
    );
    0
}

/// Non-out-of-frame standard traceback path. Sequence views retain original
/// query/subject sequence-block pointer origins and surrounding sentinels.
pub fn blast_gapped_alignment_with_traceback(
    _program: EBlastProgramType,
    query: &DynamicSequence<'_>,
    subject: &DynamicSequence<'_>,
    gap_align: Option<&mut BlastGapAlignStruct<'_>>,
    score_params: &BlastScoringParameters,
    q_start: i32,
    s_start: i32,
    query_length: i32,
    subject_length: i32,
    mut fence_hit: Option<&mut bool>,
) -> i16 {
    assert!(
        !score_params.options.is_ooframe,
        "out-of-frame alignment is outside standard nucleotide scope"
    );
    let Some(gap_align) = gap_align else {
        return -1;
    };
    gap_prelim_edit_block_reset(Some(&mut gap_align.fwd_prelim_tback));
    gap_prelim_edit_block_reset(Some(&mut gap_align.rev_prelim_tback));
    let mut private_q_length = 0;
    let mut private_s_length = 0;
    let mut score_left = align_ex(
        query,
        subject,
        q_start + 1,
        s_start + 1,
        &mut private_q_length,
        &mut private_s_length,
        &mut gap_align.rev_prelim_tback,
        &mut gap_align.dynamic,
        score_params,
        q_start,
        false,
        true,
        fence_hit.as_deref_mut(),
    );
    gap_align.dynamic.query_start = q_start - private_q_length + 1;
    gap_align.dynamic.subject_start = s_start - private_s_length + 1;
    let mut score_right = 0;
    let mut found_end = false;
    if fence_hit.as_deref().copied() != Some(true)
        && q_start < query_length
        && s_start < subject_length
    {
        found_end = true;
        score_right = align_ex(
            &query.shifted(q_start),
            &subject.shifted(s_start),
            query_length - q_start - 1,
            subject_length - s_start - 1,
            &mut private_q_length,
            &mut private_s_length,
            &mut gap_align.fwd_prelim_tback,
            &mut gap_align.dynamic,
            score_params,
            q_start,
            false,
            false,
            fence_hit.as_deref_mut(),
        );
        gap_align.dynamic.query_stop = q_start + private_q_length + 1;
        gap_align.dynamic.subject_stop = s_start + private_s_length + 1;
    }
    if !found_end {
        gap_align.dynamic.query_stop = q_start - 1;
        gap_align.dynamic.subject_stop = s_start - 1;
    }
    gap_align.edit_script = blast_prelim_edit_block_to_gap_edit_script(
        Some(&gap_align.rev_prelim_tback),
        Some(&gap_align.fwd_prelim_tback),
    );
    if let Some(esp) = &mut gap_align.edit_script {
        while esp.size > 0 && esp.op_type[0] != EGapAlignOpType::Sub {
            score_left += score_params.gap_open + esp.num[0] * score_params.gap_extend;
            if esp.op_type[0] == EGapAlignOpType::Del {
                gap_align.dynamic.subject_start += esp.num[0];
            } else {
                gap_align.dynamic.query_start += esp.num[0];
            }
            for i in 1..esp.size as usize {
                esp.op_type[i - 1] = esp.op_type[i];
                esp.num[i - 1] = esp.num[i];
            }
            esp.size -= 1;
        }
        let mut i = esp.size;
        while esp.size > 0 && esp.op_type[(i - 1) as usize] != EGapAlignOpType::Sub {
            score_right +=
                score_params.gap_open + esp.num[(i - 1) as usize] * score_params.gap_extend;
            if esp.op_type[(i - 1) as usize] == EGapAlignOpType::Del {
                gap_align.dynamic.subject_stop -= esp.num[(i - 1) as usize];
            } else {
                gap_align.dynamic.query_stop -= esp.num[(i - 1) as usize];
            }
            esp.size -= 1;
            i -= 1;
            assert_eq!(esp.size, i);
        }
    }
    gap_align.dynamic.score = score_right + score_left;
    0
}
