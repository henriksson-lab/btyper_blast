// NCBI BLAST+ 2.17.0+, scoped default non-OOF tblastn blast_kappa.c callbacks.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_composition_adjustment::{ECompoAdjustModes, EMatrixAdjustRule};
use super::blast_composition_redo::*;
use super::blast_dynamic_align::{BlastDynamicAlign, DynamicSequence};
use super::blast_gapalign::{blast_gapped_alignment_with_traceback, BlastGapAlignStruct};
use super::blast_hits::{BlastHSP, BlastHSPList};
use super::blast_parameters::{BlastExtensionParameters, BlastScoringParameters};
use super::blast_program::EBlastProgramType;
use super::blast_seg::{seg_parameters_new_aa, seq_buffer_seg};
use super::blast_stat::BlastScoreBlk;
use super::blast_util::blast_get_partial_translation;
use super::gapinfo::{gap_prelim_edit_block_new, GapEditScript};
use super::ncbi_math::blast_nint;
use std::{collections::VecDeque, sync::Arc};

pub const BLASTP_MASK_RESIDUE: u8 = 21;
pub const BLASTP_MASK_INSTRUCTIONS: &str = "S 10 1.8 2.1";
pub const NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;
/// The source reuses void* for two distinct contexts. The enum preserves the
/// original HSP seed on input and owns the resulting edit script on output.
#[derive(Debug)]
pub enum BlastKappaAlignmentContext {
    Preliminary(BlastHSP),
    Traceback(GapEditScript),
}
/// Borrowed NCBI4na subject beginning at the first nucleotide, as expected by
/// Blast_GetPartialTranslation. Genetic-code entries are NCBIstdaa residues.
pub struct BlastKappaSubject<'a> {
    pub sequence_start: &'a [u8],
    pub length: i32,
    pub genetic_code: &'a [u8; 64],
}
pub struct BlastKappaCallbacks<'a> {
    pub subjects: &'a [BlastKappaSubject<'a>],
    pub scoring: &'a BlastScoringParameters,
    pub gap_align: Option<BlastGapAlignStruct<'static>>,
    /// Original Kappa thread-local matrix rows survive individual alignments.
    pub matrix_rows: &'a mut Vec<Vec<i32>>,
}
pub fn s_extend_right(
    query: &[u8],
    query_len: i32,
    subject: &[u8],
    subject_len: i32,
    max_shift: i32,
    query_ext_len: &mut i32,
    subject_ext_len: &mut i32,
    align_len: &mut i32,
) -> i32 {
    let (mut q, mut s, mut identical, mut gaps_q, mut gaps_s) = (0, 0, 0, 0, 0);
    while q < query_len && s < subject_len {
        while q < query_len && s < subject_len && query[q as usize] == subject[s as usize] {
            identical += 1;
            q += 1;
            s += 1;
        }
        let (mut n, mut matched) = (1, false);
        while n < max_shift && q + n + 1 < query_len && s + n + 1 < subject_len && !matched {
            if query[(q + n) as usize] == subject[(s + n) as usize]
                && query[(q + n + 1) as usize] == subject[(s + n + 1) as usize]
            {
                q += n + 2;
                s += n + 2;
                identical += 2;
                matched = true;
            }
            if !matched
                && query[(q + n) as usize] == subject[s as usize]
                && query[(q + n + 1) as usize] == subject[(s + 1) as usize]
            {
                q += n + 2;
                s += 2;
                identical += 2;
                gaps_s += n;
                matched = true;
            }
            if !matched
                && query[q as usize] == subject[(s + n) as usize]
                && query[(q + 1) as usize] == subject[(s + n + 1) as usize]
            {
                q += 2;
                s += n + 2;
                identical += 2;
                gaps_q += n;
                matched = true;
            }
            n += 1;
        }
        if matched {
            continue;
        }
        break;
    }
    *query_ext_len = q;
    *subject_ext_len = s;
    *align_len = if q > s { q + gaps_q } else { s + gaps_s };
    identical
}
pub fn s_extend_left(
    query: &[u8],
    query_len: i32,
    subject: &[u8],
    subject_len: i32,
    max_shift: i32,
    query_ext_len: &mut i32,
    subject_ext_len: &mut i32,
    align_len: &mut i32,
) -> i32 {
    let (mut q, mut s, mut identical, mut gaps_q, mut gaps_s) =
        (query_len - 1, subject_len - 1, 0, 0, 0);
    while q >= 0 && s >= 0 {
        while q > 0 && s > 0 && query[q as usize] == subject[s as usize] {
            identical += 1;
            q -= 1;
            s -= 1;
        }
        let (mut n, mut matched) = (1, false);
        while n < max_shift && q - n - 1 > 0 && s - n - 1 > 0 && !matched {
            if query[(q - n) as usize] == subject[(s - n) as usize]
                && query[(q - n - 1) as usize] == subject[(s - n - 1) as usize]
            {
                q -= n + 2;
                s -= n + 2;
                identical += 2;
                matched = true;
            }
            if !matched
                && query[(q - n) as usize] == subject[s as usize]
                && query[(q - n - 1) as usize] == subject[(s - 1) as usize]
            {
                q -= n + 2;
                s -= 2;
                identical += 2;
                gaps_s += n;
                matched = true;
            }
            if !matched
                && query[q as usize] == subject[(s - n) as usize]
                && query[(q - 1) as usize] == subject[(s - n - 1) as usize]
            {
                q -= 2;
                s -= n + 2;
                identical += 2;
                gaps_q += n;
                matched = true;
            }
            n += 1;
        }
        if matched {
            continue;
        }
        break;
    }
    *query_ext_len = query_len - q - 1;
    *subject_ext_len = subject_len - s - 1;
    *align_len += if *query_ext_len > *subject_ext_len {
        *query_ext_len + gaps_q
    } else {
        *subject_ext_len + gaps_s
    };
    identical
}
pub fn s_get_hash(data: &[u8], word_size: i32) -> u64 {
    let mut hash = 0u64;
    for k in 0..word_size {
        hash = hash.wrapping_shl(5).wrapping_add(data[k as usize] as u64);
    }
    hash
}
pub fn s_create_word_array(
    sequence: Option<&[u8]>,
    seq_len: i32,
    words: Option<&mut Vec<u64>>,
) -> i32 {
    let (Some(seq), Some(words)) = (sequence, words) else {
        return -1;
    };
    if seq_len < 8 {
        return -1;
    }
    *words = vec![0; (seq_len - 8 + 1) as usize];
    words[0] = s_get_hash(seq, 8);
    for i in 1..seq_len - 8 {
        words[i as usize] =
            ((words[i as usize - 1] << 5) & 0xffffffffff) + seq[(i + 7) as usize] as u64;
    }
    0
}
pub fn s_find_num_identical(
    query: &[u8],
    hashes: &[u64],
    query_len: i32,
    subject: &[u8],
    subject_len: i32,
    max_shift: i32,
) -> i32 {
    if query_len < 8 || subject_len < 8 || hashes.is_empty() {
        return 0;
    }
    let (mut hash, mut query_from, mut subject_from, mut s, mut identical, mut matched) =
        (0u64, 0, 0, 0, 0, false);
    while s < subject_len - 8 {
        hash = if s == 0 || matched {
            s_get_hash(&subject[s as usize..], 8)
        } else {
            ((hash << 5) & 0xffffffffff) + subject[(s + 7) as usize] as u64
        };
        let mut q = query_from;
        while q < query_len - 8 {
            if hashes[q as usize] == hash {
                break;
            }
            q += 1;
        }
        if q < query_len - 8 {
            let (mut ql, mut qr, mut sl, mut sr, mut al, mut ar) = (0, 0, 0, 0, 0, 0);
            matched = true;
            identical += 8;
            identical += s_extend_left(
                &query[query_from as usize..],
                q - query_from,
                &subject[subject_from as usize..],
                s - subject_from,
                max_shift,
                &mut ql,
                &mut sl,
                &mut al,
            );
            identical += s_extend_right(
                &query[(q + 8) as usize..],
                query_len - q - 8,
                &subject[(s + 8) as usize..],
                subject_len - s - 8,
                max_shift,
                &mut qr,
                &mut sr,
                &mut ar,
            );
            query_from = q + 8 + qr;
            subject_from = s + 8 + sr;
            s = subject_from - 1;
        } else {
            matched = false;
        }
        s += 1;
    }
    identical
}
pub fn s_test_near_identical<T>(
    seq: &BlastCompoSequenceData,
    seq_offset: i32,
    query: &BlastCompoSequenceData,
    query_offset: i32,
    words: &[u64],
    align: &BlastCompoAlignment<T>,
) -> bool {
    let qs = align.query_start - query_offset;
    let ss = align.match_start - seq_offset;
    let qlen = align.query_end - query_offset - qs;
    let slen = align.match_end - seq_offset - ss;
    let len = qlen.min(slen);
    let q = &query.data()[qs as usize..];
    let s = &seq.data()[ss as usize..];
    let (mut ql, mut qr, mut sl, mut sr, mut al, mut ar) = (0, 0, 0, 0, 0, 0);
    let mut identical = s_extend_right(q, qlen, s, slen, 8, &mut qr, &mut sr, &mut ar);
    if qr >= qlen || sr >= slen {
        let fraction = identical as f64 / len as f64;
        assert!(fraction - 1.0 < 1e-10);
        return fraction > 0.95;
    }
    identical += s_extend_left(
        &q[qr as usize..],
        qlen - qr,
        &s[sr as usize..],
        slen - sr,
        8,
        &mut ql,
        &mut sl,
        &mut al,
    );
    if ql + qr >= qlen || sl + sr >= slen {
        let fraction = identical as f64 / len as f64;
        assert!(fraction - 1.0 < 1e-10);
        return fraction > 0.95;
    }
    identical += s_find_num_identical(
        &q[qr as usize..],
        words.get((qs + qr) as usize..).unwrap_or(&[]),
        qlen - ql - qr,
        &s[sr as usize..],
        slen - sl - sr,
        8,
    );
    let fraction = identical as f64 / len as f64;
    assert!(fraction - 1.0 < 1e-10);
    fraction > 0.95
}
pub fn s_do_seg_sequence_data(data: &mut BlastCompoSequenceData, biased: Option<&mut bool>) -> i32 {
    let mut params = seg_parameters_new_aa();
    params.window = 10;
    params.locut = 1.8;
    params.hicut = 2.1;
    let mut masks = VecDeque::new();
    let status = seq_buffer_seg(
        &data.buffer[data.data_offset..],
        data.length,
        0,
        Some(&mut params),
        &mut masks,
    );
    if let Some(biased) = biased {
        *biased = !masks.is_empty();
    }
    if status == 0 {
        for mask in masks {
            for i in mask.left..=mask.right {
                data.data_mut()[i as usize] = BLASTP_MASK_RESIDUE;
            }
        }
    }
    status as i32
}
pub fn s_sequence_get_translated_range<T>(
    subject: &BlastKappaSubject<'_>,
    range: &BlastCompoSequenceRange,
    data: &mut BlastCompoSequenceData,
    q_range: &BlastCompoSequenceRange,
    query: &BlastCompoSequenceData,
    words: &[u64],
    align: &BlastCompoAlignment<T>,
    test: bool,
    mode: ECompoAdjustModes,
    biased: &mut bool,
) -> i32 {
    *data = BlastCompoSequenceData::default();
    let frame = range.context;
    let start = if frame > 0 {
        3 * range.begin
    } else {
        subject.length - 3 * range.end + frame + 1
    };
    let nucleotides = 3 * (range.end - range.begin) + frame.abs() - 1;
    let mut length = 0;
    let mut buffer = Vec::new();
    let mut status = blast_get_partial_translation(
        &subject.sequence_start[start as usize..],
        nucleotides,
        frame as i16,
        subject.genetic_code,
        Some(&mut buffer),
        Some(&mut length),
        None,
    );
    if status == 0 {
        *data = BlastCompoSequenceData {
            buffer,
            data_offset: 1,
            length,
        };
        if mode as i32 != 0
            && *biased
            && (!test
                || !s_test_near_identical(data, range.begin, query, q_range.begin, words, align))
        {
            status = s_do_seg_sequence_data(data, Some(biased));
            if status != 0 {
                *data = BlastCompoSequenceData::default();
            }
        }
    }
    status
}
pub fn s_sequence_get_range<T>(
    subject: &BlastKappaSubject<'_>,
    range: &BlastCompoSequenceRange,
    data: &mut BlastCompoSequenceData,
    query: &BlastCompoSequenceData,
    q_range: &BlastCompoSequenceRange,
    q_data: &mut BlastCompoSequenceData,
    words: &[u64],
    align: &BlastCompoAlignment<T>,
    test: bool,
    mode: ECompoAdjustModes,
    biased: &mut bool,
) -> i32 {
    *q_data = BlastCompoSequenceData::from_sequence(
        &query.data()[q_range.begin as usize..q_range.end as usize],
    );
    for residue in q_data.data_mut() {
        if *residue == 24 {
            *residue = 3;
        }
    }
    s_sequence_get_translated_range(
        subject, range, data, q_range, q_data, words, align, test, mode, biased,
    )
}
/// Rebind only the matrix borrow, preserving all dynamic traceback scratch.
fn bind_gap_align<'a>(
    gap: BlastGapAlignStruct<'_>,
    matrix: &'a [Vec<i32>],
) -> BlastGapAlignStruct<'a> {
    let d = gap.dynamic;
    BlastGapAlignStruct {
        dynamic: BlastDynamicAlign {
            matrix,
            pssm: None,
            position_based: false,
            gap_x_dropoff: d.gap_x_dropoff,
            dp_mem: d.dp_mem,
            dp_mem_alloc: d.dp_mem_alloc,
            state_struct: d.state_struct,
            query_start: d.query_start,
            query_stop: d.query_stop,
            subject_start: d.subject_start,
            subject_stop: d.subject_stop,
            score: d.score,
        },
        greedy_align_mem: gap.greedy_align_mem,
        edit_script: gap.edit_script,
        fwd_prelim_tback: gap.fwd_prelim_tback,
        rev_prelim_tback: gap.rev_prelim_tback,
        greedy_query_seed_start: gap.greedy_query_seed_start,
        greedy_subject_seed_start: gap.greedy_subject_seed_start,
    }
}
fn empty_gap_align() -> BlastGapAlignStruct<'static> {
    BlastGapAlignStruct {
        dynamic: BlastDynamicAlign::default(),
        greedy_align_mem: None,
        edit_script: None,
        fwd_prelim_tback: gap_prelim_edit_block_new(),
        rev_prelim_tback: gap_prelim_edit_block_new(),
        greedy_query_seed_start: 0,
        greedy_subject_seed_start: 0,
    }
}
pub fn s_new_alignment_from_gap_align(
    gap: &mut BlastGapAlignStruct<'_>,
    query: &BlastCompoSequenceRange,
    subject: &BlastCompoSequenceRange,
    rule: EMatrixAdjustRule,
) -> BlastCompoAlignment<BlastKappaAlignmentContext> {
    blast_compo_alignment_new(
        gap.dynamic.score,
        rule,
        gap.dynamic.query_start + query.begin,
        gap.dynamic.query_stop + query.begin,
        query.context,
        gap.dynamic.subject_start + subject.begin,
        gap.dynamic.subject_stop + subject.begin,
        subject.context,
        gap.edit_script
            .take()
            .map(|script| Arc::new(BlastKappaAlignmentContext::Traceback(script))),
    )
}
pub fn s_redo_one_alignment(
    callbacks: &mut BlastKappaCallbacks<'_>,
    input: &BlastCompoAlignment<BlastKappaAlignmentContext>,
    rule: EMatrixAdjustRule,
    query: &BlastCompoSequenceData,
    q_range: &BlastCompoSequenceRange,
    _ccat_length: i32,
    subject: &BlastCompoSequenceData,
    s_range: &BlastCompoSequenceRange,
    _subject_length: i32,
    gapping: &BlastCompoGappingParams,
    matrix: &[i32],
    alphsize: usize,
) -> Option<BlastCompoAlignment<BlastKappaAlignmentContext>> {
    let Some(BlastKappaAlignmentContext::Preliminary(hsp)) = input.context.as_deref() else {
        return None;
    };
    let q_start = hsp.query.gapped_start - q_range.begin;
    let s_start = hsp.subject.gapped_start - s_range.begin;
    assert_eq!(matrix.len(), alphsize * alphsize);
    let rows = &mut *callbacks.matrix_rows;
    rows.resize_with(alphsize, Vec::new);
    for (row, values) in rows.iter_mut().zip(matrix.chunks_exact(alphsize)) {
        row.resize(alphsize, 0);
        row.copy_from_slice(values);
    }
    let mut gap = bind_gap_align(
        callbacks.gap_align.take().unwrap_or_else(empty_gap_align),
        rows,
    );
    gap.dynamic.gap_x_dropoff = gapping.x_dropoff;
    let mut fence = false;
    let status = blast_gapped_alignment_with_traceback(
        EBlastProgramType::Tblastn,
        &DynamicSequence {
            data: &query.buffer,
            offset: query.data_offset as i32,
        },
        &DynamicSequence {
            data: &subject.buffer,
            offset: subject.data_offset as i32,
        },
        Some(&mut gap),
        callbacks.scoring,
        q_start,
        s_start,
        query.length,
        subject.length,
        Some(&mut fence),
    );
    let output = if status == 0 {
        Some(s_new_alignment_from_gap_align(
            &mut gap, q_range, s_range, rule,
        ))
    } else {
        None
    };
    callbacks.gap_align = Some(bind_gap_align(gap, &[]));
    output
}
impl BlastRedoAlignCallbacks<BlastKappaAlignmentContext> for BlastKappaCallbacks<'_> {
    fn get_range(
        &mut self,
        sequence: &BlastCompoMatchingSequence,
        range: &BlastCompoSequenceRange,
        data: &mut BlastCompoSequenceData,
        query: &BlastCompoSequenceData,
        q_range: &BlastCompoSequenceRange,
        q_data: &mut BlastCompoSequenceData,
        words: &[u64],
        align: &BlastCompoAlignment<BlastKappaAlignmentContext>,
        test: bool,
        mode: ECompoAdjustModes,
        is_sw: bool,
        biased: &mut bool,
    ) -> i32 {
        assert!(
            !is_sw,
            "default tblastn uses X-drop redo, not Smith-Waterman"
        );
        let Some(subject) = usize::try_from(sequence.index)
            .ok()
            .and_then(|i| self.subjects.get(i))
        else {
            return -1;
        };
        s_sequence_get_range(
            subject, range, data, query, q_range, q_data, words, align, test, mode, biased,
        )
    }
    fn redo_one_alignment(
        &mut self,
        input: &BlastCompoAlignment<BlastKappaAlignmentContext>,
        rule: EMatrixAdjustRule,
        query: &BlastCompoSequenceData,
        q_range: &BlastCompoSequenceRange,
        ccat: i32,
        subject: &BlastCompoSequenceData,
        s_range: &BlastCompoSequenceRange,
        full_length: i32,
        gapping: &BlastCompoGappingParams,
        matrix: &[i32],
        alphsize: usize,
    ) -> Option<BlastCompoAlignment<BlastKappaAlignmentContext>> {
        s_redo_one_alignment(
            self,
            input,
            rule,
            query,
            q_range,
            ccat,
            subject,
            s_range,
            full_length,
            gapping,
            matrix,
            alphsize,
        )
    }
}
pub fn s_rescale_search(
    sbp: &mut BlastScoreBlk,
    scoring: &mut BlastScoringParameters,
    num_queries: i32,
    scale_factor: f64,
) {
    if let Some(blocks) = &mut sbp.kbp_gap {
        for block in &mut blocks[..num_queries as usize] {
            if let Some(block) = block {
                block.lambda /= scale_factor;
                block.log_k = block.k.ln();
            }
        }
    }
    scoring.gap_open = blast_nint(scoring.gap_open as f64 * scale_factor) as i32;
    scoring.gap_extend = blast_nint(scoring.gap_extend as f64 * scale_factor) as i32;
    scoring.scale_factor = scale_factor;
}
pub fn s_gapping_params_new(
    sbp: &BlastScoreBlk,
    scoring: &BlastScoringParameters,
    ext: &BlastExtensionParameters,
    num_queries: i32,
) -> BlastCompoGappingParams {
    let mut min_lambda = f64::MAX;
    if let Some(blocks) = &sbp.kbp_gap {
        for block in &blocks[..num_queries as usize] {
            if let Some(block) = block {
                if block.lambda < min_lambda {
                    min_lambda = block.lambda;
                }
            }
        }
    }
    let bits = ext.options.gap_x_dropoff_final * std::f64::consts::LN_2 / min_lambda;
    BlastCompoGappingParams {
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        x_dropoff: if bits > ext.gap_x_dropoff_final as f64 {
            bits
        } else {
            ext.gap_x_dropoff_final as f64
        } as i32,
        decline_align: 0,
    }
}
pub fn s_hsp_list_normalize_scores(
    list: &mut BlastHSPList,
    lambda: f64,
    log_k: f64,
    score_divisor: f64,
) {
    for hsp in list.hsp_array[..list.hspcnt as usize].iter_mut().flatten() {
        hsp.score = blast_nint(hsp.score as f64 / score_divisor) as i32;
        hsp.bit_score =
            (hsp.score as f64 * lambda * score_divisor - log_k) / std::f64::consts::LN_2;
    }
}
pub fn s_get_query_info(
    query_data: &[u8],
    info: &super::blast_query_info::BlastQueryInfo,
    skip: bool,
) -> Vec<BlastCompoQueryInfo> {
    use super::blast_composition_adjustment::{
        blast_read_aa_composition, BlastAminoAcidComposition,
    };
    let mut output = Vec::with_capacity((info.last_context + 1) as usize);
    for context in &info.contexts[..(info.last_context + 1) as usize] {
        let seq = BlastCompoSequenceData::from_sequence(
            &query_data[context.query_offset as usize
                ..(context.query_offset + context.query_length) as usize],
        );
        let mut words = Vec::new();
        s_create_word_array(Some(seq.data()), seq.length, Some(&mut words));
        let mut composition = BlastAminoAcidComposition::default();
        if !skip {
            blast_read_aa_composition(&mut composition, 28, seq.data(), seq.length as usize);
        }
        output.push(BlastCompoQueryInfo {
            origin: context.query_offset,
            seq,
            composition,
            eff_search_space: context.eff_searchsp as f64,
            words,
        });
    }
    output
}
pub fn s_get_align_params(
    sbp: &BlastScoreBlk,
    scoring: &BlastScoringParameters,
    query_info: &super::blast_query_info::BlastQueryInfo,
    hit: &super::blast_parameters::BlastHitSavingParameters,
    ext: &BlastExtensionParameters,
    scale: f64,
) -> Option<BlastRedoAlignParams> {
    use super::blast_composition_adjustment::{blast_matrix_info_new, s_matrix_info_init};
    let mode = match ext.options.composition_based_stats {
        0 => ECompoAdjustModes::NoCompositionBasedStats,
        1 => ECompoAdjustModes::CompositionBasedStats,
        2 => ECompoAdjustModes::CompositionMatrixAdjust,
        3 => ECompoAdjustModes::CompoForceFullMatrixAdjust,
        _ => panic!("invalid composition mode"),
    };
    let mut near_cutoff = 0.0;
    for i in query_info.first_context..=query_info.last_context {
        if query_info.contexts[i as usize].is_valid {
            near_cutoff = NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2
                / sbp.kbp_gap.as_ref()?[i as usize].as_ref()?.lambda;
            break;
        }
    }
    let cutoff = if hit.do_sum_stats {
        assert!(hit.link_hsp_params.is_some());
        (hit.cutoff_score_min as f64 * scale) as i32
    } else {
        1
    };
    let mut matrix = blast_matrix_info_new(28, 28, false);
    if s_matrix_info_init(
        &mut matrix,
        sbp.kbp_ideal.as_ref()?.lambda,
        scale,
        b"BLOSUM62",
    ) != 0
    {
        return None;
    }
    Some(blast_redo_align_params_new(
        matrix,
        s_gapping_params_new(sbp, scoring, ext, query_info.last_context + 1),
        mode,
        false,
        false,
        true,
        query_info.max_length as i32,
        cutoff,
        hit.options.as_ref()?.expect_value,
        hit.do_sum_stats,
        near_cutoff,
    ))
}
/// Consumes the preliminary HSP owners: the C version retains borrowed HSP
/// pointers until redo finishes. Transferring ownership keeps those seeds live.
pub fn s_result_hsp_to_distinct_align(
    hsps: &mut [Option<Box<BlastHSP>>],
    hspcnt: i32,
    initial_context: i32,
    scale: f64,
) -> [Vec<BlastCompoAlignment<BlastKappaAlignmentContext>>; 6] {
    let mut output: [Vec<_>; 6] = std::array::from_fn(|_| Vec::new());
    for slot in &mut hsps[..hspcnt as usize] {
        let hsp = *slot.take().expect("live HSP prefix");
        let frame = hsp.context - initial_context;
        assert!((0..6).contains(&frame));
        output[frame as usize].push(blast_compo_alignment_new(
            (hsp.score as f64 * scale) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            hsp.query.offset,
            hsp.query.end,
            hsp.context,
            hsp.subject.offset,
            hsp.subject.end,
            hsp.subject.frame as i32,
            Some(Arc::new(BlastKappaAlignmentContext::Preliminary(hsp))),
        ));
    }
    output
}
/// As in C, returned alignment script contexts must have exclusive ownership.
pub fn s_hsp_list_from_distinct_alignments(
    list: Option<&mut BlastHSPList>,
    alignments: &mut Vec<BlastCompoAlignment<BlastKappaAlignmentContext>>,
    oid: i32,
    frame: i16,
) -> i32 {
    use super::blast_hits::{
        blast_hsp_init, blast_hsp_list_save_hsp, blast_hsp_list_sort_by_score,
    };
    let Some(list) = list else {
        return -1;
    };
    list.oid = oid;
    let mut status = 0;
    for align in alignments.iter_mut() {
        let mut script = align.context.take().map(|context| {
            match Arc::try_unwrap(context).expect("exclusive output edit-script owner") {
                BlastKappaAlignmentContext::Traceback(script) => script,
                _ => panic!("redo output requires traceback context"),
            }
        });
        let mut hsp = None;
        status = blast_hsp_init(
            align.query_start,
            align.query_end,
            align.match_start,
            align.match_end,
            0,
            0,
            align.query_index,
            frame,
            align.frame as i16,
            align.score,
            Some(&mut script),
            Some(&mut hsp),
        ) as i32;
        if status != 0 {
            break;
        }
        hsp.as_mut().unwrap().comp_adjustment_method = match align.matrix_adjust_rule {
            EMatrixAdjustRule::DontAdjustMatrix => 0,
            EMatrixAdjustRule::CompoScaleOldMatrix => 1,
            _ => 2,
        };
        hsp.as_mut().unwrap().num_ident = 0;
        status = blast_hsp_list_save_hsp(list, hsp) as i32;
        if status != 0 {
            break;
        }
    }
    if status == 0 {
        alignments.clear();
        blast_hsp_list_sort_by_score(Some(list));
    }
    // C returns zero even when it frees the destination on allocation failure.
    // Rust preserves the caller's list owner; allocation failure cannot arise in
    // HSPInit, and SaveHSP retains its defined bounded-insertion error behavior.
    0
}
