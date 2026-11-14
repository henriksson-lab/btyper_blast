// NCBI BLAST+ 2.17.0+, default redo_alignment.c dependency closure.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_composition_adjustment::{
    blast_adjust_scores, blast_get_composition_range, blast_read_aa_composition, s_calc_lambda,
    BlastAminoAcidComposition, BlastCompositionWorkspace, BlastMatrixInfo, ECompoAdjustModes,
    EMatrixAdjustRule,
};
use std::sync::Arc;
pub const LOCAL_LN2: f64 = 0.69314718055994530941723212145818;
pub const EVALUE_STRETCH: i32 = 5;
pub const K_WINDOW_BORDER: i32 = 200;
pub const K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS: i32 = 20;
pub const KAPPA_BIT_TOL: f64 = 2.0;
pub const MINIMUM_LENGTH_NEAR_IDENTICAL: i32 = 50;
#[derive(Debug)]
pub struct BlastCompoAlignment<T> {
    pub score: i32,
    pub matrix_adjust_rule: EMatrixAdjustRule,
    pub query_index: i32,
    pub query_start: i32,
    pub query_end: i32,
    pub match_start: i32,
    pub match_end: i32,
    pub frame: i32,
    pub context: Option<Arc<T>>,
}
impl<T> Clone for BlastCompoAlignment<T> {
    fn clone(&self) -> Self {
        Self {
            score: self.score,
            matrix_adjust_rule: self.matrix_adjust_rule,
            query_index: self.query_index,
            query_start: self.query_start,
            query_end: self.query_end,
            match_start: self.match_start,
            match_end: self.match_end,
            frame: self.frame,
            context: self.context.clone(),
        }
    }
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastCompoSequenceRange {
    pub begin: i32,
    pub end: i32,
    pub context: i32,
}
#[derive(Clone, Debug, Default)]
pub struct BlastCompoSequenceData {
    pub buffer: Vec<u8>,
    pub data_offset: usize,
    pub length: i32,
}
impl BlastCompoSequenceData {
    pub(crate) fn data(&self) -> &[u8] {
        &self.buffer[self.data_offset..self.data_offset + self.length as usize]
    }
    pub(crate) fn data_mut(&mut self) -> &mut [u8] {
        &mut self.buffer[self.data_offset..self.data_offset + self.length as usize]
    }
    pub(crate) fn from_sequence(sequence: &[u8]) -> Self {
        let mut buffer = Vec::with_capacity(sequence.len() + 2);
        buffer.push(0);
        buffer.extend_from_slice(sequence);
        buffer.push(0);
        Self {
            buffer,
            data_offset: 1,
            length: sequence.len() as i32,
        }
    }
}
#[derive(Debug)]
pub struct BlastCompoQueryInfo {
    pub origin: i32,
    pub seq: BlastCompoSequenceData,
    pub composition: BlastAminoAcidComposition,
    pub eff_search_space: f64,
    pub words: Vec<u64>,
}
#[derive(Clone, Copy, Debug)]
pub struct BlastCompoMatchingSequence {
    pub length: i32,
    pub index: i32,
}
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastCompoGappingParams {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub decline_align: i32,
    pub x_dropoff: i32,
}
#[derive(Debug)]
pub struct BlastRedoAlignParams {
    pub matrix_info: BlastMatrixInfo,
    pub gapping_params: BlastCompoGappingParams,
    pub compo_adjust_mode: ECompoAdjustModes,
    pub position_based: bool,
    pub re_pseudocounts: i32,
    pub query_is_translated: bool,
    pub subject_is_translated: bool,
    pub ccat_query_length: i32,
    pub cutoff_s: i32,
    pub cutoff_e: f64,
    pub do_link_hsps: bool,
    pub near_identical_cutoff: f64,
}
/// The source callback's opaque state is held by the implementing object.
/// The adjusted matrix is explicit instead of hidden in gapping_params.context.
pub trait BlastRedoAlignCallbacks<T> {
    fn calc_lambda(probs: &[f64], min_score: i32, max_score: i32, lambda0: f64) -> f64 {
        s_calc_lambda(probs, min_score, max_score, lambda0)
    }
    fn get_range(
        &mut self,
        sequence: &BlastCompoMatchingSequence,
        range: &BlastCompoSequenceRange,
        data: &mut BlastCompoSequenceData,
        orig_query: &BlastCompoSequenceData,
        q_range: &BlastCompoSequenceRange,
        q_data: &mut BlastCompoSequenceData,
        query_words: &[u64],
        align: &BlastCompoAlignment<T>,
        should_test_identical: bool,
        compo_adjust_mode: ECompoAdjustModes,
        is_smith_waterman: bool,
        subject_maybe_biased: &mut bool,
    ) -> i32;
    fn redo_one_alignment(
        &mut self,
        in_align: &BlastCompoAlignment<T>,
        matrix_adjust_rule: EMatrixAdjustRule,
        query_data: &BlastCompoSequenceData,
        query_range: &BlastCompoSequenceRange,
        ccat_query_length: i32,
        subject_data: &BlastCompoSequenceData,
        subject_range: &BlastCompoSequenceRange,
        full_subject_length: i32,
        gapping_params: &BlastCompoGappingParams,
        matrix: &[i32],
        alphsize: usize,
    ) -> Option<BlastCompoAlignment<T>>;
}
#[derive(Debug)]
pub struct SWindowInfo<T> {
    pub query_range: BlastCompoSequenceRange,
    pub subject_range: BlastCompoSequenceRange,
    pub align: Vec<BlastCompoAlignment<T>>,
    pub hspcnt: i32,
}
pub fn blast_compo_alignment_new<T>(
    score: i32,
    matrix_adjust_rule: EMatrixAdjustRule,
    query_start: i32,
    query_end: i32,
    query_index: i32,
    match_start: i32,
    match_end: i32,
    frame: i32,
    context: Option<Arc<T>>,
) -> BlastCompoAlignment<T> {
    BlastCompoAlignment {
        score,
        matrix_adjust_rule,
        query_start,
        query_end,
        query_index,
        match_start,
        match_end,
        frame,
        context,
    }
}
pub fn s_alignments_rev<T>(list: &mut [BlastCompoAlignment<T>]) {
    list.reverse();
}
pub fn s_alignment_cmp<T>(a: &BlastCompoAlignment<T>, b: &BlastCompoAlignment<T>) -> i32 {
    let result = b
        .score
        .cmp(&a.score)
        .then(a.match_start.cmp(&b.match_start))
        .then(b.match_end.cmp(&a.match_end))
        .then(a.query_start.cmp(&b.query_start))
        .then(b.query_end.cmp(&a.query_end));
    match result {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}
pub fn s_alignments_are_sorted<T>(alignments: &[BlastCompoAlignment<T>]) -> bool {
    alignments
        .windows(2)
        .all(|pair| pair[1].score <= pair[0].score)
}
pub fn s_distinct_alignments_length<T>(list: &[BlastCompoAlignment<T>]) -> i32 {
    list.len() as i32
}
/// Preserve source recursive mergesort, including right-before-left exact ties.
pub fn s_distinct_alignments_sort<T>(list: &mut Vec<BlastCompoAlignment<T>>, hspcnt: i32) {
    if hspcnt > 1 {
        let leftcnt = hspcnt / 2;
        let rightcnt = hspcnt - leftcnt;
        let mut right = list.split_off(leftcnt as usize);
        let mut left = std::mem::take(list);
        if leftcnt > 1 {
            s_distinct_alignments_sort(&mut left, leftcnt);
        }
        if rightcnt > 1 {
            s_distinct_alignments_sort(&mut right, rightcnt);
        }
        let mut left = left.into_iter().peekable();
        let mut right = right.into_iter().peekable();
        while left.peek().is_some() || right.peek().is_some() {
            let take_left = match (left.peek(), right.peek()) {
                (None, _) => false,
                (_, None) => true,
                (Some(a), Some(b)) => s_alignment_cmp(a, b) < 0,
            };
            list.push(if take_left {
                left.next().unwrap()
            } else {
                right.next().unwrap()
            });
        }
    }
}
pub fn s_alignment_copy<T>(align: &BlastCompoAlignment<T>) -> BlastCompoAlignment<T> {
    align.clone()
}
pub fn s_is_same_end_point<T>(
    new_align: &BlastCompoAlignment<T>,
    align: &BlastCompoAlignment<T>,
) -> bool {
    assert_eq!(new_align.frame, align.frame);
    (align.query_start == new_align.query_start && align.match_start == new_align.match_start)
        || (align.query_end == new_align.query_end && align.match_end == new_align.match_end)
}
fn contained(a: i32, b: i32, c: i32, d: i32, e: i32, f: i32) -> bool {
    (a <= c && b >= c) && (d <= f && e >= f)
}
pub fn s_is_similar_end_point<T>(
    new_align: &BlastCompoAlignment<T>,
    align: &BlastCompoAlignment<T>,
) -> bool {
    let start_contained = contained(
        align.query_start,
        align.query_end,
        new_align.query_start,
        align.match_start,
        align.match_end,
        new_align.match_start,
    );
    let end_contained = contained(
        align.query_start,
        align.query_end,
        new_align.query_end,
        align.match_start,
        align.match_end,
        new_align.match_end,
    );
    (start_contained
        && new_align.query_start - new_align.match_start == align.query_start - align.match_start)
        || (end_contained
            && new_align.query_end - new_align.match_end == align.query_end - align.match_end)
}
pub fn s_with_distinct_ends<T>(
    p_new_align: &mut Option<BlastCompoAlignment<T>>,
    old_alignments: &mut Vec<BlastCompoAlignment<T>>,
    is_same_adjustment: bool,
) {
    let new_align = p_new_align.take().unwrap();
    let mut include_new_align = true;
    for align in old_alignments.iter() {
        if align.frame == new_align.frame
            && ((is_same_adjustment && s_is_same_end_point(&new_align, align))
                || (!is_same_adjustment && s_is_similar_end_point(&new_align, align)))
            && new_align.score <= align.score
        {
            include_new_align = false;
            break;
        }
    }
    if include_new_align {
        old_alignments.retain(|align| {
            !(align.frame == new_align.frame
                && ((align.query_start == new_align.query_start
                    && align.match_start == new_align.match_start)
                    || (align.query_end == new_align.query_end
                        && align.match_end == new_align.match_end)))
        });
        old_alignments.insert(0, new_align);
    }
}
pub fn s_sequence_data_release(sequence: &mut BlastCompoSequenceData) {
    sequence.buffer = Vec::new();
    sequence.data_offset = 0;
}
pub fn s_window_info_new<T>(
    begin: i32,
    end: i32,
    context: i32,
    query_origin: i32,
    query_length: i32,
    query_index: i32,
    align: Vec<BlastCompoAlignment<T>>,
) -> SWindowInfo<T> {
    SWindowInfo {
        subject_range: BlastCompoSequenceRange {
            begin,
            end,
            context,
        },
        query_range: BlastCompoSequenceRange {
            begin: query_origin,
            end: query_origin + query_length,
            context: query_index,
        },
        hspcnt: align.len() as i32,
        align,
    }
}
pub fn s_window_swap_range<T>(window: &mut SWindowInfo<T>) {
    std::mem::swap(&mut window.query_range, &mut window.subject_range);
}
pub fn s_window_info_join<T>(win1: &mut SWindowInfo<T>, win2: SWindowInfo<T>) {
    assert_eq!(win1.subject_range.context, win2.subject_range.context);
    assert_eq!(win1.query_range.context, win2.query_range.context);
    win1.subject_range.begin = win1.subject_range.begin.min(win2.subject_range.begin);
    win1.subject_range.end = win1.subject_range.end.max(win2.subject_range.end);
    win1.hspcnt += win2.hspcnt;
    win1.align.extend(win2.align);
}
fn cmp_result(value: std::cmp::Ordering) -> i32 {
    match value {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}
pub fn s_location_compare_windows<T>(w1: &SWindowInfo<T>, w2: &SWindowInfo<T>) -> i32 {
    let (sr1, sr2, qr1, qr2) = (
        &w1.subject_range,
        &w2.subject_range,
        &w1.query_range,
        &w2.query_range,
    );
    cmp_result(
        qr1.context
            .cmp(&qr2.context)
            .then(sr1.context.cmp(&sr2.context))
            .then(sr1.begin.cmp(&sr2.begin))
            .then(sr1.end.cmp(&sr2.end))
            .then(qr1.begin.cmp(&qr2.begin))
            .then(qr1.end.cmp(&qr2.end)),
    )
}
pub fn s_subject_compare_windows<T>(w1: &SWindowInfo<T>, w2: &SWindowInfo<T>) -> i32 {
    let (sr1, sr2, qr1, qr2) = (
        &w1.subject_range,
        &w2.subject_range,
        &w1.query_range,
        &w2.query_range,
    );
    cmp_result(
        sr1.begin
            .cmp(&sr2.begin)
            .then(sr1.end.cmp(&sr2.end))
            .then(sr1.context.cmp(&sr2.context))
            .then(qr1.begin.cmp(&qr2.begin))
            .then(qr1.end.cmp(&qr2.end))
            .then(qr1.context.cmp(&qr2.context)),
    )
}
pub fn s_get_translated_length(length: i32, frame: i32, is_pos_based: bool) -> i32 {
    assert!(!is_pos_based, "PSSM outside default tblastn");
    (length - frame.abs() + 1) / 3
}
pub fn s_windows_from_translated_aligns<T>(
    alignments: &[BlastCompoAlignment<T>],
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    border: i32,
    sequence_length: i32,
    subject_is_translated: bool,
    is_pos_based: bool,
) -> Vec<SWindowInfo<T>> {
    assert_eq!(alignments.len(), hspcnt as usize);
    let mut windows = Vec::with_capacity(hspcnt as usize);
    for align in alignments {
        let frame = align.frame;
        let query_index = align.query_index;
        let query_length = query_info[query_index as usize].seq.length;
        let translated_length = s_get_translated_length(sequence_length, frame, is_pos_based);
        let align_copy = s_alignment_copy(align);
        windows.push(if subject_is_translated {
            s_window_info_new(
                0.max(align.match_start - border),
                translated_length.min(align.match_end + border),
                frame,
                0,
                query_length,
                query_index,
                vec![align_copy],
            )
        } else {
            s_window_info_new(
                0.max(align.query_start - border),
                query_length.min(align.query_end + border),
                query_index,
                0,
                sequence_length,
                0,
                vec![align_copy],
            )
        });
    }
    windows.sort_unstable_by(|a, b| s_location_compare_windows(a, b).cmp(&0));
    let mut joined = Vec::with_capacity(windows.len());
    let mut windows = windows.into_iter().peekable();
    while let Some(window) = windows.next() {
        if let Some(next) = windows.peek_mut() {
            if window.subject_range.context == next.subject_range.context
                && window.query_range.context == next.query_range.context
                && window.subject_range.end >= next.subject_range.begin
            {
                s_window_info_join(next, window);
                continue;
            }
        }
        joined.push(window);
    }
    if !subject_is_translated {
        for window in &mut joined {
            s_window_swap_range(window);
        }
    }
    for window in &mut joined {
        s_distinct_alignments_sort(&mut window.align, window.hspcnt);
    }
    joined.sort_unstable_by(|a, b| s_subject_compare_windows(a, b).cmp(&0));
    joined
}
pub fn s_windows_from_protein_aligns<T>(
    alignments: &[BlastCompoAlignment<T>],
    query_info: &[BlastCompoQueryInfo],
    num_queries: i32,
    sequence_length: i32,
) -> Vec<SWindowInfo<T>> {
    let mut windows: Vec<Option<SWindowInfo<T>>> = (0..num_queries).map(|_| None).collect();
    for align in alignments {
        let query_index = align.query_index as usize;
        let query_length = query_info[query_index].seq.length;
        let window = windows[query_index].get_or_insert_with(|| {
            s_window_info_new(
                0,
                sequence_length,
                0,
                0,
                query_length,
                query_index as i32,
                vec![],
            )
        });
        window.align.insert(0, s_alignment_copy(align));
        window.hspcnt += 1;
    }
    let mut windows: Vec<_> = windows.into_iter().flatten().collect();
    for window in &mut windows {
        s_alignments_rev(&mut window.align);
    }
    windows.sort_unstable_by(|a, b| s_subject_compare_windows(a, b).cmp(&0));
    windows
}
pub fn s_windows_from_aligns<T>(
    alignments: &[BlastCompoAlignment<T>],
    query_info: &[BlastCompoQueryInfo],
    hspcnt: i32,
    num_queries: i32,
    border: i32,
    sequence_length: i32,
    query_is_translated: bool,
    subject_is_translated: bool,
    is_pos_based: bool,
) -> Vec<SWindowInfo<T>> {
    if subject_is_translated || query_is_translated {
        s_windows_from_translated_aligns(
            alignments,
            query_info,
            hspcnt,
            border,
            sequence_length,
            subject_is_translated,
            is_pos_based,
        )
    } else {
        s_windows_from_protein_aligns(alignments, query_info, num_queries, sequence_length)
    }
}
pub fn s_get_composition<T>(
    composition: &mut BlastAminoAcidComposition,
    alphsize: usize,
    sequence: &BlastCompoSequenceData,
    range: &BlastCompoSequenceRange,
    align: &BlastCompoAlignment<T>,
    query_is_translated: bool,
    subject_is_translated: bool,
) {
    let length = range.end - range.begin;
    let (mut left, mut right) = (0, length);
    if query_is_translated || subject_is_translated {
        let start = if query_is_translated {
            align.query_start
        } else {
            align.match_start
        } - range.begin;
        let end = if query_is_translated {
            align.query_end
        } else {
            align.match_end
        } - range.begin;
        blast_get_composition_range(&mut left, &mut right, sequence.data(), length, start, end);
    }
    blast_read_aa_composition(
        composition,
        alphsize,
        &sequence.data()[left as usize..],
        (right - left) as usize,
    );
}
pub fn s_evalue_from_score(score: i32, lambda: f64, log_k: f64, searchsp: f64) -> f64 {
    searchsp * (-(lambda * score as f64) + log_k).exp()
}
pub fn s_is_contained<T>(
    in_align: &BlastCompoAlignment<T>,
    alignments: &[BlastCompoAlignment<T>],
    lambda: f64,
) -> bool {
    let score_thresh = in_align.score as f64 + KAPPA_BIT_TOL * LOCAL_LN2 / lambda;
    for align in alignments {
        if in_align.frame.signum() == align.frame.signum()
            && contained(
                align.query_start,
                align.query_end,
                in_align.query_start,
                align.match_start,
                align.match_end,
                in_align.match_start,
            )
            && contained(
                align.query_start,
                align.query_end,
                in_align.query_end,
                align.match_start,
                align.match_end,
                in_align.match_end,
            )
            && score_thresh <= align.score as f64
        {
            return true;
        }
    }
    false
}
pub fn blast_redo_align_params_new(
    matrix_info: BlastMatrixInfo,
    gapping_params: BlastCompoGappingParams,
    compo_adjust_mode: ECompoAdjustModes,
    position_based: bool,
    query_is_translated: bool,
    subject_is_translated: bool,
    ccat_query_length: i32,
    cutoff_s: i32,
    cutoff_e: f64,
    do_link_hsps: bool,
    near_identical_cutoff: f64,
) -> BlastRedoAlignParams {
    BlastRedoAlignParams {
        matrix_info,
        gapping_params,
        compo_adjust_mode,
        position_based,
        re_pseudocounts: K_RE_MATRIX_ADJUSTMENT_PSEUDOCOUNTS,
        query_is_translated,
        subject_is_translated,
        ccat_query_length,
        cutoff_s,
        cutoff_e,
        do_link_hsps,
        near_identical_cutoff,
    }
}
pub fn s_preliminary_test_near_identical<T>(
    query_info: &[BlastCompoQueryInfo],
    window: &SWindowInfo<T>,
    align: &BlastCompoAlignment<T>,
    cutoff: f64,
) -> bool {
    let query_length = query_info[align.query_index as usize].seq.length;
    if cutoff > 0.0 {
        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
        let align_len =
            (align.query_end - align.query_start).min(align.match_end - align.match_start);
        if align.score as f64 / (align_len as f64) < cutoff {
            return false;
        }
    } else {
        if window.hspcnt > 1 || window.hspcnt < 1 {
            return false;
        }
        if align.query_end - align.query_start != align.match_end - align.match_start {
            return false;
        }
        if align.match_end - align.match_start + 1 < query_length.min(MINIMUM_LENGTH_NEAR_IDENTICAL)
        {
            return false;
        }
    }
    true
}
pub fn blast_redo_one_match<T, C: BlastRedoAlignCallbacks<T>>(
    alignments: &mut [Vec<BlastCompoAlignment<T>>],
    params: &BlastRedoAlignParams,
    incoming_aligns: &[BlastCompoAlignment<T>],
    hspcnt: i32,
    lambda: f64,
    matching_seq: &BlastCompoMatchingSequence,
    ccat_query_length: i32,
    query_info: &mut [BlastCompoQueryInfo],
    num_queries: i32,
    matrix: &mut [i32],
    alphsize: usize,
    workspace: &mut BlastCompositionWorkspace,
    pvalue_for_this_pair: &mut f64,
    composition_test_index: i32,
    lambda_ratio: &mut f64,
    callbacks: &mut C,
) -> i32 {
    let mut status = 0;
    let mut matrix_adjust_rule = EMatrixAdjustRule::DontAdjustMatrix;
    assert!((params.compo_adjust_mode as i32) < 2 || !params.position_based);
    for query_index in 0..num_queries as usize {
        alignments[query_index].clear();
    }
    let windows = s_windows_from_aligns(
        incoming_aligns,
        query_info,
        hspcnt,
        num_queries,
        K_WINDOW_BORDER,
        matching_seq.length,
        params.query_is_translated,
        params.subject_is_translated,
        params.position_based,
    );
    for window in &windows {
        let mut subject = BlastCompoSequenceData::default();
        let mut query = BlastCompoSequenceData::default();
        let mut near_identical_status = true;
        let mut old_near_identical_status = false;
        let mut subject_maybe_biased = true;
        let mut num_adjustments = 0;
        let query_index = window.align[0].query_index as usize;
        for (hsp_index, in_align) in window.align.iter().enumerate() {
            if hsp_index == 0 || subject_maybe_biased {
                near_identical_status = s_preliminary_test_near_identical(
                    query_info,
                    window,
                    in_align,
                    params.near_identical_cutoff,
                );
            }
            if hsp_index == 0
                || (subject_maybe_biased && near_identical_status != old_near_identical_status)
            {
                s_sequence_data_release(&mut subject);
                s_sequence_data_release(&mut query);
                status = callbacks.get_range(
                    matching_seq,
                    &window.subject_range,
                    &mut subject,
                    &query_info[query_index].seq,
                    &window.query_range,
                    &mut query,
                    &query_info[query_index].words,
                    &window.align[0],
                    near_identical_status,
                    params.compo_adjust_mode,
                    false,
                    &mut subject_maybe_biased,
                );
                if status != 0 {
                    break;
                }
            }
            if params.query_is_translated {
                s_get_composition(
                    &mut query_info[query_index].composition,
                    alphsize,
                    &query,
                    &window.query_range,
                    in_align,
                    true,
                    false,
                );
            }
            if !s_is_contained(in_align, &alignments[query_index], lambda) {
                let mut adjust_search_failed = 0;
                if params.compo_adjust_mode != ECompoAdjustModes::NoCompositionBasedStats
                    && (params.subject_is_translated
                        || hsp_index == 0
                        || near_identical_status != old_near_identical_status)
                {
                    let mut subject_composition = BlastAminoAcidComposition::default();
                    s_get_composition(
                        &mut subject_composition,
                        alphsize,
                        &subject,
                        &window.subject_range,
                        in_align,
                        false,
                        params.subject_is_translated,
                    );
                    adjust_search_failed = blast_adjust_scores(
                        matrix,
                        &query_info[query_index].composition,
                        query.length,
                        &subject_composition,
                        subject.length,
                        &params.matrix_info,
                        params.compo_adjust_mode,
                        params.re_pseudocounts,
                        workspace,
                        &mut matrix_adjust_rule,
                        C::calc_lambda,
                        pvalue_for_this_pair,
                        composition_test_index,
                        lambda_ratio,
                    );
                    if adjust_search_failed < 0 {
                        status = adjust_search_failed;
                        break;
                    }
                    num_adjustments += 1;
                }
                if adjust_search_failed == 0 {
                    let mut new_align = callbacks.redo_one_alignment(
                        in_align,
                        matrix_adjust_rule,
                        &query,
                        &window.query_range,
                        ccat_query_length,
                        &subject,
                        &window.subject_range,
                        matching_seq.length,
                        &params.gapping_params,
                        matrix,
                        alphsize,
                    );
                    if new_align
                        .as_ref()
                        .is_some_and(|a| a.score >= params.cutoff_s)
                    {
                        s_with_distinct_ends(
                            &mut new_align,
                            &mut alignments[query_index],
                            num_adjustments == 1,
                        );
                    }
                }
            }
            old_near_identical_status = near_identical_status;
        }
        if status != 0 {
            break;
        }
    }
    if status != 0 {
        for query_index in 0..num_queries as usize {
            alignments[query_index].clear();
        }
    }
    status
}
pub fn blast_compo_early_termination<T>(
    evalue: f64,
    significant_matches: &[super::compo_heap::BlastCompoHeap<T>],
    num_queries: i32,
) -> bool {
    for i in 0..num_queries as usize {
        if super::compo_heap::blast_compo_heap_filled_to_cutoff(&significant_matches[i]) {
            if evalue <= EVALUE_STRETCH as f64 * significant_matches[i].ecutoff {
                return false;
            }
        } else {
            return false;
        }
    }
    true
}
