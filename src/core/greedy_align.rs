// NCBI BLAST+ 2.17.0+, greedy_align.c/h and selected scratch allocator.
// Author: Ilya Dondoshansky. See LICENSE-NCBI.
use super::blast_util::FENCE_SENTRY;
use super::gapinfo::{gap_prelim_edit_block_add, EGapAlignOpType, GapPrelimEditBlock};
use super::ncbi_math::blast_gdb3;

pub const GREEDY_MAX_COST_FRACTION: i32 = 2;
pub const GREEDY_MAX_COST: i32 = 1000;
pub const K_INVALID_OFFSET: i32 = -2;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct SGreedyOffset {
    pub insert_off: i32,
    pub match_off: i32,
    pub delete_off: i32,
}
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct SGreedySeed {
    pub start_q: i32,
    pub start_s: i32,
    pub match_length: i32,
}

/// Stable pool buffers preserve the original million-cell block size.
#[derive(Debug)]
pub struct SMBSpace {
    space_array: Vec<SGreedyOffset>,
    space_allocated: i32,
    space_used: i32,
    next: Option<Box<SMBSpace>>,
}
pub fn mb_space_new(num_space_arrays: i32) -> Box<SMBSpace> {
    let count = num_space_arrays.max(1_000_000);
    Box::new(SMBSpace {
        space_array: vec![SGreedyOffset::default(); count as usize],
        space_allocated: count,
        space_used: 0,
        next: None,
    })
}
pub fn s_refresh_mb_space(mut space: Option<&mut SMBSpace>) {
    while let Some(current) = space {
        current.space_used = 0;
        space = current.next.as_deref_mut();
    }
}
pub fn mb_space_free(mut space: Option<Box<SMBSpace>>) {
    while let Some(mut current) = space {
        space = current.next.take();
    }
}
pub fn s_get_mb_space(mut pool: &mut SMBSpace, num_alloc: i32) -> Option<*mut SGreedyOffset> {
    if num_alloc < 0 {
        return None;
    }
    while pool.space_used + num_alloc > pool.space_allocated {
        if pool.next.is_none() {
            pool.next = Some(mb_space_new(num_alloc));
        }
        pool = pool.next.as_deref_mut().unwrap();
    }
    let out = pool
        .space_array
        .as_mut_ptr()
        .wrapping_offset(pool.space_used as isize);
    pool.space_used += num_alloc;
    Some(out)
}

/// A biased original row pointer with its allocated logical bounds. Buffers
/// live in the owning scratch allocation or pool; biased addresses are formed
/// with wrapping_offset and never dereferenced before the bounds check.
#[derive(Clone, Copy, Debug)]
struct GreedyRow<T> {
    ptr: *mut T,
    low: i32,
    high: i32,
}
impl<T: Copy> GreedyRow<T> {
    fn empty() -> Self {
        Self {
            ptr: std::ptr::null_mut(),
            low: 1,
            high: 0,
        }
    }
    fn get(&self, index: i32) -> T {
        assert!(index >= self.low && index <= self.high);
        unsafe { *self.ptr.wrapping_offset(index as isize) }
    }
    fn set(&mut self, index: i32, value: T) {
        assert!(index >= self.low && index <= self.high);
        unsafe {
            *self.ptr.wrapping_offset(index as isize) = value;
        }
    }
    fn field(&mut self, index: i32) -> &mut T {
        assert!(index >= self.low && index <= self.high);
        unsafe { &mut *self.ptr.wrapping_offset(index as isize) }
    }
}
#[derive(Debug)]
pub struct SGreedyAlignMem {
    pub max_dist: i32,
    pub xdrop: i32,
    last_seq2_off: Vec<GreedyRow<i32>>,
    max_score: Vec<i32>,
    last_seq2_off_affine: Vec<GreedyRow<SGreedyOffset>>,
    diag_bounds: Vec<i32>,
    space: Option<Box<SMBSpace>>,
    _initial_non_affine: Vec<i32>,
    _initial_affine: Vec<SGreedyOffset>,
    configured_scores: (i32, i32, i32, i32),
}
// Row pointers refer only to exclusively owned stable buffers. A move between
// threads does not move their allocations; mutation requires &mut scratch.
unsafe impl Send for SGreedyAlignMem {}
impl Drop for SGreedyAlignMem {
    fn drop(&mut self) {
        mb_space_free(self.space.take());
    }
}

#[derive(Clone, Copy, Debug)]
pub struct GreedyScoringParameters {
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}
#[derive(Clone, Copy, Debug)]
pub struct GreedyExtensionParameters {
    pub gap_x_dropoff: i32,
    pub gap_x_dropoff_final: i32,
}

pub fn s_blast_greedy_align_mem_alloc(
    score_params: Option<&GreedyScoringParameters>,
    ext_params: Option<&GreedyExtensionParameters>,
    mut max_d: i32,
    mut xdrop: i32,
) -> Option<SGreedyAlignMem> {
    let score_params = score_params?;
    if ext_params.is_none() && xdrop == 0 {
        return None;
    }
    assert!(
        max_d > 0
            && score_params.reward > 0
            && score_params.penalty <= 0
            && score_params.gap_open >= 0
            && score_params.gap_extend >= 0
    );
    let (reward, penalty, gap_open, mut gap_extend) = if score_params.reward % 2 == 1 {
        if xdrop == 0 {
            let e = ext_params.unwrap();
            xdrop = 2 * e.gap_x_dropoff.max(e.gap_x_dropoff_final);
        }
        (
            2 * score_params.reward,
            -2 * score_params.penalty,
            2 * score_params.gap_open,
            2 * score_params.gap_extend,
        )
    } else {
        if xdrop == 0 {
            let e = ext_params.unwrap();
            xdrop = e.gap_x_dropoff.max(e.gap_x_dropoff_final);
        }
        (
            score_params.reward,
            -score_params.penalty,
            score_params.gap_open,
            score_params.gap_extend,
        )
    };
    let configured_scores = (reward, penalty, gap_open, gap_extend);
    if gap_open == 0 && gap_extend == 0 {
        gap_extend = reward / 2 + penalty;
    }
    let width = 2 * max_d + 6;
    let mut out = SGreedyAlignMem {
        max_dist: max_d,
        xdrop,
        last_seq2_off: vec![],
        max_score: vec![],
        last_seq2_off_affine: vec![],
        diag_bounds: vec![],
        space: None,
        _initial_non_affine: vec![],
        _initial_affine: vec![],
        configured_scores,
    };
    let d_diff;
    if score_params.gap_open == 0 && score_params.gap_extend == 0 {
        d_diff = (xdrop + reward / 2) / (penalty + reward) + 1;
        out.last_seq2_off = vec![GreedyRow::empty(); (max_d + 2) as usize];
        out._initial_non_affine = vec![0; (width * 2) as usize];
        let ptr = out._initial_non_affine.as_mut_ptr();
        out.last_seq2_off[0] = GreedyRow {
            ptr,
            low: 0,
            high: width - 1,
        };
        out.last_seq2_off[1] = GreedyRow {
            ptr: ptr.wrapping_offset(width as isize),
            low: 0,
            high: width - 1,
        };
    } else {
        let mut mis_cost = reward + penalty;
        let mut ge_cost = gap_extend + reward / 2;
        let mut open = gap_open;
        max_d *= ge_cost;
        let max_cost = mis_cost.max(open + ge_cost);
        let gd = blast_gdb3(&mut mis_cost, &mut open, &mut ge_cost);
        d_diff = (xdrop + reward / 2) / gd + 1;
        out.diag_bounds = vec![0; (2 * (max_d + 1 + max_cost)) as usize];
        out.last_seq2_off_affine = vec![GreedyRow::empty(); (max_d.max(max_cost) + 2) as usize];
        out._initial_affine = vec![SGreedyOffset::default(); (width * (max_cost + 1)) as usize];
        for i in 0..=max_cost {
            out.last_seq2_off_affine[i as usize] = GreedyRow {
                ptr: out
                    ._initial_affine
                    .as_mut_ptr()
                    .wrapping_offset((i * width) as isize),
                low: 0,
                high: width - 1,
            };
        }
    }
    out.max_score = vec![0; (max_d + 1 + d_diff) as usize];
    out.space = Some(mb_space_new(0));
    Some(out)
}

pub fn s_find_first_mismatch(
    seq1: &[u8],
    seq2: &[u8],
    len1: i32,
    len2: i32,
    mut seq1_index: i32,
    mut seq2_index: i32,
    mut fence_hit: Option<&mut bool>,
    reverse: bool,
    rem: u8,
) -> i32 {
    assert!(
        len1 >= 0
            && len2 >= 0
            && len1 as usize <= seq1.len()
            && seq1_index >= 0
            && seq2_index >= 0
            && rem <= 4
    );
    let tmp = seq1_index;
    if reverse {
        if rem == 4 {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[(len1 - 1 - seq1_index) as usize] < 4
                && seq1[(len1 - 1 - seq1_index) as usize] == seq2[(len2 - 1 - seq2_index) as usize]
            {
                seq1_index += 1;
                seq2_index += 1;
            }
            if seq2_index < len2 && seq2[(len2 - 1 - seq2_index) as usize] == FENCE_SENTRY {
                *fence_hit
                    .as_deref_mut()
                    .expect("fence output required for a sentry") = true;
            }
        } else {
            while seq1_index < len1
                && seq2_index < len2
                && seq1[(len1 - 1 - seq1_index) as usize]
                    == ((seq2[((len2 - 1 - seq2_index) / 4) as usize]
                        >> (2 * (3 - (len2 - 1 - seq2_index) % 4)))
                        & 3)
            {
                seq1_index += 1;
                seq2_index += 1;
            }
        }
    } else if rem == 4 {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[seq1_index as usize] < 4
            && seq1[seq1_index as usize] == seq2[seq2_index as usize]
        {
            seq1_index += 1;
            seq2_index += 1;
        }
        if seq2_index < len2 && seq2[seq2_index as usize] == FENCE_SENTRY {
            *fence_hit
                .as_deref_mut()
                .expect("fence output required for a sentry") = true;
        }
    } else {
        while seq1_index < len1
            && seq2_index < len2
            && seq1[seq1_index as usize]
                == ((seq2[((seq2_index + rem as i32) / 4) as usize]
                    >> (2 * (3 - (seq2_index + rem as i32) % 4)))
                    & 3)
        {
            seq1_index += 1;
            seq2_index += 1;
        }
    }
    seq1_index - tmp
}

pub fn s_get_next_non_affine_tback(
    aux: &SGreedyAlignMem,
    d: i32,
    diag: i32,
    seq2_index: &mut i32,
) -> i32 {
    let row = &aux.last_seq2_off[(d - 1) as usize];
    if row.get(diag - 1) > row.get(diag).max(row.get(diag + 1)) {
        *seq2_index = row.get(diag - 1);
        return diag - 1;
    }
    if row.get(diag) > row.get(diag + 1) {
        *seq2_index = row.get(diag);
        return diag;
    }
    *seq2_index = row.get(diag + 1);
    diag + 1
}

pub fn s_get_next_affine_tback_from_match(
    aux: &SGreedyAlignMem,
    diag_lower: &[i32],
    diag_upper: &[i32],
    bounds_offset: i32,
    d: &mut i32,
    diag: i32,
    op_cost: i32,
    seq2_index: &mut i32,
) -> EGapAlignOpType {
    let old = *d - op_cost;
    if diag >= diag_lower[(old + bounds_offset) as usize]
        && diag <= diag_upper[(old + bounds_offset) as usize]
    {
        let new = aux.last_seq2_off_affine[old as usize].get(diag).match_off;
        let cell = aux.last_seq2_off_affine[*d as usize].get(diag);
        if new >= cell.insert_off.max(cell.delete_off) {
            *d -= op_cost;
            *seq2_index = new;
            return EGapAlignOpType::Sub;
        }
    }
    let cell = aux.last_seq2_off_affine[*d as usize].get(diag);
    if cell.insert_off > cell.delete_off {
        *seq2_index = cell.insert_off;
        EGapAlignOpType::Ins
    } else {
        *seq2_index = cell.delete_off;
        EGapAlignOpType::Del
    }
}
pub fn s_get_next_affine_tback_from_indel(
    aux: &SGreedyAlignMem,
    diag_lower: &[i32],
    diag_upper: &[i32],
    bounds_offset: i32,
    d: &mut i32,
    diag: i32,
    gap_open: i32,
    gap_extend: i32,
    i_or_d: EGapAlignOpType,
) -> EGapAlignOpType {
    let new_diag = if i_or_d == EGapAlignOpType::Ins {
        diag - 1
    } else {
        diag + 1
    };
    let mut last_d = *d - gap_extend;
    let new_index = if new_diag >= diag_lower[(last_d + bounds_offset) as usize]
        && new_diag <= diag_upper[(last_d + bounds_offset) as usize]
    {
        let c = aux.last_seq2_off_affine[last_d as usize].get(new_diag);
        if i_or_d == EGapAlignOpType::Ins {
            c.insert_off
        } else {
            c.delete_off
        }
    } else {
        K_INVALID_OFFSET
    };
    last_d = *d - gap_open - gap_extend;
    if new_diag >= diag_lower[(last_d + bounds_offset) as usize]
        && new_diag <= diag_upper[(last_d + bounds_offset) as usize]
        && new_index
            < aux.last_seq2_off_affine[last_d as usize]
                .get(new_diag)
                .match_off
    {
        *d -= gap_open + gap_extend;
        return EGapAlignOpType::Sub;
    }
    assert_ne!(new_index, K_INVALID_OFFSET);
    *d -= gap_extend;
    i_or_d
}

pub fn blast_greedy_align(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    reverse: bool,
    xdrop_threshold: i32,
    match_cost: i32,
    mismatch_cost: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    aux: &mut SGreedyAlignMem,
    mut edit_block: Option<&mut GapPrelimEditBlock>,
    rem: u8,
    mut fence_hit: Option<&mut bool>,
    seed: &mut SGreedySeed,
) -> i32 {
    assert!(match_cost >= 2 && mismatch_cost >= 0 && xdrop_threshold >= 0);
    assert_eq!(aux.configured_scores, (match_cost, mismatch_cost, 0, 0));
    let max_dist = aux.max_dist;
    let diag_origin = max_dist + 2;
    let xdrop_offset = (xdrop_threshold + match_cost / 2) / (match_cost + mismatch_cost) + 1;
    assert!(aux.max_score.len() >= (max_dist + 1 + xdrop_offset) as usize);
    let mut best_dist = 0;
    let mut best_diag = 0;
    let initial = s_find_first_mismatch(
        seq1,
        seq2,
        len1,
        len2,
        0,
        0,
        fence_hit.as_deref_mut(),
        reverse,
        rem,
    );
    *seq1_align_len = initial;
    *seq2_align_len = initial;
    seed.start_q = 0;
    seed.start_s = 0;
    seed.match_length = initial;
    let mut longest_match_run = initial;
    if initial == len1 || initial == len2 {
        if let Some(edit) = edit_block.as_deref_mut() {
            gap_prelim_edit_block_add(edit, EGapAlignOpType::Sub, initial);
        }
        return 0;
    }
    if edit_block.is_some() {
        if aux.space.is_none() {
            aux.space = Some(mb_space_new(0));
        } else {
            s_refresh_mb_space(aux.space.as_deref_mut());
        }
    }
    for i in 0..xdrop_offset {
        aux.max_score[i as usize] = 0;
    }
    aux.last_seq2_off[0].set(diag_origin, initial);
    aux.max_score[xdrop_offset as usize] = initial * match_cost;
    let mut diag_lower = diag_origin - 1;
    let mut diag_upper = diag_origin + 1;
    let mut end1_reached = false;
    let mut end2_reached = false;
    let mut converged = false;
    for d in 1..=max_dist {
        let mut curr_extent = 0;
        let mut curr_seq2_index = 0;
        let mut curr_diag = 0;
        let tmp_diag_lower = diag_lower;
        let tmp_diag_upper = diag_upper;
        for k in [diag_lower - 1, diag_lower, diag_upper, diag_upper + 1] {
            aux.last_seq2_off[(d - 1) as usize].set(k, K_INVALID_OFFSET);
        }
        let xdrop_score =
            aux.max_score[d as usize] + (match_cost + mismatch_cost) * d - xdrop_threshold;
        let xdrop_score = (xdrop_score as f64 / (match_cost / 2) as f64).ceil() as i32;
        for k in tmp_diag_lower..=tmp_diag_upper {
            let previous = aux.last_seq2_off[(d - 1) as usize];
            let mut seq2_index =
                (previous.get(k + 1).max(previous.get(k)) + 1).max(previous.get(k - 1));
            let mut seq1_index = seq2_index + k - diag_origin;
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == diag_lower {
                    diag_lower += 1;
                } else {
                    aux.last_seq2_off[d as usize].set(k, K_INVALID_OFFSET);
                }
                continue;
            }
            diag_upper = k;
            let index = s_find_first_mismatch(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index,
                seq2_index,
                fence_hit.as_deref_mut(),
                reverse,
                rem,
            );
            if fence_hit.as_deref().copied() == Some(true) {
                return 0;
            }
            if index > longest_match_run {
                seed.start_q = seq1_index;
                seed.start_s = seq2_index;
                seed.match_length = index;
                longest_match_run = index;
            }
            seq1_index += index;
            seq2_index += index;
            aux.last_seq2_off[d as usize].set(k, seq2_index);
            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }
            if seq2_index == len2 {
                diag_lower = k + 1;
                end2_reached = true;
            }
            if seq1_index == len1 {
                diag_upper = k - 1;
                end1_reached = true;
            }
        }
        let curr_score = curr_extent * (match_cost / 2) - d * (match_cost + mismatch_cost);
        if curr_score >= aux.max_score[(d - 1 + xdrop_offset) as usize] {
            aux.max_score[(d + xdrop_offset) as usize] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            aux.max_score[(d + xdrop_offset) as usize] =
                aux.max_score[(d - 1 + xdrop_offset) as usize];
        }
        if diag_lower > diag_upper {
            converged = true;
            break;
        }
        if !end2_reached {
            diag_lower -= 1;
        }
        if !end1_reached {
            diag_upper += 1;
        }
        if edit_block.is_none() {
            aux.last_seq2_off[(d + 1) as usize] = aux.last_seq2_off[(d - 1) as usize];
        } else {
            let count = (diag_upper - diag_lower + 7) / 3;
            let ptr = s_get_mb_space(aux.space.as_deref_mut().unwrap(), count).unwrap() as *mut i32;
            let low = diag_lower - 2;
            aux.last_seq2_off[(d + 1) as usize] = GreedyRow {
                ptr: ptr.wrapping_offset((-low) as isize),
                low,
                high: low + 3 * count - 1,
            };
        }
    }
    if !converged {
        return -1;
    }
    let Some(edit) = edit_block else {
        return best_dist;
    };
    let mut d = best_dist;
    let mut seq2_index = *seq2_align_len;
    if fence_hit.as_deref().copied() != Some(true) {
        while d > 0 {
            let mut new_seq2_index = 0;
            let new_diag = s_get_next_non_affine_tback(aux, d, best_diag, &mut new_seq2_index);
            if new_diag == best_diag {
                if seq2_index - new_seq2_index > 0 {
                    gap_prelim_edit_block_add(
                        edit,
                        EGapAlignOpType::Sub,
                        seq2_index - new_seq2_index,
                    );
                }
            } else if new_diag < best_diag {
                if seq2_index - new_seq2_index > 0 {
                    gap_prelim_edit_block_add(
                        edit,
                        EGapAlignOpType::Sub,
                        seq2_index - new_seq2_index,
                    );
                }
                gap_prelim_edit_block_add(edit, EGapAlignOpType::Ins, 1);
            } else {
                if seq2_index - new_seq2_index - 1 > 0 {
                    gap_prelim_edit_block_add(
                        edit,
                        EGapAlignOpType::Sub,
                        seq2_index - new_seq2_index - 1,
                    );
                }
                gap_prelim_edit_block_add(edit, EGapAlignOpType::Del, 1);
            }
            d -= 1;
            best_diag = new_diag;
            seq2_index = new_seq2_index;
        }
    }
    gap_prelim_edit_block_add(
        edit,
        EGapAlignOpType::Sub,
        aux.last_seq2_off[0].get(diag_origin),
    );
    best_dist
}

pub fn blast_affine_greedy_align(
    seq1: &[u8],
    len1: i32,
    seq2: &[u8],
    len2: i32,
    reverse: bool,
    mut xdrop_threshold: i32,
    mut match_score: i32,
    mut mismatch_score: i32,
    mut in_gap_open: i32,
    mut in_gap_extend: i32,
    seq1_align_len: &mut i32,
    seq2_align_len: &mut i32,
    aux: &mut SGreedyAlignMem,
    mut edit_block: Option<&mut GapPrelimEditBlock>,
    rem: u8,
    mut fence_hit: Option<&mut bool>,
    seed: &mut SGreedySeed,
) -> i32 {
    if match_score % 2 == 1 {
        match_score *= 2;
        mismatch_score *= 2;
        xdrop_threshold *= 2;
        in_gap_open *= 2;
        in_gap_extend *= 2;
    }
    assert!(
        match_score >= 2
            && mismatch_score >= 0
            && in_gap_open >= 0
            && in_gap_extend >= 0
            && xdrop_threshold >= 0
    );
    assert_eq!(
        aux.configured_scores,
        (match_score, mismatch_score, in_gap_open, in_gap_extend)
    );
    if in_gap_open == 0 && in_gap_extend == 0 {
        return blast_greedy_align(
            seq1,
            len1,
            seq2,
            len2,
            reverse,
            xdrop_threshold,
            match_score,
            mismatch_score,
            seq1_align_len,
            seq2_align_len,
            aux,
            edit_block,
            rem,
            fence_hit,
            seed,
        );
    }
    let mut best_dist = 0;
    let mut best_diag = 0;
    let match_score_half = match_score / 2;
    let mut op_cost = match_score + mismatch_score;
    let mut gap_open = in_gap_open;
    let mut gap_extend = in_gap_extend + match_score_half;
    let score_common_factor = blast_gdb3(&mut op_cost, &mut gap_open, &mut gap_extend);
    let gap_open_extend = gap_open + gap_extend;
    let max_penalty = op_cost.max(gap_open_extend);
    let max_dist = aux.max_dist;
    let scaled_max_dist = max_dist * gap_extend;
    let diag_origin = max_dist + 2;
    let xdrop_offset = (xdrop_threshold + match_score_half) / score_common_factor + 1;
    assert!(aux.max_score.len() >= (scaled_max_dist + xdrop_offset + 1) as usize);
    let initial = s_find_first_mismatch(
        seq1,
        seq2,
        len1,
        len2,
        0,
        0,
        fence_hit.as_deref_mut(),
        reverse,
        rem,
    );
    if fence_hit.as_deref().copied() == Some(true) {
        return -1;
    }
    *seq1_align_len = initial;
    *seq2_align_len = initial;
    seed.start_q = 0;
    seed.start_s = 0;
    seed.match_length = initial;
    let mut longest_match_run = initial;
    if initial == len1 || initial == len2 {
        if let Some(edit) = edit_block.as_deref_mut() {
            gap_prelim_edit_block_add(edit, EGapAlignOpType::Sub, initial);
        }
        return initial * match_score;
    }
    if edit_block.is_some() {
        if aux.space.is_none() {
            aux.space = Some(mb_space_new(0));
        } else {
            s_refresh_mb_space(aux.space.as_deref_mut());
        }
    }
    for i in 0..xdrop_offset {
        aux.max_score[i as usize] = 0;
    }
    let upper_base = scaled_max_dist + 1 + max_penalty;
    assert!(aux.diag_bounds.len() >= (2 * upper_base) as usize);
    macro_rules! lo {
        ($distance:expr) => {
            aux.diag_bounds[(($distance) + max_penalty) as usize]
        };
    }
    macro_rules! hi {
        ($distance:expr) => {
            aux.diag_bounds[(upper_base + ($distance) + max_penalty) as usize]
        };
    }
    const INVALID_DIAG: i32 = 100_000_000;
    for i in 0..max_penalty {
        aux.diag_bounds[i as usize] = INVALID_DIAG;
        aux.diag_bounds[(upper_base + i) as usize] = -INVALID_DIAG;
    }
    aux.last_seq2_off_affine[0].set(
        diag_origin,
        SGreedyOffset {
            match_off: initial,
            insert_off: K_INVALID_OFFSET,
            delete_off: K_INVALID_OFFSET,
        },
    );
    aux.max_score[xdrop_offset as usize] = initial * match_score;
    lo!(0) = diag_origin;
    hi!(0) = diag_origin;
    let mut curr_diag_lower = diag_origin - 1;
    let mut curr_diag_upper = diag_origin + 1;
    let mut end1_diag = 0;
    let mut end2_diag = 0;
    let mut num_nonempty_dist = 1;
    let mut converged = false;
    let mut d = 1;
    while d <= scaled_max_dist {
        let mut curr_extent = 0;
        let mut curr_seq2_index = 0;
        let mut curr_diag = 0;
        let tmp_diag_lower = curr_diag_lower;
        let tmp_diag_upper = curr_diag_upper;
        let xdrop_score = aux.max_score[d as usize] + score_common_factor * d - xdrop_threshold;
        let xdrop_score = ((xdrop_score as f64 / match_score_half as f64).ceil() as i32).max(0);
        for k in tmp_diag_lower..=tmp_diag_upper {
            let mut seq2_index = K_INVALID_OFFSET;
            if k + 1 <= hi!(d - gap_open_extend) && k + 1 >= lo!(d - gap_open_extend) {
                seq2_index = aux.last_seq2_off_affine[(d - gap_open_extend) as usize]
                    .get(k + 1)
                    .match_off;
            }
            if k + 1 <= hi!(d - gap_extend)
                && k + 1 >= lo!(d - gap_extend)
                && seq2_index
                    < aux.last_seq2_off_affine[(d - gap_extend) as usize]
                        .get(k + 1)
                        .delete_off
            {
                seq2_index = aux.last_seq2_off_affine[(d - gap_extend) as usize]
                    .get(k + 1)
                    .delete_off;
            }
            aux.last_seq2_off_affine[d as usize].field(k).delete_off =
                if seq2_index == K_INVALID_OFFSET {
                    K_INVALID_OFFSET
                } else {
                    seq2_index + 1
                };
            seq2_index = K_INVALID_OFFSET;
            if k - 1 <= hi!(d - gap_open_extend) && k - 1 >= lo!(d - gap_open_extend) {
                seq2_index = aux.last_seq2_off_affine[(d - gap_open_extend) as usize]
                    .get(k - 1)
                    .match_off;
            }
            if k - 1 <= hi!(d - gap_extend)
                && k - 1 >= lo!(d - gap_extend)
                && seq2_index
                    < aux.last_seq2_off_affine[(d - gap_extend) as usize]
                        .get(k - 1)
                        .insert_off
            {
                seq2_index = aux.last_seq2_off_affine[(d - gap_extend) as usize]
                    .get(k - 1)
                    .insert_off;
            }
            aux.last_seq2_off_affine[d as usize].field(k).insert_off = seq2_index;
            let cell = aux.last_seq2_off_affine[d as usize].get(k);
            seq2_index = cell.insert_off.max(cell.delete_off);
            if k <= hi!(d - op_cost) && k >= lo!(d - op_cost) {
                seq2_index = seq2_index.max(
                    aux.last_seq2_off_affine[(d - op_cost) as usize]
                        .get(k)
                        .match_off
                        + 1,
                );
            }
            let mut seq1_index = seq2_index + k - diag_origin;
            if seq2_index < 0 || seq1_index + seq2_index < xdrop_score {
                if k == curr_diag_lower {
                    curr_diag_lower += 1;
                } else {
                    aux.last_seq2_off_affine[d as usize].field(k).match_off = K_INVALID_OFFSET;
                }
                continue;
            }
            curr_diag_upper = k;
            let index = s_find_first_mismatch(
                seq1,
                seq2,
                len1,
                len2,
                seq1_index,
                seq2_index,
                fence_hit.as_deref_mut(),
                reverse,
                rem,
            );
            if fence_hit.as_deref().copied() == Some(true) {
                return -1;
            }
            if index > longest_match_run {
                seed.start_q = seq1_index;
                seed.start_s = seq2_index;
                seed.match_length = index;
                longest_match_run = index;
            }
            seq1_index += index;
            seq2_index += index;
            aux.last_seq2_off_affine[d as usize].field(k).match_off = seq2_index;
            if seq1_index + seq2_index > curr_extent {
                curr_extent = seq1_index + seq2_index;
                curr_seq2_index = seq2_index;
                curr_diag = k;
            }
            if seq1_index == len1 {
                curr_diag_upper = k;
                end1_diag = k - 1;
            }
            if seq2_index == len2 {
                curr_diag_lower = k;
                end2_diag = k + 1;
            }
        }
        let curr_score = curr_extent * match_score_half - d * score_common_factor;
        if curr_score > aux.max_score[(d - 1 + xdrop_offset) as usize] {
            aux.max_score[(d + xdrop_offset) as usize] = curr_score;
            best_dist = d;
            best_diag = curr_diag;
            *seq2_align_len = curr_seq2_index;
            *seq1_align_len = curr_seq2_index + best_diag - diag_origin;
        } else {
            aux.max_score[(d + xdrop_offset) as usize] =
                aux.max_score[(d - 1 + xdrop_offset) as usize];
        }
        if curr_diag_lower <= curr_diag_upper {
            num_nonempty_dist += 1;
            lo!(d) = curr_diag_lower;
            hi!(d) = curr_diag_upper;
        } else {
            lo!(d) = INVALID_DIAG;
            hi!(d) = -INVALID_DIAG;
        }
        if lo!(d - max_penalty) <= hi!(d - max_penalty) {
            num_nonempty_dist -= 1;
        }
        if num_nonempty_dist == 0 {
            converged = true;
            break;
        }
        d += 1;
        curr_diag_lower =
            (lo!(d - gap_open_extend).min(lo!(d - gap_extend)) - 1).min(lo!(d - op_cost));
        if end2_diag > 0 {
            curr_diag_lower = curr_diag_lower.max(end2_diag);
        }
        curr_diag_upper =
            (hi!(d - gap_open_extend).max(hi!(d - gap_extend)) + 1).max(hi!(d - op_cost));
        if end1_diag > 0 {
            curr_diag_upper = curr_diag_upper.min(end1_diag);
        }
        if d > max_penalty {
            if edit_block.is_none() {
                aux.last_seq2_off_affine[d as usize] =
                    aux.last_seq2_off_affine[(d - max_penalty - 1) as usize];
            } else {
                let count = curr_diag_upper - curr_diag_lower + 1;
                aux.last_seq2_off_affine[d as usize] =
                    match s_get_mb_space(aux.space.as_deref_mut().unwrap(), count) {
                        Some(ptr) => GreedyRow {
                            ptr: ptr.wrapping_offset((-curr_diag_lower) as isize),
                            low: curr_diag_lower,
                            high: curr_diag_upper,
                        },
                        None => GreedyRow::empty(),
                    };
            }
        }
    }
    if !converged {
        return -1;
    }
    if let Some(edit) = edit_block {
        d = best_dist;
        let mut seq2_index = *seq2_align_len;
        let mut state = EGapAlignOpType::Sub;
        let (diag_lower, diag_upper) = aux.diag_bounds.split_at(upper_base as usize);
        while d > 0 {
            if state == EGapAlignOpType::Sub {
                let mut new_seq2_index = 0;
                state = s_get_next_affine_tback_from_match(
                    aux,
                    diag_lower,
                    diag_upper,
                    max_penalty,
                    &mut d,
                    best_diag,
                    op_cost,
                    &mut new_seq2_index,
                );
                assert!(seq2_index > new_seq2_index);
                gap_prelim_edit_block_add(edit, EGapAlignOpType::Sub, seq2_index - new_seq2_index);
                seq2_index = new_seq2_index;
            } else if state == EGapAlignOpType::Ins {
                gap_prelim_edit_block_add(edit, EGapAlignOpType::Ins, 1);
                state = s_get_next_affine_tback_from_indel(
                    aux,
                    diag_lower,
                    diag_upper,
                    max_penalty,
                    &mut d,
                    best_diag,
                    gap_open,
                    gap_extend,
                    EGapAlignOpType::Ins,
                );
                best_diag -= 1;
            } else {
                gap_prelim_edit_block_add(edit, EGapAlignOpType::Del, 1);
                state = s_get_next_affine_tback_from_indel(
                    aux,
                    diag_lower,
                    diag_upper,
                    max_penalty,
                    &mut d,
                    best_diag,
                    gap_open,
                    gap_extend,
                    EGapAlignOpType::Del,
                );
                best_diag += 1;
                seq2_index -= 1;
            }
        }
        gap_prelim_edit_block_add(
            edit,
            EGapAlignOpType::Sub,
            aux.last_seq2_off_affine[0].get(diag_origin).match_off,
        );
    }
    aux.max_score[(best_dist + xdrop_offset) as usize]
}
