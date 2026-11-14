// NCBI BLAST+ 2.17.0+, selected blast_extend.h storage definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
pub const DIAGHASH_NUM_BUCKETS: u32 = 512;

/// Packed logical representation of C's signed 31-bit level and one-bit flag.
/// Together with the other fields, the original 16-byte cell size is retained.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct DiagHashCell {
    pub diag: i32,
    pub level_and_saved: u32,
    pub hit_len: i32,
    pub next: u32,
}
impl DiagHashCell {
    #[inline]
    pub(crate) fn level(&self) -> i32 {
        (self.level_and_saved.wrapping_shl(1) as i32) >> 1
    }
    #[inline]
    pub(crate) fn hit_saved(&self) -> i32 {
        (self.level_and_saved >> 31) as i32
    }
    #[inline]
    pub(crate) fn set_level(&mut self, value: i32) {
        self.level_and_saved = (self.level_and_saved & 0x80000000) | (value as u32 & 0x7fffffff);
    }
    #[inline]
    pub(crate) fn set_hit_saved(&mut self, value: i32) {
        self.level_and_saved = (self.level_and_saved & 0x7fffffff) | ((value as u32 & 1) << 31);
    }
}

#[derive(Debug, Default)]
pub struct BlastDiagHash {
    pub num_buckets: u32,
    pub occupancy: u32,
    pub capacity: u32,
    pub backbone: Vec<u32>,
    pub chain: Vec<DiagHashCell>,
    pub offset: i32,
    pub window: i32,
}

use super::blast_nascan::BlastOffsetPair;
use super::na_ungapped::BlastUngappedData;
use std::cmp::Ordering;

pub const DIAGHASH_CHAIN_LENGTH: u32 = 256;
pub const MIN_INIT_HITLIST_SIZE: i32 = 100;

#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct DiagStruct {
    pub last_hit_and_flag: u32,
}
impl DiagStruct {
    pub(crate) fn last_hit(&self) -> i32 {
        (self.last_hit_and_flag.wrapping_shl(1) as i32) >> 1
    }
    pub(crate) fn flag(&self) -> u32 {
        self.last_hit_and_flag >> 31
    }
    pub(crate) fn set_last_hit(&mut self, value: i32) {
        self.last_hit_and_flag =
            (self.last_hit_and_flag & 0x80000000) | (value as u32 & 0x7fffffff);
    }
    pub(crate) fn set_flag(&mut self, value: u32) {
        self.last_hit_and_flag = (self.last_hit_and_flag & 0x7fffffff) | ((value & 1) << 31);
    }
}
#[derive(Debug, Default)]
pub struct BlastDiagTable {
    pub hit_level_array: Vec<DiagStruct>,
    pub hit_len_array: Option<Vec<u8>>,
    pub diag_array_length: i32,
    pub diag_mask: i32,
    pub offset: i32,
    pub window: i32,
    pub multiple_hits: bool,
    pub actual_window: i32,
}
#[derive(Debug, Default)]
pub struct BlastExtendWord {
    pub diag_table: Option<BlastDiagTable>,
    pub hash_table: Option<BlastDiagHash>,
}

pub fn s_blast_diag_table_new(qlen: i32, multiple_hits: bool, window_size: i32) -> BlastDiagTable {
    let mut diag_array_length = 1i32;
    while diag_array_length < qlen + window_size {
        diag_array_length = diag_array_length
            .checked_shl(1)
            .filter(|n| *n > 0)
            .expect("diagonal array exceeds signed offset range");
    }
    BlastDiagTable {
        diag_array_length,
        diag_mask: diag_array_length - 1,
        multiple_hits,
        offset: window_size,
        window: window_size,
        ..BlastDiagTable::default()
    }
}
pub fn s_blast_diag_clear(diag: Option<&mut BlastDiagTable>) -> i32 {
    let Some(diag) = diag else { return 0 };
    diag.offset = diag.window;
    for i in 0..diag.diag_array_length as usize {
        diag.hit_level_array[i].set_flag(0);
        diag.hit_level_array[i].set_last_hit(-diag.window);
        if let Some(lengths) = &mut diag.hit_len_array {
            lengths[i] = 0;
        }
    }
    0
}
/// Only the two fields read from BlastInitialWordParameters are passed here.
pub fn blast_extend_word_new(
    query_length: u32,
    use_hash: bool,
    window_size: i32,
) -> BlastExtendWord {
    let mut ewp = BlastExtendWord::default();
    if use_hash {
        ewp.hash_table = Some(BlastDiagHash {
            num_buckets: DIAGHASH_NUM_BUCKETS,
            backbone: vec![0; DIAGHASH_NUM_BUCKETS as usize],
            capacity: DIAGHASH_CHAIN_LENGTH,
            chain: vec![DiagHashCell::default(); DIAGHASH_CHAIN_LENGTH as usize],
            occupancy: 1,
            window: window_size,
            offset: window_size,
        });
    } else {
        let mut diag = s_blast_diag_table_new(query_length as i32, window_size > 0, window_size);
        diag.hit_level_array = vec![DiagStruct::default(); diag.diag_array_length as usize];
        if window_size != 0 {
            diag.hit_len_array = Some(vec![0; diag.diag_array_length as usize]);
        }
        ewp.diag_table = Some(diag);
    }
    ewp
}
pub fn blast_extend_word_exit(ewp: Option<&mut BlastExtendWord>, subject_length: i32) -> i16 {
    let Some(ewp) = ewp else { return -1 };
    if let Some(diag) = &mut ewp.diag_table {
        if diag.offset >= i32::MAX / 4 {
            diag.offset = diag.window;
            s_blast_diag_clear(Some(diag));
        } else {
            diag.offset = diag
                .offset
                .wrapping_add(subject_length)
                .wrapping_add(diag.window);
        }
    } else if let Some(hash) = &mut ewp.hash_table {
        if hash.offset >= i32::MAX / 4 {
            hash.occupancy = 1;
            hash.offset = hash.window;
            hash.backbone[..hash.num_buckets as usize].fill(0);
        } else {
            hash.offset = hash
                .offset
                .wrapping_add(subject_length)
                .wrapping_add(hash.window);
        }
    }
    0
}

#[derive(Debug)]
pub struct BlastInitHSP {
    pub offsets: BlastOffsetPair,
    pub ungapped_data: Option<BlastUngappedData>,
}
#[derive(Debug)]
pub struct BlastInitHitList {
    pub total: i32,
    pub allocated: i32,
    pub init_hsp_array: Vec<BlastInitHSP>,
    pub do_not_reallocate: bool,
}
pub fn blast_init_hit_list_new() -> BlastInitHitList {
    BlastInitHitList {
        total: 0,
        allocated: MIN_INIT_HITLIST_SIZE,
        init_hsp_array: Vec::with_capacity(MIN_INIT_HITLIST_SIZE as usize),
        do_not_reallocate: false,
    }
}
pub fn blast_init_hit_list_reset(list: &mut BlastInitHitList) {
    list.init_hsp_array.clear();
    list.total = 0;
}
pub fn blast_init_hit_list_move(dst: &mut BlastInitHitList, src: &mut BlastInitHitList) {
    assert!(!dst.do_not_reallocate);
    *dst = BlastInitHitList {
        total: src.total,
        allocated: src.allocated,
        init_hsp_array: std::mem::take(&mut src.init_hsp_array),
        do_not_reallocate: src.do_not_reallocate,
    };
    src.total = 0;
    src.allocated = 0;
}
pub fn score_compare_match(h1: &BlastInitHSP, h2: &BlastInitHSP) -> Ordering {
    match (&h1.ungapped_data, &h2.ungapped_data) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Greater,
        (Some(_), None) => Ordering::Less,
        (Some(h1), Some(h2)) => {
            let result = h2
                .score
                .cmp(&h1.score)
                .then(h1.s_start.cmp(&h2.s_start))
                .then(h2.length.cmp(&h1.length))
                .then(h1.q_start.cmp(&h2.q_start));
            if result == Ordering::Equal {
                h2.length.cmp(&h1.length)
            } else {
                result
            }
        }
    }
}
pub fn blast_init_hit_list_sort_by_score(list: &mut BlastInitHitList) {
    list.init_hsp_array[..list.total as usize].sort_by(score_compare_match);
}
pub fn blast_init_hit_list_is_sorted_by_score(list: &BlastInitHitList) -> bool {
    for pair in list.init_hsp_array[..list.total as usize].windows(2) {
        if score_compare_match(&pair[0], &pair[1]) == Ordering::Greater {
            return false;
        }
    }
    true
}
pub fn blast_save_initial_hit(
    list: &mut BlastInitHitList,
    q_off: i32,
    s_off: i32,
    ungapped_data: Option<BlastUngappedData>,
) -> bool {
    let num = list.total;
    let mut num_avail = list.allocated;
    if num >= num_avail {
        if list.do_not_reallocate {
            return false;
        }
        num_avail *= 2;
        // A moved-from list has no storage, as in C; it must be replaced before reuse.
        assert!(num_avail > num);
        list.init_hsp_array
            .reserve_exact(num_avail as usize - list.init_hsp_array.len());
        list.allocated = num_avail;
    }
    list.init_hsp_array.push(BlastInitHSP {
        offsets: BlastOffsetPair { q_off, s_off },
        ungapped_data,
    });
    list.total += 1;
    true
}
pub fn blast_save_init_hsp(
    list: &mut BlastInitHitList,
    q_start: i32,
    s_start: i32,
    q_off: i32,
    s_off: i32,
    len: i32,
    score: i32,
) {
    let data = BlastUngappedData {
        q_start,
        s_start,
        length: len,
        score,
    };
    blast_save_initial_hit(list, q_off, s_off, Some(data));
}
