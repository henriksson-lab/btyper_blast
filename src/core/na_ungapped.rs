// NCBI BLAST+ 2.17.0+, selected na_ungapped.c and blast_extend.h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_extend::{BlastDiagHash, DIAGHASH_NUM_BUCKETS};
use super::blast_lookup::SSeqRange;
use super::blast_nalookup::{
    BlastMBLookupTable, BlastNaLookupTable, BlastSmallNaLookupTable, NA_HITS_PER_CELL,
};
use super::blast_query_info::{b_search_context_info, BlastQueryInfo};
use super::blast_stat::NucleotideScoreMatrix;

/// Matches C `s_BlastDiagHashRetrieve`. Misses preserve caller output values.
#[inline]
pub fn s_blast_diag_hash_retrieve(
    table: &BlastDiagHash,
    diag: i32,
    level: &mut i32,
    hit_len: &mut i32,
    hit_saved: &mut i32,
) -> i32 {
    let bucket = (diag as u32).wrapping_mul(0x9e370001) % DIAGHASH_NUM_BUCKETS;
    let mut index = table.backbone[bucket as usize];
    while index != 0 {
        let cell = &table.chain[index as usize];
        if cell.diag == diag {
            *level = cell.level();
            *hit_len = cell.hit_len;
            *hit_saved = cell.hit_saved();
            return 1;
        }
        index = cell.next;
    }
    0
}

/// Matches C `s_BlastDiagHashInsert`: update matching cells, reuse the first
/// stale collision, then grow/insert at the chain head only when needed.
#[inline]
pub fn s_blast_diag_hash_insert(
    table: &mut BlastDiagHash,
    diag: i32,
    level: i32,
    len: i32,
    hit_saved: i32,
    s_off: i32,
    window_size: i32,
) -> i32 {
    let bucket = (diag as u32).wrapping_mul(0x9e370001) % DIAGHASH_NUM_BUCKETS;
    let mut index = table.backbone[bucket as usize];
    while index != 0 {
        let cell = &mut table.chain[index as usize];
        if cell.diag == diag {
            cell.set_level(level);
            cell.hit_len = len;
            cell.set_hit_saved(hit_saved);
            return 1;
        } else if s_off - cell.level() > window_size {
            cell.diag = diag;
            cell.set_level(level);
            cell.hit_len = len;
            cell.set_hit_saved(hit_saved);
            return 1;
        }
        index = cell.next;
    }
    if table.occupancy == table.capacity {
        table.capacity = table.capacity.wrapping_mul(2);
        table
            .chain
            .resize(table.capacity as usize, Default::default());
    }
    let cell = &mut table.chain[table.occupancy as usize];
    cell.diag = diag;
    cell.set_level(level);
    cell.hit_len = len;
    cell.set_hit_saved(hit_saved);
    cell.next = table.backbone[bucket as usize];
    table.backbone[bucket as usize] = table.occupancy;
    table.occupancy = table.occupancy.wrapping_add(1);
    1
}

/// Borrowed, typed nucleotide members of the C LookupTableWrap's void* field.
#[derive(Clone, Copy, Debug)]
pub enum NucleotideLookup<'a> {
    Megablast(&'a BlastMBLookupTable),
    Small(&'a BlastSmallNaLookupTable),
    Ordinary(&'a BlastNaLookupTable),
}
pub type TLookupCallback = fn(&NucleotideLookupWrap<'_>, i32, i32) -> bool;
#[derive(Clone, Copy, Debug)]
pub struct NucleotideLookupWrap<'a> {
    pub lut: NucleotideLookup<'a>,
    pub lookup_callback: TLookupCallback,
}

/// Matches C `s_MBLookup`.
pub fn s_mb_lookup(lookup_wrap: &NucleotideLookupWrap<'_>, mut index: i32, mut q_pos: i32) -> bool {
    let NucleotideLookup::Megablast(mb_lt) = lookup_wrap.lut else {
        panic!("wrong lookup type");
    };
    index = (i64::from(index) & (mb_lt.hashsize - 1)) as i32;
    q_pos += 1;
    if mb_lt.pv_array[(index >> mb_lt.pv_array_bts) as usize] & (1_u32 << (index & 31)) == 0 {
        return false;
    }
    let mut q_off = mb_lt.hashtable[index as usize];
    while q_off != 0 {
        if q_off == q_pos {
            return true;
        }
        q_off = mb_lt.next_pos[q_off as usize];
    }
    false
}

/// Matches C `s_SmallNaLookup`, retaining the initial direct equality test.
pub fn s_small_na_lookup(lookup_wrap: &NucleotideLookupWrap<'_>, index: i32, q_pos: i32) -> bool {
    let NucleotideLookup::Small(lookup) = lookup_wrap.lut else {
        panic!("wrong lookup type");
    };
    let mut index = i32::from(lookup.final_backbone[(index & lookup.mask) as usize]);
    if index == q_pos {
        return true;
    }
    if index == -1 || index >= 0 {
        return false;
    }
    let mut src_off = -index;
    index = i32::from(lookup.overflow[src_off as usize]);
    src_off += 1;
    loop {
        if index == q_pos {
            return true;
        }
        index = i32::from(lookup.overflow[src_off as usize]);
        src_off += 1;
        if index < 0 {
            break;
        }
    }
    false
}

/// Matches C `s_NaLookup`.
pub fn s_na_lookup(lookup_wrap: &NucleotideLookupWrap<'_>, mut index: i32, q_pos: i32) -> bool {
    let NucleotideLookup::Ordinary(lookup) = lookup_wrap.lut else {
        panic!("wrong lookup type");
    };
    index &= lookup.mask;
    if lookup.pv[index as usize >> 5] & (1_u32 << (index & 31)) == 0 {
        return false;
    }
    let cell = &lookup.thick_backbone[index as usize];
    let num_hits = cell.num_used;
    let lookup_pos = if num_hits <= NA_HITS_PER_CELL as i32 {
        &cell.payload[..]
    } else {
        &lookup.overflow[cell.payload[0] as usize..]
    };
    for i in 0..num_hits {
        if lookup_pos[i as usize] == q_pos {
            return true;
        }
    }
    false
}

/// Matches C `s_IsSeedMasked`. Retains the width/alignment-specific fetches
/// and the stored lookup callback rather than substituting interval searches.
#[inline]
pub fn s_is_seed_masked(
    lookup_wrap: &NucleotideLookupWrap<'_>,
    subject: &[u8],
    s_off: i32,
    lut_word_length: i32,
    q_pos: i32,
) -> bool {
    let s = (s_off / 4) as usize;
    let shift = 2 * (16 - s_off % 4 - lut_word_length);
    let index = match shift {
        8 | 10 | 12 | 14 => {
            (i32::from(subject[s]) << 24
                | i32::from(subject[s + 1]) << 16
                | i32::from(subject[s + 2]) << 8)
                >> shift
        }
        16 | 18 | 20 | 22 => {
            (i32::from(subject[s]) << 24 | i32::from(subject[s + 1]) << 16) >> shift
        }
        24 => i32::from(subject[s]),
        _ => {
            (i32::from(subject[s]) << 24
                | i32::from(subject[s + 1]) << 16
                | i32::from(subject[s + 2]) << 8
                | i32::from(subject[s + 3]))
                >> shift
        }
    };
    !(lookup_wrap.lookup_callback)(lookup_wrap, index, q_pos)
}

/// Matches C `s_TypeOfWord`. Query bytes are unused upstream and omitted.
/// Locations are tested for presence; the original seed-membership callbacks
/// perform the mask checks. Output offsets and extended length are retained.
pub fn s_type_of_word(
    subject: &[u8],
    q_off: &mut i32,
    s_off: &mut i32,
    locations: &[SSeqRange],
    query_info: &BlastQueryInfo,
    s_range: u32,
    word_length: u32,
    lut_word_length: u32,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    check_double: bool,
    extended: &mut i32,
) -> i32 {
    let mut q_end = (*q_off as u32).wrapping_add(word_length) as i32;
    let mut s_end = (*s_off as u32).wrapping_add(word_length) as i32;
    *extended = 0;
    if word_length == lut_word_length {
        return 1;
    }
    let context = b_search_context_info(q_end, query_info);
    let q_range = query_info.contexts[context as usize].query_offset
        + query_info.contexts[context as usize].query_length;
    if !locations.is_empty() {
        if s_is_seed_masked(
            lookup_wrap,
            subject,
            s_end.wrapping_sub(lut_word_length as i32),
            lut_word_length as i32,
            q_end.wrapping_sub(lut_word_length as i32),
        ) {
            return 0;
        }
        loop {
            if !s_is_seed_masked(lookup_wrap, subject, *s_off, lut_word_length as i32, *q_off) {
                break;
            }
            *s_off += 1;
            *q_off += 1;
        }
    }
    let mut ext_to = word_length.wrapping_sub((q_end - *q_off) as u32) as i32;
    let mut ext_max = ((q_range - q_end) as u32).min(s_range.wrapping_sub(s_end as u32)) as i32;
    if ext_to != 0 || !locations.is_empty() {
        if ext_to > ext_max {
            return 0;
        }
        q_end += ext_to;
        s_end += ext_to;
        let mut s_pos = s_end.wrapping_sub(lut_word_length as i32);
        let mut q_pos = q_end.wrapping_sub(lut_word_length as i32);
        while s_pos > *s_off {
            if s_is_seed_masked(lookup_wrap, subject, s_pos, lut_word_length as i32, q_pos) {
                return 0;
            }
            s_pos = s_pos.wrapping_sub(lut_word_length as i32);
            q_pos = q_pos.wrapping_sub(lut_word_length as i32);
        }
        *extended = ext_to;
    }
    if !check_double {
        return 1;
    }
    ext_to = (ext_to as u32).wrapping_add(word_length) as i32;
    ext_max = ext_max.min(ext_to);
    let mut s_pos = s_end;
    let mut q_pos = q_end;
    while (*extended as u32).wrapping_add(lut_word_length) <= ext_max as u32 {
        if s_is_seed_masked(lookup_wrap, subject, s_pos, lut_word_length as i32, q_pos) {
            break;
        }
        s_pos = s_pos.wrapping_add(lut_word_length as i32);
        q_pos = q_pos.wrapping_add(lut_word_length as i32);
        *extended = (*extended as u32).wrapping_add(lut_word_length) as i32;
    }
    s_pos -= lut_word_length as i32 - 1;
    q_pos -= lut_word_length as i32 - 1;
    while *extended < ext_max {
        if s_is_seed_masked(lookup_wrap, subject, s_pos, lut_word_length as i32, q_pos) {
            return 1;
        }
        *extended += 1;
        s_pos += 1;
        q_pos += 1;
    }
    if ext_max == ext_to {
        2
    } else {
        1
    }
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
#[repr(C)]
pub struct BlastUngappedData {
    pub q_start: i32,
    pub s_start: i32,
    pub length: i32,
    pub score: i32,
}

/// Matches C `s_NuclUngappedExtendExact`. Query is BLASTna including internal
/// context sentinels; subject is packed NCBI2na. X is the negative dropoff.
pub fn s_nucl_ungapped_extend_exact(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    q_off: i32,
    s_off: i32,
    x: i32,
    ungapped_data: &mut BlastUngappedData,
) {
    let mut base = (3 - (s_off % 4)) as i16;
    let q_avail = query.len() as i32 - q_off;
    let s_avail = subject_length - s_off;
    let mut q = q_off;
    let mut q_beg = q;
    let mut q_end = q;
    let mut s = s_off / 4;
    let (start, mut remainder) = if q_off < s_off {
        ((s_off - q_off) / 4, (3 - ((s_off - q_off) % 4)) as i16)
    } else {
        (0, 3_i16)
    };
    let mut score = 0;
    let mut sum = 0;
    while s > start || (s == start && base < remainder) {
        if base == 3 {
            s -= 1;
            base = 0;
        } else {
            base += 1;
        }
        let ch = subject[s as usize];
        q -= 1;
        sum += matrix[query[q as usize] as usize][((ch >> (2 * base)) & 3) as usize];
        if sum > 0 {
            q_beg = q;
            score += sum;
            sum = 0;
        } else if sum < x {
            break;
        }
    }
    ungapped_data.q_start = q_beg;
    ungapped_data.s_start = s_off - (q_off - ungapped_data.q_start);
    let sf = if q_avail < s_avail {
        remainder = (3 - ((s_off + q_avail) % 4)) as i16;
        (s_off + q_avail) / 4
    } else {
        remainder = (3 - (subject_length % 4)) as i16;
        subject_length / 4
    };
    q = q_off;
    s = s_off / 4;
    sum = 0;
    base = (3 - (s_off % 4)) as i16;
    let mut x_current = x;
    while s < sf || (s == sf && base > remainder) {
        let ch = subject[s as usize];
        sum += matrix[query[q as usize] as usize][((ch >> (2 * base)) & 3) as usize];
        q += 1;
        if sum > 0 {
            q_end = q;
            score += sum;
            x_current = (-score).max(x);
            sum = 0;
        } else if sum < x_current {
            break;
        }
        if base == 0 {
            base = 3;
            s += 1;
        } else {
            base -= 1;
        }
    }
    ungapped_data.length = q_end - q_beg;
    ungapped_data.score = score;
}

/// Matches C `s_NuclUngappedExtend`. Retains the approximate four-base loops
/// and exact fallback. Callers must supply a valid full-word seed so that
/// aligning its start to the next four-base boundary stays inside both inputs.
pub fn s_nucl_ungapped_extend(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    q_off: i32,
    s_match_end: i32,
    s_off: i32,
    x: i32,
    ungapped_data: &mut BlastUngappedData,
    score_table: &[i32; 256],
    reduced_cutoff: i32,
) {
    let mut len = (4 - (s_off % 4)) % 4;
    let q_ext = q_off + len;
    let s_ext = s_off + len;
    let mut q = q_ext;
    let mut s = s_ext / 4;
    len = q_ext.min(s_ext) / 4;
    let mut score = 0;
    let mut sum = 0;
    let mut new_q = q;
    for _ in 0..len {
        let s_byte = subject[(s - 1) as usize];
        let q_byte = ((i32::from(query[(q - 4) as usize]) << 6)
            | (i32::from(query[(q - 3) as usize]) << 4)
            | (i32::from(query[(q - 2) as usize]) << 2)
            | i32::from(query[(q - 1) as usize])) as u8;
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q - 4;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        s -= 1;
        q -= 4;
    }
    ungapped_data.q_start = new_q;
    ungapped_data.s_start = s_ext - (q_ext - ungapped_data.q_start);
    q = q_ext;
    s = s_ext / 4;
    len = (query.len() as i32 - q_ext).min(subject_length - s_ext) / 4;
    sum = 0;
    new_q = q;
    for _ in 0..len {
        let s_byte = subject[s as usize];
        let q_byte = ((i32::from(query[q as usize]) << 6)
            | (i32::from(query[q as usize + 1]) << 4)
            | (i32::from(query[q as usize + 2]) << 2)
            | i32::from(query[q as usize + 3])) as u8;
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q + 3;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        s += 1;
        q += 4;
    }
    if score >= reduced_cutoff {
        s_nucl_ungapped_extend_exact(
            query,
            subject,
            subject_length,
            matrix,
            q_off,
            s_off,
            x,
            ungapped_data,
        );
    } else {
        ungapped_data.score = score;
        ungapped_data.length =
            (s_match_end - ungapped_data.s_start).max(new_q - ungapped_data.q_start + 1);
    }
}

use super::blast_extend::{blast_save_initial_hit, BlastDiagTable, BlastInitHitList};
use super::blast_parameters::BlastInitialWordParameters;
use super::blast_program::EBlastProgramType;

pub fn s_blastn_diag_table_extend_initial_hit(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    mut q_off: i32,
    mut s_off: i32,
    query_mask: &[SSeqRange],
    query_info: &BlastQueryInfo,
    s_range: i32,
    word_length: i32,
    lut_word_length: i32,
    lut: &NucleotideLookupWrap<'_>,
    word_params: &BlastInitialWordParameters,
    matrix: &NucleotideScoreMatrix,
    diag_table: &mut BlastDiagTable,
    init_hitlist: &mut BlastInitHitList,
    check_masks: bool,
) -> i32 {
    let window_size = word_params.options.window_size;
    let mut hit_ready = 1;
    let two_hits = window_size > 0;
    let mut off_found = false;
    let mut delta_limit = word_params
        .options
        .scan_range
        .min(window_size - word_length);
    let diag = s_off + diag_table.diag_array_length - q_off;
    let real_diag = (diag & diag_table.diag_mask) as usize;
    let last_hit = diag_table.hit_level_array[real_diag].last_hit();
    let hit_saved = diag_table.hit_level_array[real_diag].flag();
    let mut s_end = s_off + word_length;
    let s_off_pos = s_off + diag_table.offset;
    let mut s_end_pos = s_end + diag_table.offset;
    if s_off_pos < last_hit {
        return 0;
    }
    let mut extended = 0;
    if two_hits && (hit_saved != 0 || s_end_pos > last_hit + window_size) {
        let word_type = s_type_of_word(
            subject,
            &mut q_off,
            &mut s_off,
            query_mask,
            query_info,
            s_range as u32,
            word_length as u32,
            lut_word_length as u32,
            lut,
            true,
            &mut extended,
        );
        if word_type == 0 {
            return 0;
        }
        s_end += extended;
        s_end_pos += extended;
        if word_type == 1 {
            let orig_diag = real_diag as i32 + diag_table.diag_array_length;
            let s_a = s_off_pos + word_length - window_size;
            let s_b = s_end_pos - 2 * word_length;
            if delta_limit < 0 {
                delta_limit = 0;
            }
            for delta in 1..=delta_limit {
                let off_diag = ((orig_diag + delta) & diag_table.diag_mask) as usize;
                let off_s_end = diag_table.hit_level_array[off_diag].last_hit();
                let off_s_l = diag_table.hit_len_array.as_ref().unwrap()[off_diag] as i32;
                if off_s_l != 0 && off_s_end - delta >= s_a && off_s_end - off_s_l <= s_b {
                    off_found = true;
                    break;
                }
                let off_diag = ((orig_diag - delta) & diag_table.diag_mask) as usize;
                let off_s_end = diag_table.hit_level_array[off_diag].last_hit();
                let off_s_l = diag_table.hit_len_array.as_ref().unwrap()[off_diag] as i32;
                if off_s_l != 0 && off_s_end >= s_a && off_s_end - off_s_l + delta <= s_b {
                    off_found = true;
                    break;
                }
            }
            if !off_found {
                hit_ready = 0;
            }
        }
    } else if check_masks {
        if s_type_of_word(
            subject,
            &mut q_off,
            &mut s_off,
            query_mask,
            query_info,
            s_range as u32,
            word_length as u32,
            lut_word_length as u32,
            lut,
            false,
            &mut extended,
        ) == 0
        {
            return 0;
        }
        s_end += extended;
        s_end_pos += extended;
    }
    if hit_ready != 0 {
        if word_params.ungapped_extension {
            let context = b_search_context_info(q_off, query_info);
            let cutoffs = &word_params.cutoffs[context as usize];
            let mut ungapped_data = BlastUngappedData::default();
            if word_params.options.program_number == EBlastProgramType::Blastn
                && (word_params.matrix_only_scoring || word_length < 11)
            {
                s_nucl_ungapped_extend_exact(
                    query,
                    subject,
                    subject_length,
                    matrix,
                    q_off,
                    s_off,
                    -cutoffs.x_dropoff,
                    &mut ungapped_data,
                );
            } else {
                s_nucl_ungapped_extend(
                    query,
                    subject,
                    subject_length,
                    matrix,
                    q_off,
                    s_end,
                    s_off,
                    -cutoffs.x_dropoff,
                    &mut ungapped_data,
                    &word_params.nucl_score_table,
                    cutoffs.reduced_nucl_cutoff_score,
                );
            }
            if off_found || ungapped_data.score >= cutoffs.cutoff_score {
                blast_save_initial_hit(init_hitlist, q_off, s_off, Some(ungapped_data));
                s_end_pos = ungapped_data.length + ungapped_data.s_start + diag_table.offset;
            } else {
                hit_ready = 0;
            }
        } else {
            blast_save_initial_hit(init_hitlist, q_off, s_off, None);
        }
    }
    diag_table.hit_level_array[real_diag].set_last_hit(s_end_pos);
    diag_table.hit_level_array[real_diag].set_flag(hit_ready as u32);
    if two_hits {
        diag_table.hit_len_array.as_mut().unwrap()[real_diag] = if hit_ready != 0 {
            0
        } else {
            (s_end_pos - s_off_pos) as u8
        };
    }
    hit_ready
}

pub fn s_blastn_diag_hash_extend_initial_hit(
    query: &[u8],
    subject: &[u8],
    subject_length: i32,
    mut q_off: i32,
    mut s_off: i32,
    query_mask: &[SSeqRange],
    query_info: &BlastQueryInfo,
    s_range: i32,
    word_length: i32,
    lut_word_length: i32,
    lut: &NucleotideLookupWrap<'_>,
    word_params: &BlastInitialWordParameters,
    matrix: &NucleotideScoreMatrix,
    hash_table: &mut BlastDiagHash,
    init_hitlist: &mut BlastInitHitList,
    check_masks: bool,
) -> i32 {
    let window_size = word_params.options.window_size;
    let mut hit_ready = 1;
    let two_hits = window_size > 0;
    let mut off_found = false;
    let mut delta_limit = word_params
        .options
        .scan_range
        .min(window_size - word_length);
    let diag = s_off - q_off;
    let mut s_end = s_off + word_length;
    let s_off_pos = s_off + hash_table.offset;
    let mut s_end_pos = s_end + hash_table.offset;
    let mut last_hit = 0;
    let mut hit_saved = 0;
    let mut s_l = 0;
    let rc = s_blast_diag_hash_retrieve(hash_table, diag, &mut last_hit, &mut s_l, &mut hit_saved);
    if rc == 0 {
        last_hit = 0;
    }
    if s_off_pos < last_hit {
        return 0;
    }
    let mut extended = 0;
    if two_hits && (hit_saved != 0 || s_end_pos > last_hit + window_size) {
        let word_type = s_type_of_word(
            subject,
            &mut q_off,
            &mut s_off,
            query_mask,
            query_info,
            s_range as u32,
            word_length as u32,
            lut_word_length as u32,
            lut,
            true,
            &mut extended,
        );
        if word_type == 0 {
            return 0;
        }
        s_end += extended;
        s_end_pos += extended;
        if word_type == 1 {
            let s_a = s_off_pos + word_length - window_size;
            let s_b = s_end_pos - 2 * word_length;
            if delta_limit < 0 {
                delta_limit = 0;
            }
            for delta in 1..=delta_limit {
                let mut off_s_end = 0;
                let mut off_s_l = 0;
                let mut off_hit_saved = 0;
                let off_rc = s_blast_diag_hash_retrieve(
                    hash_table,
                    diag + delta,
                    &mut off_s_end,
                    &mut off_s_l,
                    &mut off_hit_saved,
                );
                if off_rc != 0
                    && off_s_l != 0
                    && off_s_end - delta >= s_a
                    && off_s_end - off_s_l <= s_b
                {
                    off_found = true;
                    break;
                }
                let off_rc = s_blast_diag_hash_retrieve(
                    hash_table,
                    diag - delta,
                    &mut off_s_end,
                    &mut off_s_l,
                    &mut off_hit_saved,
                );
                if off_rc != 0
                    && off_s_l != 0
                    && off_s_end >= s_a
                    && off_s_end - off_s_l + delta <= s_b
                {
                    off_found = true;
                    break;
                }
            }
            if !off_found {
                hit_ready = 0;
            }
        }
    } else if check_masks {
        if s_type_of_word(
            subject,
            &mut q_off,
            &mut s_off,
            query_mask,
            query_info,
            s_range as u32,
            word_length as u32,
            lut_word_length as u32,
            lut,
            false,
            &mut extended,
        ) == 0
        {
            return 0;
        }
        s_end += extended;
        s_end_pos += extended;
    }
    if hit_ready != 0 {
        if word_params.ungapped_extension {
            let context = b_search_context_info(q_off, query_info);
            let cutoffs = &word_params.cutoffs[context as usize];
            let mut ungapped_data = BlastUngappedData::default();
            if word_params.options.program_number == EBlastProgramType::Blastn
                && (word_params.matrix_only_scoring || word_length < 11)
            {
                s_nucl_ungapped_extend_exact(
                    query,
                    subject,
                    subject_length,
                    matrix,
                    q_off,
                    s_off,
                    -cutoffs.x_dropoff,
                    &mut ungapped_data,
                );
            } else {
                s_nucl_ungapped_extend(
                    query,
                    subject,
                    subject_length,
                    matrix,
                    q_off,
                    s_end,
                    s_off,
                    -cutoffs.x_dropoff,
                    &mut ungapped_data,
                    &word_params.nucl_score_table,
                    cutoffs.reduced_nucl_cutoff_score,
                );
            }
            if off_found || ungapped_data.score >= cutoffs.cutoff_score {
                blast_save_initial_hit(init_hitlist, q_off, s_off, Some(ungapped_data));
                s_end_pos = ungapped_data.length + ungapped_data.s_start + hash_table.offset;
            } else {
                hit_ready = 0;
            }
        } else {
            blast_save_initial_hit(init_hitlist, q_off, s_off, None);
        }
    }
    s_blast_diag_hash_insert(
        hash_table,
        diag,
        s_end_pos,
        if hit_ready != 0 {
            0
        } else {
            s_end_pos - s_off_pos
        },
        hit_ready,
        s_off_pos,
        window_size + delta_limit + 1,
    );
    hit_ready
}

use super::blast_extend::BlastExtendWord;
use super::blast_nascan::BlastOffsetPair;

/// Endpoint sentinels and the rolling compressed query's three-byte prefix
/// remain explicit, so the original negative-relative reads stay in bounds.
#[derive(Clone, Copy)]
pub struct NucleotideExtensionQuery<'a> {
    pub sequence_start: &'a [u8],
    pub compressed_nuc_seq_start: Option<&'a [u8]>,
}
pub type TNaExtendFunction = fn(
    &[BlastOffsetPair],
    &BlastInitialWordParameters,
    &NucleotideLookupWrap<'_>,
    &NucleotideExtensionQuery<'_>,
    &[u8],
    i32,
    &NucleotideScoreMatrix,
    &BlastQueryInfo,
    &mut BlastExtendWord,
    &mut BlastInitHitList,
    u32,
) -> i32;
pub const S_EXACT_MATCH_EXTEND_LEFT: [u8; 256] = [
    4, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
];
pub const S_EXACT_MATCH_EXTEND_RIGHT: [u8; 256] = [
    4, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];
pub fn s_blast_na_extend_direct(
    offset_pairs: &[BlastOffsetPair],
    word_params: &BlastInitialWordParameters,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    query: &NucleotideExtensionQuery<'_>,
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    query_info: &BlastQueryInfo,
    ewp: &mut BlastExtendWord,
    init_hitlist: &mut BlastInitHitList,
    s_range: u32,
) -> i32 {
    let mut hits_extended = 0;
    let (word_length, check_masks) = match lookup_wrap.lut {
        NucleotideLookup::Megablast(lut) => {
            assert_eq!(lut.word_length, lut.lut_word_length);
            (lut.word_length, !lut.stride)
        }
        NucleotideLookup::Small(lut) => (lut.word_length, true),
        NucleotideLookup::Ordinary(lut) => (lut.word_length, true),
    };
    let lut_word_length = word_length;
    let masked_locations = &[];
    for pair in offset_pairs {
        let q_offset = pair.q_off;
        let s_offset = pair.s_off;
        if let Some(hash) = &mut ewp.hash_table {
            hits_extended += s_blastn_diag_hash_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                hash,
                init_hitlist,
                check_masks,
            );
        } else {
            hits_extended += s_blastn_diag_table_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                ewp.diag_table.as_mut().unwrap(),
                init_hitlist,
                check_masks,
            );
        }
    }
    hits_extended
}
pub fn s_blast_na_extend(
    offset_pairs: &[BlastOffsetPair],
    word_params: &BlastInitialWordParameters,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    query: &NucleotideExtensionQuery<'_>,
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    query_info: &BlastQueryInfo,
    ewp: &mut BlastExtendWord,
    init_hitlist: &mut BlastInitHitList,
    s_range: u32,
) -> i32 {
    let mut hits_extended = 0;
    let (word_length, lut_word_length, masked_locations, check_masks) = match lookup_wrap.lut {
        NucleotideLookup::Megablast(lut) => (
            lut.word_length,
            lut.lut_word_length,
            lut.masked_locations.as_slice(),
            !lut.stride,
        ),
        NucleotideLookup::Ordinary(lut) => (
            lut.word_length,
            lut.lut_word_length,
            lut.masked_locations.as_slice(),
            true,
        ),
        _ => panic!("ordinary or contiguous megablast lookup required"),
    };
    let ext_to = word_length - lut_word_length;
    for pair in offset_pairs {
        let mut s_offset = pair.s_off;
        let mut q_offset = pair.q_off;
        let mut ext_left = 0;
        let mut s_off = s_offset;
        let mut q = q_offset;
        let mut s = s_off / 4;
        while ext_left < ext_to.min(s_offset) {
            s_off -= 1;
            q -= 1;
            if s_off % 4 == 3 {
                s -= 1;
            }
            if ((subject[s as usize] << (2 * (s_off % 4))) >> 6)
                != query.sequence_start[(q + 1) as usize]
            {
                break;
            }
            ext_left += 1;
        }
        if ext_left < ext_to {
            let mut ext_right = 0;
            s_off = s_offset + lut_word_length;
            if (s_off + ext_to - ext_left) as u32 > s_range {
                continue;
            }
            q = q_offset + lut_word_length;
            s = s_off / 4;
            while ext_right < ext_to - ext_left {
                if ((subject[s as usize] << (2 * (s_off % 4))) >> 6)
                    != query.sequence_start[(q + 1) as usize]
                {
                    break;
                }
                s_off += 1;
                q += 1;
                if s_off % 4 == 0 {
                    s += 1;
                }
                ext_right += 1;
            }
            if ext_left + ext_right < ext_to {
                continue;
            }
        }
        q_offset -= ext_left;
        s_offset -= ext_left;
        if let Some(hash) = &mut ewp.hash_table {
            hits_extended += s_blastn_diag_hash_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                hash,
                init_hitlist,
                check_masks,
            );
        } else {
            hits_extended += s_blastn_diag_table_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                ewp.diag_table.as_mut().unwrap(),
                init_hitlist,
                check_masks,
            );
        }
    }
    hits_extended
}
pub fn s_blast_na_extend_aligned(
    offset_pairs: &[BlastOffsetPair],
    word_params: &BlastInitialWordParameters,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    query: &NucleotideExtensionQuery<'_>,
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    query_info: &BlastQueryInfo,
    ewp: &mut BlastExtendWord,
    init_hitlist: &mut BlastInitHitList,
    s_range: u32,
) -> i32 {
    let mut hits_extended = 0;
    let (word_length, lut_word_length, masked_locations, check_masks) = match lookup_wrap.lut {
        NucleotideLookup::Megablast(lut) => (
            lut.word_length,
            lut.lut_word_length,
            lut.masked_locations.as_slice(),
            !lut.stride,
        ),
        NucleotideLookup::Ordinary(lut) => (
            lut.word_length,
            lut.lut_word_length,
            lut.masked_locations.as_slice(),
            true,
        ),
        _ => panic!("ordinary or contiguous megablast lookup required"),
    };
    let ext_to = word_length - lut_word_length;
    for pair in offset_pairs {
        let mut s_offset = pair.s_off;
        let mut q_offset = pair.q_off;
        let mut ext_left = 0;
        let mut ext_max = ext_to.min(s_offset);
        let mut q = q_offset;
        let mut s = s_offset / 4;
        while ext_left < ext_max {
            let byte = subject[(s - 1) as usize];
            if (byte & 3) != query.sequence_start[q as usize] {
                break;
            }
            ext_left += 1;
            if ext_left == ext_max {
                break;
            }
            if ((byte >> 2) & 3) != query.sequence_start[(q - 1) as usize] {
                break;
            }
            ext_left += 1;
            if ext_left == ext_max {
                break;
            }
            if ((byte >> 4) & 3) != query.sequence_start[(q - 2) as usize] {
                break;
            }
            ext_left += 1;
            if ext_left == ext_max {
                break;
            }
            if (byte >> 6) != query.sequence_start[(q - 3) as usize] {
                break;
            }
            s -= 1;
            q -= 4;
            ext_left += 1;
        }
        if ext_left < ext_to {
            let mut ext_right = 0;
            ext_max = ext_to - ext_left;
            if (s_offset + lut_word_length + ext_max) as u32 > s_range {
                continue;
            }
            q = q_offset + lut_word_length;
            s = (s_offset + lut_word_length) / 4;
            while ext_right < ext_max {
                let byte = subject[s as usize];
                if (byte >> 6) != query.sequence_start[(q + 1) as usize] {
                    break;
                }
                ext_right += 1;
                if ext_right == ext_max {
                    break;
                }
                if ((byte >> 4) & 3) != query.sequence_start[(q + 2) as usize] {
                    break;
                }
                ext_right += 1;
                if ext_right == ext_max {
                    break;
                }
                if ((byte >> 2) & 3) != query.sequence_start[(q + 3) as usize] {
                    break;
                }
                ext_right += 1;
                if ext_right == ext_max {
                    break;
                }
                if (byte & 3) != query.sequence_start[(q + 4) as usize] {
                    break;
                }
                s += 1;
                q += 4;
                ext_right += 1;
            }
            if ext_left + ext_right < ext_to {
                continue;
            }
        }
        q_offset -= ext_left;
        s_offset -= ext_left;
        if let Some(hash) = &mut ewp.hash_table {
            hits_extended += s_blastn_diag_hash_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                hash,
                init_hitlist,
                check_masks,
            );
        } else {
            hits_extended += s_blastn_diag_table_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                ewp.diag_table.as_mut().unwrap(),
                init_hitlist,
                check_masks,
            );
        }
    }
    hits_extended
}
pub fn s_blast_small_na_extend_aligned_one_byte(
    offset_pairs: &[BlastOffsetPair],
    word_params: &BlastInitialWordParameters,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    query: &NucleotideExtensionQuery<'_>,
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    query_info: &BlastQueryInfo,
    ewp: &mut BlastExtendWord,
    init_hitlist: &mut BlastInitHitList,
    s_range: u32,
) -> i32 {
    let mut hits_extended = 0;
    let NucleotideLookup::Small(lut) = lookup_wrap.lut else {
        panic!("small lookup required")
    };
    let word_length = lut.word_length;
    let lut_word_length = lut.lut_word_length;
    let masked_locations = lut.masked_locations.as_slice();
    let check_masks = true;
    let q = query
        .compressed_nuc_seq_start
        .expect("rolling compressed query required");
    for pair in offset_pairs {
        let mut s_offset = pair.s_off;
        let mut q_offset = pair.q_off;
        let ext_to = word_length - lut_word_length;
        let mut ext_left = 0;
        let context = b_search_context_info(q_offset, query_info);
        let q_start = query_info.contexts[context as usize].query_offset;
        let q_range = q_start + query_info.contexts[context as usize].query_length;
        if s_offset > 0 && q_offset > 0 {
            let q_byte = q[(q_offset - 4 + 3) as usize];
            let s_byte = subject[(s_offset / 4 - 1) as usize];
            ext_left = S_EXACT_MATCH_EXTEND_LEFT[(q_byte ^ s_byte) as usize] as i32;
            ext_left = ext_left.min(ext_to).min(q_offset - q_start);
        }
        if ext_left < ext_to && q_offset + lut_word_length < query.sequence_start.len() as i32 - 2 {
            let q_byte = q[(q_offset + lut_word_length + 3) as usize];
            let s_byte = subject[((s_offset + lut_word_length) / 4) as usize];
            let mut ext_right = S_EXACT_MATCH_EXTEND_RIGHT[(q_byte ^ s_byte) as usize] as i32;
            // C's nested MIN promotes the first result to unsigned via s_range.
            ext_right = (ext_right as u32)
                .min(s_range.wrapping_sub((s_offset + lut_word_length) as u32))
                .min((q_range - (q_offset + lut_word_length)) as u32)
                as i32;
            if ext_left + ext_right < ext_to {
                continue;
            }
        }
        q_offset -= ext_left;
        s_offset -= ext_left;
        if let Some(hash) = &mut ewp.hash_table {
            hits_extended += s_blastn_diag_hash_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                hash,
                init_hitlist,
                check_masks,
            );
        } else {
            hits_extended += s_blastn_diag_table_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                ewp.diag_table.as_mut().unwrap(),
                init_hitlist,
                check_masks,
            );
        }
    }
    hits_extended
}
pub fn s_blast_small_na_extend(
    offset_pairs: &[BlastOffsetPair],
    word_params: &BlastInitialWordParameters,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    query: &NucleotideExtensionQuery<'_>,
    subject: &[u8],
    subject_length: i32,
    matrix: &NucleotideScoreMatrix,
    query_info: &BlastQueryInfo,
    ewp: &mut BlastExtendWord,
    init_hitlist: &mut BlastInitHitList,
    s_range: u32,
) -> i32 {
    let mut hits_extended = 0;
    let NucleotideLookup::Small(lut) = lookup_wrap.lut else {
        panic!("small lookup required")
    };
    let word_length = lut.word_length;
    let lut_word_length = lut.lut_word_length;
    let masked_locations = lut.masked_locations.as_slice();
    let check_masks = true;
    let q = query
        .compressed_nuc_seq_start
        .expect("rolling compressed query required");
    for pair in offset_pairs {
        let mut s_offset = pair.s_off;
        let mut q_offset = pair.q_off;
        let mut ext_left = 0;
        let mut ext_right = 0;
        let context = b_search_context_info(q_offset, query_info);
        let q_start = query_info.contexts[context as usize].query_offset;
        let q_range = q_start + query_info.contexts[context as usize].query_length;
        let mut ext_max = (word_length - lut_word_length)
            .min(s_offset)
            .min(q_offset - q_start);
        let rsdl = 4 - s_offset % 4;
        s_offset += rsdl;
        q_offset += rsdl;
        ext_max += rsdl;
        let mut s_off = s_offset;
        let mut q_off = q_offset;
        while ext_left < ext_max {
            let q_byte = q[(q_off - 4 + 3) as usize];
            let s_byte = subject[(s_off / 4 - 1) as usize];
            let bases = S_EXACT_MATCH_EXTEND_LEFT[(q_byte ^ s_byte) as usize];
            ext_left += bases as i32;
            if bases < 4 {
                break;
            }
            q_off -= 4;
            s_off -= 4;
        }
        ext_left = ext_left.min(ext_max);
        s_off = s_offset;
        q_off = q_offset;
        ext_max = ((word_length - ext_left) as u32)
            .min(s_range.wrapping_sub(s_off as u32))
            .min((q_range - q_off) as u32) as i32;
        while ext_right < ext_max {
            let q_byte = q[(q_off + 3) as usize];
            let s_byte = subject[(s_off / 4) as usize];
            let bases = S_EXACT_MATCH_EXTEND_RIGHT[(q_byte ^ s_byte) as usize];
            ext_right += bases as i32;
            if bases < 4 {
                break;
            }
            q_off += 4;
            s_off += 4;
        }
        ext_right = ext_right.min(ext_max);
        if ext_left + ext_right < word_length {
            continue;
        }
        q_offset -= ext_left;
        s_offset -= ext_left;
        if let Some(hash) = &mut ewp.hash_table {
            hits_extended += s_blastn_diag_hash_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                hash,
                init_hitlist,
                check_masks,
            );
        } else {
            hits_extended += s_blastn_diag_table_extend_initial_hit(
                &query.sequence_start[1..query.sequence_start.len() - 1],
                subject,
                subject_length,
                q_offset,
                s_offset,
                masked_locations,
                query_info,
                s_range as i32,
                word_length,
                lut_word_length,
                lookup_wrap,
                word_params,
                matrix,
                ewp.diag_table.as_mut().unwrap(),
                init_hitlist,
                check_masks,
            );
        }
    }
    hits_extended
}

/// Return the typed callback in place of writing a void* into the table.
pub fn blast_choose_na_extend(lookup_wrap: &mut NucleotideLookupWrap<'_>) -> TNaExtendFunction {
    match lookup_wrap.lut {
        NucleotideLookup::Megablast(lut) => {
            lookup_wrap.lookup_callback = s_mb_lookup;
            if lut.lut_word_length == lut.word_length {
                s_blast_na_extend_direct
            } else if lut.lut_word_length % 4 == 0 && lut.scan_step % 4 == 0 {
                s_blast_na_extend_aligned
            } else {
                s_blast_na_extend
            }
        }
        NucleotideLookup::Small(lut) => {
            lookup_wrap.lookup_callback = s_small_na_lookup;
            if lut.lut_word_length == lut.word_length {
                s_blast_na_extend_direct
            } else if lut.lut_word_length % 4 == 0
                && lut.scan_step % 4 == 0
                && lut.word_length - lut.lut_word_length <= 4
            {
                s_blast_small_na_extend_aligned_one_byte
            } else {
                s_blast_small_na_extend
            }
        }
        NucleotideLookup::Ordinary(lut) => {
            lookup_wrap.lookup_callback = s_na_lookup;
            if lut.lut_word_length == lut.word_length {
                s_blast_na_extend_direct
            } else if lut.lut_word_length % 4 == 0 && lut.scan_step % 4 == 0 {
                s_blast_na_extend_aligned
            } else {
                s_blast_na_extend
            }
        }
    }
}

use super::blast_diagnostics::{blast_ungapped_stats_update, BlastUngappedStats};
use super::blast_extend::{blast_extend_word_exit, blast_init_hit_list_sort_by_score};
use super::blast_nascan::{ESubjectMaskingType, NaScanSubject, PackedNucleotideSubject};
use super::masksubj::s_determine_scanning_offsets;

pub type SmallNaScanSubject =
    fn(&BlastSmallNaLookupTable, &[u8], &mut [BlastOffsetPair], i32, &mut [i32; 2]) -> i32;
pub type MbScanSubject = fn(
    &BlastMBLookupTable,
    &PackedNucleotideSubject<'_>,
    &mut [BlastOffsetPair],
    i32,
    &mut [i32; 2],
) -> i32;
/// Typed representation of the original void* scan callbacks.
#[derive(Clone, Copy, Debug)]
pub enum NucleotideScanCallback {
    Small(SmallNaScanSubject),
    Ordinary(NaScanSubject),
    Megablast(MbScanSubject),
}

/// The default path has no database masking. Callbacks are supplied explicitly
/// in place of the C table's untyped stored callbacks; their original selection
/// must be performed before calling, with a sufficient output buffer.
pub fn blast_na_word_finder(
    subject: &PackedNucleotideSubject<'_>,
    subject_length: i32,
    seq_ranges: &[SSeqRange],
    query: &NucleotideExtensionQuery<'_>,
    query_info: &BlastQueryInfo,
    lookup_wrap: &NucleotideLookupWrap<'_>,
    matrix: &NucleotideScoreMatrix,
    word_params: &BlastInitialWordParameters,
    ewp: &mut BlastExtendWord,
    offset_pairs: &mut [BlastOffsetPair],
    max_hits: i32,
    init_hitlist: &mut BlastInitHitList,
    ungapped_stats: Option<&mut BlastUngappedStats>,
    scansub: NucleotideScanCallback,
    extend: TNaExtendFunction,
) -> i16 {
    assert_eq!(
        subject.mask_type,
        ESubjectMaskingType::NoSubjMasking,
        "database masking is outside the standard-preset scope"
    );
    assert!(max_hits >= 0 && max_hits as usize <= offset_pairs.len());
    let mut total_hits = 0i32;
    let mut hits_extended = 0i32;
    let (word_length, lut_word_length) = match lookup_wrap.lut {
        NucleotideLookup::Small(t) => (t.word_length, t.lut_word_length),
        NucleotideLookup::Megablast(t) => (t.word_length, t.lut_word_length),
        NucleotideLookup::Ordinary(t) => (t.word_length, t.lut_word_length),
    };
    let mut scan_range = [0, 0, subject_length - lut_word_length];
    while s_determine_scanning_offsets(seq_ranges, word_length, lut_word_length, &mut scan_range) {
        let mut window = [scan_range[1], scan_range[2]];
        let hitsfound = match (scansub, lookup_wrap.lut) {
            (NucleotideScanCallback::Small(scan), NucleotideLookup::Small(t)) => {
                scan(t, subject.sequence, offset_pairs, max_hits, &mut window)
            }
            (NucleotideScanCallback::Ordinary(scan), NucleotideLookup::Ordinary(t)) => {
                scan(t, subject, offset_pairs, max_hits, &mut window)
            }
            (NucleotideScanCallback::Megablast(scan), NucleotideLookup::Megablast(t)) => {
                scan(t, subject, offset_pairs, max_hits, &mut window)
            }
            _ => panic!("scanner and lookup types must agree"),
        };
        scan_range[1] = window[0];
        scan_range[2] = window[1];
        if hitsfound == 0 {
            continue;
        }
        total_hits = total_hits.wrapping_add(hitsfound);
        hits_extended = hits_extended.wrapping_add(extend(
            &offset_pairs[..hitsfound as usize],
            word_params,
            lookup_wrap,
            query,
            subject.sequence,
            subject_length,
            matrix,
            query_info,
            ewp,
            init_hitlist,
            (scan_range[2] + lut_word_length) as u32,
        ));
    }
    blast_extend_word_exit(Some(ewp), subject_length);
    blast_ungapped_stats_update(
        ungapped_stats,
        total_hits,
        hits_extended,
        init_hitlist.total,
    );
    if word_params.ungapped_extension {
        blast_init_hit_list_sort_by_score(init_hitlist);
    }
    0
}
