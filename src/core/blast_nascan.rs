// NCBI BLAST+ 2.17.0+, blast_nascan.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

use super::blast_nalookup::BlastMBLookupTable;
use super::blast_nalookup::{BlastNaLookupTable, BlastSmallNaLookupTable, NA_HITS_PER_CELL};

/// Matches C `s_BlastMBLookupHasHits`. Compressed PV indexing retains the
/// low five index bits for the bit position, even when the array shift is >5.
#[inline]
pub fn s_blast_mb_lookup_has_hits(lookup: &BlastMBLookupTable, index: i64) -> i32 {
    if lookup.pv_array[(index >> lookup.pv_array_bts) as usize] & (1_u32 << (index & 31)) != 0 {
        1
    } else {
        0
    }
}

/// Matches C `s_BlastMBLookupRetrieve`. The reversed query-position chain
/// order and one-based-to-zero-based conversion are retained.
#[inline]
pub fn s_blast_mb_lookup_retrieve(
    lookup: &BlastMBLookupTable,
    index: i64,
    offset_pairs: &mut [BlastOffsetPair],
    s_off: i32,
) -> i32 {
    let mut i = 0;
    let mut q_off = lookup.hashtable[index as usize];
    while q_off != 0 {
        offset_pairs[i as usize] = BlastOffsetPair {
            q_off: q_off - 1,
            s_off,
        };
        i += 1;
        q_off = lookup.next_pos[q_off as usize];
    }
    i
}

// Original MB_ACCESS_HITS: unlike the small table, equality stops scanning.
macro_rules! mb_access_hits {
    ($lookup:expr, $pairs:expr, $range:expr, $index:expr, $total:ident, $max:expr) => {
        if s_blast_mb_lookup_has_hits($lookup, i64::from($index)) != 0 {
            if $total >= $max {
                return $total;
            }
            $total += s_blast_mb_lookup_retrieve(
                $lookup,
                i64::from($index),
                &mut $pairs[$total as usize..],
                $range[0],
            );
        }
    };
}

/// Matches C `s_MBScanSubject_10_2`, retaining its two explicit phases.
pub fn s_mb_scan_subject_10_2(
    mb_lt: &BlastMBLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(mb_lt.lut_word_length, 10);
    assert_eq!(mb_lt.scan_step, 2);
    let mut s = (scan_range[0] / 4) as usize;
    let mut total_hits = 0;
    max_hits -= mb_lt.longest_chain;
    let mut init_index = 0;
    let mut entry_base_2 = scan_range[0] % 4 == 2;
    if entry_base_2 {
        init_index = i32::from(subject.sequence[s]) << 16
            | i32::from(subject.sequence[s + 1]) << 8
            | i32::from(subject.sequence[s + 2]);
    }
    while scan_range[0] <= scan_range[1] {
        if !entry_base_2 {
            init_index = i32::from(subject.sequence[s]) << 16
                | i32::from(subject.sequence[s + 1]) << 8
                | i32::from(subject.sequence[s + 2]);
            let index = init_index >> 4;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += 2;
        }
        if scan_range[0] > scan_range[1] {
            break;
        }
        let index = init_index & 0xfffff;
        s += 1;
        mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
        scan_range[0] += 2;
        entry_base_2 = false;
    }
    total_hits
}

/// Matches C `s_MBScanSubject_11_1Mod4`, with four explicit alignment phases.
pub fn s_mb_scan_subject_11_1_mod4(
    mb_lt: &BlastMBLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(mb_lt.lut_word_length, 11);
    assert_eq!(mb_lt.scan_step % 4, 1);
    let mut s = (scan_range[0] / 4) as usize;
    let scan_step = mb_lt.scan_step;
    let scan_step_byte = (scan_step / 4) as usize;
    let mut total_hits = 0;
    max_hits -= mb_lt.longest_chain;
    let mut entry = scan_range[0] % 4;
    while scan_range[0] <= scan_range[1] {
        if entry == 0 {
            let mut index = i32::from(subject.sequence[s]) << 16
                | i32::from(subject.sequence[s + 1]) << 8
                | i32::from(subject.sequence[s + 2]);
            index >>= 2;
            s += scan_step_byte;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
        if entry <= 1 {
            if scan_range[0] > scan_range[1] {
                break;
            }
            let mut index = i32::from(subject.sequence[s]) << 16
                | i32::from(subject.sequence[s + 1]) << 8
                | i32::from(subject.sequence[s + 2]);
            index &= 0x3fffff;
            s += scan_step_byte;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
        if entry <= 2 {
            if scan_range[0] > scan_range[1] {
                break;
            }
            let mut index = i32::from(subject.sequence[s]) << 24
                | i32::from(subject.sequence[s + 1]) << 16
                | i32::from(subject.sequence[s + 2]) << 8
                | i32::from(subject.sequence[s + 3]);
            index = (index >> 6) & 0x3fffff;
            s += scan_step_byte;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
        if scan_range[0] > scan_range[1] {
            break;
        }
        let mut index = i32::from(subject.sequence[s]) << 24
            | i32::from(subject.sequence[s + 1]) << 16
            | i32::from(subject.sequence[s + 2]) << 8
            | i32::from(subject.sequence[s + 3]);
        index = (index >> 4) & 0x3fffff;
        s += scan_step_byte + 1;
        mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
        scan_range[0] += scan_step;
        entry = 0;
    }
    total_hits
}

/// Matches C `s_MBScanSubject_11_2Mod4`, retaining parity-dependent shifts.
pub fn s_mb_scan_subject_11_2_mod4(
    mb_lt: &BlastMBLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(mb_lt.lut_word_length, 11);
    assert_eq!(mb_lt.scan_step % 4, 2);
    let mut s = (scan_range[0] / 4) as usize;
    let scan_step = mb_lt.scan_step;
    let scan_step_byte = (scan_step / 4) as usize;
    let mut total_hits = 0;
    max_hits -= mb_lt.longest_chain;
    let (top_shift, bottom_shift) = if scan_range[0] % 2 == 0 {
        (2, 6)
    } else {
        (0, 4)
    };
    let mut entry_base_23 = scan_range[0] % 4 == 2 || scan_range[0] % 4 == 3;
    while scan_range[0] <= scan_range[1] {
        if !entry_base_23 {
            let mut index = i32::from(subject.sequence[s]) << 16
                | i32::from(subject.sequence[s + 1]) << 8
                | i32::from(subject.sequence[s + 2]);
            index = (index >> top_shift) & 0x3fffff;
            s += scan_step_byte;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
        if scan_range[0] > scan_range[1] {
            break;
        }
        let mut index = i32::from(subject.sequence[s]) << 24
            | i32::from(subject.sequence[s + 1]) << 16
            | i32::from(subject.sequence[s + 2]) << 8
            | i32::from(subject.sequence[s + 3]);
        index = (index >> bottom_shift) & 0x3fffff;
        s += scan_step_byte + 1;
        mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
        scan_range[0] += scan_step;
        entry_base_23 = false;
    }
    total_hits
}

/// Matches C `s_MBScanSubject_Any`, preserving separate aligned, width-16,
/// and three-/four-byte fetching paths. Standard presets use widths <=12.
pub fn s_mb_scan_subject_any(
    mb_lt: &BlastMBLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mask = mb_lt.hashsize - 1;
    let lut_word_length = mb_lt.lut_word_length;
    let mut scan_step = mb_lt.scan_step;
    let mut total_hits = 0;
    assert!(matches!(lut_word_length, 9 | 10 | 11 | 12 | 16));
    max_hits -= mb_lt.longest_chain;
    if scan_step % 4 == 0 && subject.mask_type == ESubjectMaskingType::NoSubjMasking {
        let s_end = scan_range[1] / 4;
        let shift = 2 * (12 - lut_word_length);
        let mut s = scan_range[0] / 4;
        scan_step /= 4;
        while s <= s_end {
            let mut index = i64::from(
                i32::from(subject.sequence[s as usize]) << 16
                    | i32::from(subject.sequence[s as usize + 1]) << 8
                    | i32::from(subject.sequence[s as usize + 2]),
            );
            index >>= shift;
            if s_blast_mb_lookup_has_hits(mb_lt, index) != 0 {
                if total_hits >= max_hits {
                    break;
                }
                total_hits += s_blast_mb_lookup_retrieve(
                    mb_lt,
                    index,
                    &mut offset_pairs[total_hits as usize..],
                    s * 4,
                );
            }
            s += scan_step;
        }
        scan_range[0] = s * 4;
    } else if lut_word_length == 16 && scan_step == 1 {
        let mut s = (scan_range[0] / 4) as usize;
        let mut w = i64::from(subject.sequence[s]) << 24
            | i64::from(subject.sequence[s + 1]) << 16
            | i64::from(subject.sequence[s + 2]) << 8
            | i64::from(subject.sequence[s + 3]);
        if scan_range[0] % 4 != 0 {
            w = w.wrapping_shl(8) | i64::from(subject.sequence[s + 4]);
        }
        let mut shift = 2 * (4 - (scan_range[0] % 4));
        while scan_range[0] <= scan_range[1] {
            if scan_range[0] % 4 == 0 {
                s = (scan_range[0] / 4) as usize;
                w = w.wrapping_shl(8) | i64::from(subject.sequence[s + 4]);
                shift = 8;
            }
            let index = (w >> shift) & mask;
            shift -= 2;
            assert!(shift >= 0);
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
    } else if lut_word_length == 16 {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (4 - (scan_range[0] % 4));
            let s = (scan_range[0] / 4) as usize;
            let w = i64::from(subject.sequence[s]) << 32
                | i64::from(subject.sequence[s + 1]) << 24
                | i64::from(subject.sequence[s + 2]) << 16
                | i64::from(subject.sequence[s + 3]) << 8
                | i64::from(subject.sequence[s + 4]);
            let index = (w >> shift) & mask;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
    } else if lut_word_length > 9 {
        assert!(lut_word_length <= 12);
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (16 - (scan_range[0] % 4 + lut_word_length));
            let s = (scan_range[0] / 4) as usize;
            let mut index = i64::from(
                i32::from(subject.sequence[s]) << 24
                    | i32::from(subject.sequence[s + 1]) << 16
                    | i32::from(subject.sequence[s + 2]) << 8
                    | i32::from(subject.sequence[s + 3]),
            );
            index = (index >> shift) & mask;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
    } else {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (12 - (scan_range[0] % 4 + lut_word_length));
            let s = (scan_range[0] / 4) as usize;
            let mut index = i64::from(
                i32::from(subject.sequence[s]) << 16
                    | i32::from(subject.sequence[s + 1]) << 8
                    | i32::from(subject.sequence[s + 2]),
            );
            index = (index >> shift) & mask;
            mb_access_hits!(mb_lt, offset_pairs, scan_range, index, total_hits, max_hits);
            scan_range[0] += scan_step;
        }
    }
    total_hits
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum ESubjectMaskingType {
    NoSubjMasking = 0,
    SoftSubjMasking = 1,
    HardSubjMasking = 2,
}

/// Borrowed view of the sequence-block fields used by ordinary-table scans.
/// The packed byte slice must include any sentinel/padding bytes read by the
/// selected upstream scan kernel (up to a three-byte fetch at the final word).
#[derive(Clone, Copy, Debug)]
pub struct PackedNucleotideSubject<'a> {
    pub sequence: &'a [u8],
    pub mask_type: ESubjectMaskingType,
}

pub type NaScanSubject = fn(
    &BlastNaLookupTable,
    &PackedNucleotideSubject<'_>,
    &mut [BlastOffsetPair],
    i32,
    &mut [i32; 2],
) -> i32;

/// Matches C `s_BlastLookupGetNumHits`, including the presence-vector gate.
#[inline]
pub fn s_blast_lookup_get_num_hits(lookup: &BlastNaLookupTable, index: i32) -> i32 {
    let index = index as usize;
    if lookup.pv[index >> 5] & (1_u32 << (index & 31)) != 0 {
        lookup.thick_backbone[index].num_used
    } else {
        0
    }
}

/// Matches C `s_BlastLookupRetrieve`.
#[inline]
pub fn s_blast_lookup_retrieve(
    lookup: &BlastNaLookupTable,
    index: i32,
    offset_pairs: &mut [BlastOffsetPair],
    s_off: i32,
) {
    let cell = &lookup.thick_backbone[index as usize];
    let num_hits = cell.num_used;
    let lookup_pos = if num_hits <= NA_HITS_PER_CELL as i32 {
        &cell.payload[..]
    } else {
        &lookup.overflow[cell.payload[0] as usize..]
    };
    for i in 0..num_hits {
        offset_pairs[i as usize] = BlastOffsetPair {
            q_off: lookup_pos[i as usize],
            s_off,
        };
    }
}

/// Matches C `s_BlastNaScanSubject_8_4`. Retains byte-aligned scanning and its
/// whole-chain capacity check, which differs from the small-table scanner.
pub fn s_blast_na_scan_subject_8_4(
    lookup: &BlastNaLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(lookup.lut_word_length, 8);
    let mut s = scan_range[0] / 4;
    let s_end = scan_range[1] / 4;
    let mut total_hits = 0;
    while s <= s_end {
        let index = i32::from(subject.sequence[s as usize]) << 8
            | i32::from(subject.sequence[s as usize + 1]);
        let num_hits = s_blast_lookup_get_num_hits(lookup, index);
        if num_hits != 0 {
            if num_hits > max_hits - total_hits {
                break;
            }
            s_blast_lookup_retrieve(
                lookup,
                index,
                &mut offset_pairs[total_hits as usize..],
                s * 4,
            );
            total_hits += num_hits;
        }
        s += 1;
    }
    scan_range[0] = s * 4;
    total_hits
}

/// Matches C `s_BlastNaScanSubject_Any`. Retains the separate aligned,
/// unaligned three-byte, and short-word two-byte loops.
pub fn s_blast_na_scan_subject_any(
    lookup: &BlastNaLookupTable,
    subject: &PackedNucleotideSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    let mask = lookup.mask;
    let mut scan_step = lookup.scan_step;
    let lut_word_length = lookup.lut_word_length;
    let mut total_hits = 0;
    assert!(scan_step > 0);
    if lut_word_length > 5 {
        if scan_step % 4 == 0 && subject.mask_type == ESubjectMaskingType::NoSubjMasking {
            let s_end = scan_range[1] / 4;
            let shift = 2 * (8 - lut_word_length);
            let mut s = scan_range[0] / 4;
            scan_step /= 4;
            while s <= s_end {
                let mut index = i32::from(subject.sequence[s as usize]) << 8
                    | i32::from(subject.sequence[s as usize + 1]);
                index >>= shift;
                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits != 0 {
                    if num_hits > max_hits - total_hits {
                        break;
                    }
                    s_blast_lookup_retrieve(
                        lookup,
                        index,
                        &mut offset_pairs[total_hits as usize..],
                        s * 4,
                    );
                    total_hits += num_hits;
                }
                s += scan_step;
            }
            scan_range[0] = s * 4;
        } else {
            while scan_range[0] <= scan_range[1] {
                let shift = 2 * (12 - (scan_range[0] % 4 + lut_word_length));
                let s = (scan_range[0] / 4) as usize;
                let mut index = i32::from(subject.sequence[s]) << 16
                    | i32::from(subject.sequence[s + 1]) << 8
                    | i32::from(subject.sequence[s + 2]);
                index = (index >> shift) & mask;
                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits != 0 {
                    if num_hits > max_hits - total_hits {
                        break;
                    }
                    s_blast_lookup_retrieve(
                        lookup,
                        index,
                        &mut offset_pairs[total_hits as usize..],
                        scan_range[0],
                    );
                    total_hits += num_hits;
                }
                scan_range[0] += scan_step;
            }
        }
    } else {
        while scan_range[0] <= scan_range[1] {
            let shift = 2 * (8 - (scan_range[0] % 4 + lut_word_length));
            let s = (scan_range[0] / 4) as usize;
            let mut index =
                i32::from(subject.sequence[s]) << 8 | i32::from(subject.sequence[s + 1]);
            index = (index >> shift) & mask;
            let num_hits = s_blast_lookup_get_num_hits(lookup, index);
            if num_hits != 0 {
                if num_hits > max_hits - total_hits {
                    break;
                }
                s_blast_lookup_retrieve(
                    lookup,
                    index,
                    &mut offset_pairs[total_hits as usize..],
                    scan_range[0],
                );
                total_hits += num_hits;
            }
            scan_range[0] += scan_step;
        }
    }
    total_hits
}

/// Matches C `s_NaChooseScanSubject`; typed function pointers replace void*.
pub fn s_na_choose_scan_subject(lookup: &mut BlastNaLookupTable) {
    lookup.scansub_callback = Some(if lookup.lut_word_length == 8 && lookup.scan_step == 4 {
        s_blast_na_scan_subject_8_4
    } else {
        s_blast_na_scan_subject_any
    });
}

/// The nucleotide member of upstream `BlastOffsetPair`'s union.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct BlastOffsetPair {
    pub q_off: i32,
    pub s_off: i32,
}

/// Matches C `s_BlastSmallNaRetrieveHits`. Caller provides space for the
/// complete chain; the order of offsets is the original indexing order.
pub fn s_blast_small_na_retrieve_hits(
    offset_pairs: &mut [BlastOffsetPair],
    mut index: i32,
    s_off: i32,
    total_hits: i32,
    overflow: &[i16],
) -> i32 {
    if index >= 0 {
        offset_pairs[total_hits as usize] = BlastOffsetPair {
            q_off: index,
            s_off,
        };
        1
    } else {
        let mut num_hits = 0;
        let mut src_off = -index;
        index = i32::from(overflow[src_off as usize]);
        src_off += 1;
        loop {
            offset_pairs[(total_hits + num_hits) as usize] = BlastOffsetPair {
                q_off: index,
                s_off,
            };
            num_hits += 1;
            index = i32::from(overflow[src_off as usize]);
            src_off += 1;
            if index < 0 {
                break;
            }
        }
        num_hits
    }
}

// Original SMALL_NA_ACCESS_HITS: stopping occurs only at a nonempty cell,
// with a strict > comparison after reserving the longest complete chain.
macro_rules! small_na_access_hits {
    ($lookup:expr, $pairs:expr, $range:expr, $index:expr, $total:ident, $max:expr, $x:expr) => {
        if $index != -1 {
            if $total > $max {
                $range[0] += $x;
                return $total;
            }
            $total += s_blast_small_na_retrieve_hits(
                $pairs,
                $index,
                $range[0] + $x,
                $total,
                &$lookup.overflow,
            );
        }
    };
}

/// Matches C `s_BlastSmallNaScanSubject_8_4`, retaining its eight-byte unroll
/// and remainder entry points. Subject bytes are packed NCBI2na; scan range
/// is inclusive in residues and must start on a four-base boundary.
pub fn s_blast_small_na_scan_subject_8_4(
    lookup: &BlastSmallNaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(lookup.lut_word_length, 8);
    assert_eq!(lookup.scan_step, 4);
    assert_eq!(scan_range[0] % 4, 0);
    let num_words = (scan_range[1] - scan_range[0]) / 4 + 1;
    let mut s = scan_range[0] / 4;
    let mut total_hits = 0;
    max_hits -= lookup.longest_chain;
    let mut init_index = i32::from(subject[s as usize]);
    let mut entry = match num_words % 8 {
        1 => 7,
        2 => 6,
        3 => 5,
        4 => 4,
        5 => 3,
        6 => 2,
        7 => 1,
        _ => 0,
    };
    s -= entry;
    scan_range[0] -= entry * 4;
    // Signed offsets replace C pointers formed before the allocation at the
    // Duff-style entry. Only valid offsets are dereferenced.
    macro_rules! byte {
        ($phase:expr, $next:expr, $offset:expr) => {
            if entry <= $phase {
                init_index = init_index.wrapping_shl(8) | i32::from(subject[(s + $next) as usize]);
                let index = i32::from(lookup.final_backbone[(init_index & 0xffff) as usize]);
                small_na_access_hits!(
                    lookup,
                    offset_pairs,
                    scan_range,
                    index,
                    total_hits,
                    max_hits,
                    $offset
                );
            }
        };
    }
    while scan_range[0] <= scan_range[1] {
        byte!(0, 1, 0);
        byte!(1, 2, 4);
        byte!(2, 3, 8);
        byte!(3, 4, 12);
        byte!(4, 5, 16);
        byte!(5, 6, 20);
        byte!(6, 7, 24);
        init_index = init_index.wrapping_shl(8) | i32::from(subject[(s + 8) as usize]);
        s += 8;
        let index = i32::from(lookup.final_backbone[(init_index & 0xffff) as usize]);
        small_na_access_hits!(
            lookup,
            offset_pairs,
            scan_range,
            index,
            total_hits,
            max_hits,
            28
        );
        scan_range[0] += 32;
        entry = 0;
    }
    total_hits
}

/// Matches C `s_BlastSmallNaScanSubject_8_1Mod4`, retaining four explicit
/// alignment phases. This is the selected width-8 scanner for default
/// megablast word length 28 (stride 21).
pub fn s_blast_small_na_scan_subject_8_1_mod4(
    lookup: &BlastSmallNaLookupTable,
    subject: &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    mut max_hits: i32,
    scan_range: &mut [i32; 2],
) -> i32 {
    assert_eq!(lookup.lut_word_length, 8);
    assert_eq!(lookup.scan_step % 4, 1);
    let scan_step = lookup.scan_step;
    let scan_step_byte = scan_step / 4;
    let mut s = scan_range[0] / 4;
    let mut total_hits = 0;
    max_hits -= lookup.longest_chain;
    let mut entry = scan_range[0] % 4;
    while scan_range[0] <= scan_range[1] {
        if entry == 0 {
            let mut index =
                i32::from(subject[s as usize]) << 8 | i32::from(subject[s as usize + 1]);
            s += scan_step_byte;
            index = i32::from(lookup.final_backbone[index as usize]);
            small_na_access_hits!(
                lookup,
                offset_pairs,
                scan_range,
                index,
                total_hits,
                max_hits,
                0
            );
            scan_range[0] += scan_step;
        }
        if entry <= 1 {
            if scan_range[0] > scan_range[1] {
                break;
            }
            let mut index = i32::from(subject[s as usize]) << 16
                | i32::from(subject[s as usize + 1]) << 8
                | i32::from(subject[s as usize + 2]);
            s += scan_step_byte;
            index = i32::from(lookup.final_backbone[((index >> 6) & 0xffff) as usize]);
            small_na_access_hits!(
                lookup,
                offset_pairs,
                scan_range,
                index,
                total_hits,
                max_hits,
                0
            );
            scan_range[0] += scan_step;
        }
        if entry <= 2 {
            if scan_range[0] > scan_range[1] {
                break;
            }
            let mut index = i32::from(subject[s as usize]) << 16
                | i32::from(subject[s as usize + 1]) << 8
                | i32::from(subject[s as usize + 2]);
            s += scan_step_byte;
            index = i32::from(lookup.final_backbone[((index >> 4) & 0xffff) as usize]);
            small_na_access_hits!(
                lookup,
                offset_pairs,
                scan_range,
                index,
                total_hits,
                max_hits,
                0
            );
            scan_range[0] += scan_step;
        }
        if scan_range[0] > scan_range[1] {
            break;
        }
        let mut index = i32::from(subject[s as usize]) << 16
            | i32::from(subject[s as usize + 1]) << 8
            | i32::from(subject[s as usize + 2]);
        s += scan_step_byte + 1;
        index = i32::from(lookup.final_backbone[((index >> 2) & 0xffff) as usize]);
        small_na_access_hits!(
            lookup,
            offset_pairs,
            scan_range,
            index,
            total_hits,
            max_hits,
            0
        );
        scan_range[0] += scan_step;
        entry = 0;
    }
    total_hits
}
