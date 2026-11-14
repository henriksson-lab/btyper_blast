// NCBI BLAST+ 2.17.0+, selected blast_aascan.c and blast_lookup.h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_aalookup::{BlastAaLookupTable, EBoneType, AA_HITS_PER_CELL};
use super::blast_lookup::{compute_table_index, SSeqRange};
use super::blast_nascan::BlastOffsetPair;
use super::masksubj::s_determine_scanning_offsets;

#[derive(Clone, Copy, Debug)]
pub struct AminoAcidSubject<'a> {
    pub sequence: &'a [u8],
    pub seq_ranges: &'a [SSeqRange],
}
pub type AaScanSubject = fn(
    &BlastAaLookupTable,
    &AminoAcidSubject<'_>,
    &mut [BlastOffsetPair],
    i32,
    &mut [i32; 3],
) -> i32;
#[inline]
pub fn compute_table_index_incremental(
    wordsize: i32,
    charsize: i32,
    mask: i32,
    word: &[u8],
    index: i32,
) -> i32 {
    ((index << charsize) | word[wordsize as usize - 1] as i32) & mask
}

pub fn s_blast_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &AminoAcidSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    assert_eq!(lookup.bone_type, EBoneType::Backbone);
    let output = &mut offset_pairs[..array_size as usize];
    let pv = &lookup.pv;
    let backbone = &lookup.thick_backbone;
    let overflow = &lookup.overflow;
    let mut totalhits = 0i32;
    let word_length = lookup.word_length;
    let charsize = lookup.charsize;
    let mask = lookup.mask;
    assert!(mask >= 0);
    while s_determine_scanning_offsets(subject.seq_ranges, word_length, word_length, s_range) {
        let s_first = s_range[1];
        let s_last = s_range[2];
        // Validate this entire scannable span once. Its incoming-byte iterator
        // has exactly one item for each source position s_first..=s_last.
        let sequence = &subject.sequence[s_first as usize..(s_last + word_length) as usize];
        let mut index = compute_table_index(word_length - 1, charsize, sequence);
        let mut s = s_first;
        for &last in &sequence[word_length as usize - 1..] {
            index = ((index << charsize) | last as i32) & mask;
            if pv[index as usize >> 5] & (1u32 << (index & 31)) != 0 {
                let cell = &backbone[index as usize];
                let numhits = cell.num_used;
                assert_ne!(numhits, 0);
                if numhits <= array_size - totalhits {
                    let destination =
                        &mut output[totalhits as usize..(totalhits + numhits) as usize];
                    let source = if numhits <= AA_HITS_PER_CELL as i32 {
                        &cell.payload[..numhits as usize]
                    } else {
                        let cursor = cell.payload[0] as usize;
                        &overflow[cursor..cursor + numhits as usize]
                    };
                    for (pair, &q_off) in destination.iter_mut().zip(source) {
                        *pair = BlastOffsetPair { q_off, s_off: s };
                    }
                    totalhits += numhits;
                } else {
                    s_range[1] = s;
                    return totalhits;
                }
            }
            s += 1;
        }
        s_range[1] = s;
    }
    totalhits
}
pub fn s_blast_small_aa_scan_subject(
    lookup: &BlastAaLookupTable,
    subject: &AminoAcidSubject<'_>,
    offset_pairs: &mut [BlastOffsetPair],
    array_size: i32,
    s_range: &mut [i32; 3],
) -> i32 {
    assert_eq!(lookup.bone_type, EBoneType::Smallbone);
    let output = &mut offset_pairs[..array_size as usize];
    let pv = &lookup.pv;
    let backbone = &lookup.small_backbone;
    let overflow = &lookup.small_overflow;
    let mut totalhits = 0i32;
    let word_length = lookup.word_length;
    let charsize = lookup.charsize;
    let mask = lookup.mask;
    assert!(mask >= 0);
    while s_determine_scanning_offsets(subject.seq_ranges, word_length, word_length, s_range) {
        let s_first = s_range[1];
        let s_last = s_range[2];
        let sequence = &subject.sequence[s_first as usize..(s_last + word_length) as usize];
        let mut index = compute_table_index(word_length - 1, charsize, sequence);
        let mut s = s_first;
        for &last in &sequence[word_length as usize - 1..] {
            index = ((index << charsize) | last as i32) & mask;
            if pv[index as usize >> 5] & (1u32 << (index & 31)) != 0 {
                let cell = &backbone[index as usize];
                let numhits = cell.num_used as i32;
                assert_ne!(numhits, 0);
                if numhits <= array_size - totalhits {
                    let destination =
                        &mut output[totalhits as usize..(totalhits + numhits) as usize];
                    if numhits <= AA_HITS_PER_CELL as i32 {
                        let source = &cell.payload[..numhits as usize * 2];
                        for (pair, entry) in destination.iter_mut().zip(source.chunks_exact(2)) {
                            *pair = BlastOffsetPair {
                                q_off: u16::from_ne_bytes([entry[0], entry[1]]) as i32,
                                s_off: s,
                            };
                        }
                    } else {
                        let cursor = cell.overflow_cursor() as usize;
                        let source = &overflow[cursor..cursor + numhits as usize];
                        for (pair, &q_off) in destination.iter_mut().zip(source) {
                            *pair = BlastOffsetPair {
                                q_off: q_off as i32,
                                s_off: s,
                            };
                        }
                    }
                    totalhits += numhits;
                } else {
                    s_range[1] = s;
                    return totalhits;
                }
            }
            s += 1;
        }
        s_range[1] = s;
    }
    totalhits
}
pub fn blast_choose_protein_scan_subject(lookup: &mut BlastAaLookupTable) {
    lookup.scansub_callback = Some(if lookup.bone_type == EBoneType::Backbone {
        s_blast_aa_scan_subject
    } else {
        s_blast_small_aa_scan_subject
    });
}
