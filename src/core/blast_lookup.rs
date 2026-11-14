// NCBI BLAST+ 2.17.0+, blast_lookup.c and blast_lookup.h.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

/// Inclusive query interval (upstream `SSeqRange`). Linked `BlastSeqLoc`
/// nodes are supplied as an ordered slice, retaining their traversal order.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SSeqRange {
    pub left: i32,
    pub right: i32,
}

/// Each chain retains upstream's `[allocated_slots, hit_count, offsets...]`
/// layout, including the initial eight slots and exact doubling thresholds.
/// Unused allocated slots are initialized to zero in Rust.
pub type ThinBackbone = Vec<Option<Vec<i32>>>;

/// Matches C `ComputeTableIndex`; word size and character width must produce
/// a nonnegative signed 32-bit index, as on the supported nucleotide paths.
pub fn compute_table_index(wordsize: i32, charsize: i32, word: &[u8]) -> i32 {
    let mut index = 0;
    for i in 0..wordsize {
        index = (index << charsize) | i32::from(word[i as usize]);
    }
    index
}

/// Matches C `BlastLookupAddWordHit`.
pub fn blast_lookup_add_word_hit(
    backbone: &mut [Option<Vec<i32>>],
    wordsize: i32,
    charsize: i32,
    seq: &[u8],
    query_offset: i32,
) {
    let index = compute_table_index(wordsize, charsize, seq) as usize;
    if backbone[index].is_none() {
        let mut chain = vec![0; 8];
        chain[0] = 8;
        chain[1] = 0;
        backbone[index] = Some(chain);
    }
    let chain = backbone[index].as_mut().unwrap();
    let mut chain_size = chain[0];
    let hits_in_chain = chain[1];
    if hits_in_chain + 2 == chain_size {
        chain_size *= 2;
        chain.resize(chain_size as usize, 0);
        chain[0] = chain_size;
    }
    let next = chain[1] as usize + 2;
    chain[next] = query_offset;
    chain[1] += 1;
}

/// Matches C `BlastLookupIndexQueryExactMatches`. Sequence is BLASTna (or
/// the caller's alphabet with `charsize` bits); ambiguous residues reset the
/// rolling word boundary. All intervals must lie inside the sequence.
pub fn blast_lookup_index_query_exact_matches(
    backbone: &mut [Option<Vec<i32>>],
    word_length: i32,
    charsize: i32,
    lut_word_length: i32,
    query: &[u8],
    locations: &[SSeqRange],
) {
    assert!(charsize > 0 && charsize < 8);
    assert!(lut_word_length > 0 && lut_word_length <= word_length);
    let invalid_mask = (0xff_u32 << charsize) as u8;
    for loc in locations {
        let from = loc.left;
        let to = loc.right;
        if word_length > to - from + 1 {
            continue;
        }
        assert!(from >= 0 && to >= from && (to as usize) < query.len());
        let mut seq = from as usize;
        let mut word_target = seq + lut_word_length as usize;
        for offset in from..=to {
            if seq >= word_target {
                blast_lookup_add_word_hit(
                    backbone,
                    lut_word_length,
                    charsize,
                    &query[seq - lut_word_length as usize..],
                    offset - lut_word_length,
                );
            }
            if query[seq] & invalid_mask != 0 {
                word_target = seq + lut_word_length as usize + 1;
            }
            seq += 1;
        }
        if seq >= word_target {
            blast_lookup_add_word_hit(
                backbone,
                lut_word_length,
                charsize,
                &query[seq - lut_word_length as usize..],
                to + 1 - lut_word_length,
            );
        }
    }
}
