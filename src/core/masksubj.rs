// NCBI BLAST+ 2.17.0+, masksubj.inl.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_lookup::SSeqRange;

/// Sequence ranges are half-open, matching BLAST_SequenceBlk (unlike the
/// inclusive query lookup intervals using the same source range type).
pub fn s_determine_scanning_offsets(
    seq_ranges: &[SSeqRange],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    assert!(!seq_ranges.is_empty());
    while range[1] > range[2] {
        range[0] += 1;
        if range[0] >= seq_ranges.len() as i32 {
            return false;
        }
        range[1] = seq_ranges[range[0] as usize].left + word_length - lut_word_length;
        range[2] = seq_ranges[range[0] as usize].right - lut_word_length;
    }
    true
}
