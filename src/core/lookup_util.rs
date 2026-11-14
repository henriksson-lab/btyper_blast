// NCBI BLAST+ 2.17.0+, lookup_util.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

use super::blast_lookup::SSeqRange;

/// Matches C `ilog2` on its nonnegative input domain. Negative signed inputs
/// do not terminate in the original arithmetic-shift loop.
pub fn ilog2(mut x: i64) -> i32 {
    assert!(x >= 0);
    let mut lg = 0;
    if x == 0 {
        return 0;
    }
    loop {
        x >>= 1;
        if x == 0 {
            break;
        }
        lg += 1;
    }
    lg
}

/// Matches C `EstimateNumTableEntries`. Deliberately sums right minus left,
/// without adding one or subtracting the word width.
pub fn estimate_num_table_entries(location: &[SSeqRange], max_off: &mut i32) -> i32 {
    let mut num_entries = 0;
    let mut curr_max = 0;
    for loc in location {
        num_entries += loc.right - loc.left;
        curr_max = curr_max.max(loc.right);
    }
    *max_off = curr_max;
    num_entries
}
