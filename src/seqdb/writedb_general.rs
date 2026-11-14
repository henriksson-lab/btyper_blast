// NCBI BLAST+ 2.17.0+, writedb_general.hpp. Author: Kevin Bealer.
// See LICENSE-NCBI.

/// Matches C++ `s_DivideRoundUp`.
pub fn s_divide_round_up(value: i32, blocksize: i32) -> i32 {
    (value + (blocksize - 1)) / blocksize
}
