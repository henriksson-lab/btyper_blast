// NCBI BLAST+ 2.17.0+, selected matrix_frequency_data.c definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::compo_blosum62_data::{BLOSUM62_BG, BLOSUM62_JOINT_PROBS};
#[derive(Debug)]
pub struct CompoFrequencyData {
    pub name: &'static [u8],
    pub joint_probs: &'static [f64; 400],
    pub background: &'static [f64; 20],
}
pub const BLOSUM62_FREQUENCY_DATA: CompoFrequencyData = CompoFrequencyData {
    name: b"BLOSUM62",
    joint_probs: &BLOSUM62_JOINT_PROBS,
    background: &BLOSUM62_BG,
};
pub fn s_locate_frequency_data(matrix: &[u8]) -> Option<&'static CompoFrequencyData> {
    let len = matrix.iter().position(|v| *v == 0).unwrap_or(matrix.len());
    if matrix[..len].eq_ignore_ascii_case(BLOSUM62_FREQUENCY_DATA.name) {
        Some(&BLOSUM62_FREQUENCY_DATA)
    } else {
        None
    }
}
pub fn blast_frequency_data_is_available(matrix_name: &[u8]) -> bool {
    s_locate_frequency_data(matrix_name).is_some()
}
pub fn blast_get_matrix_background_freq(matrix_name: &[u8]) -> Option<&'static [f64; 20]> {
    s_locate_frequency_data(matrix_name).map(|d| d.background)
}
pub fn blast_get_joint_probs_for_matrix(
    probs: &mut [f64],
    row_sums: &mut [f64],
    col_sums: &mut [f64],
    matrix_name: &[u8],
) -> i32 {
    let Some(data) = s_locate_frequency_data(matrix_name) else {
        return -1;
    };
    for j in 0..20 {
        col_sums[j] = 0.0;
    }
    for i in 0..20 {
        row_sums[i] = 0.0;
        for j in 0..20 {
            probs[i * 20 + j] = data.joint_probs[i * 20 + j];
            row_sums[i] += probs[i * 20 + j];
            col_sums[j] += probs[i * 20 + j];
        }
    }
    0
}
