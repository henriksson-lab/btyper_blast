// NCBI BLAST+ 2.17.0+, selected blast_stat.c and raw_scoremat.c definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_encoding::{AMINOACID_TO_NCBISTDAA, BLASTAA_SIZE};
use super::blast_stat::{BlastKarlinBlk, BlastScoreBlk, BLAST_SCORE_MIN};
use super::ncbi_math::blast_nint;
#[path = "protein_matrix_tables.rs"]
mod tables;
pub(crate) use tables::BLOSUM62_VALUES;
use tables::{BLOSUM62_PREFS, S_BLOSUM62_PSM};

pub const K_NCBISTDAA: &[u8; 29] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ\0";
pub const BLAST_MATRIX_NOMINAL: i32 = 0;
pub const BLAST_MATRIX_BEST: i32 = 2;

#[derive(Clone, Copy, Debug)]
pub struct SNCBIPackedScoreMatrix<'a> {
    /// Logical symbols without the trailing C NUL; helpers model that NUL.
    pub symbols: &'a [u8],
    pub scores: &'a [i32],
    pub defscore: i32,
}
pub const NCBISM_BLOSUM62: SNCBIPackedScoreMatrix<'static> = SNCBIPackedScoreMatrix {
    symbols: b"ARNDCQEGHILKMFPSTWYVBJZX*",
    scores: &S_BLOSUM62_PSM,
    defscore: -4,
};

pub fn ncbism_get_index(sm: &SNCBIPackedScoreMatrix<'_>, mut aa: i32) -> i32 {
    if aa >= 0 && (aa as usize) < K_NCBISTDAA.len() {
        aa = K_NCBISTDAA[aa as usize] as i32;
    } else if (aa as u8).is_ascii_lowercase() {
        aa = (aa as u8).to_ascii_uppercase() as i32;
    }
    let c = aa as u8;
    let dim = sm
        .symbols
        .iter()
        .position(|s| *s == 0)
        .unwrap_or(sm.symbols.len());
    if c == 0 {
        return dim as i32;
    }
    match sm.symbols[..dim].iter().position(|s| *s == c) {
        Some(i) => i as i32,
        None => -1,
    }
}
pub fn ncbism_get_score(sm: &SNCBIPackedScoreMatrix<'_>, aa1: i32, aa2: i32) -> i32 {
    let i1 = ncbism_get_index(sm, aa1);
    let i2 = ncbism_get_index(sm, aa2);
    if i1 >= 0 && i2 >= 0 {
        sm.scores[i1 as usize
            * sm.symbols
                .iter()
                .position(|s| *s == 0)
                .unwrap_or(sm.symbols.len())
            + i2 as usize]
    } else {
        sm.defscore
    }
}
pub fn s_ncbism_starts_with(string: &[u8], prefix: &[u8]) -> bool {
    let mut i = 0;
    while i < prefix.len() && prefix[i] != 0 {
        if string.get(i).copied().unwrap_or(0).to_ascii_lowercase() != prefix[i] {
            return false;
        }
        i += 1;
    }
    true
}
/// Only the default BLOSUM62 matrix is retained from the original registry.
/// Logical byte slices are treated as C strings, including an embedded NUL.
pub fn ncbism_get_standard_matrix(name: &[u8]) -> Option<&'static SNCBIPackedScoreMatrix<'static>> {
    match name.first().copied().unwrap_or(0) {
        b'B' | b'b' => {
            if !s_ncbism_starts_with(name, b"blosum") {
                return None;
            }
            match name.get(6).copied().unwrap_or(0) {
                b'6' => {
                    let suffix = &name[6..];
                    let len = suffix.iter().position(|v| *v == 0).unwrap_or(suffix.len());
                    if &suffix[..len] == b"62" {
                        Some(&NCBISM_BLOSUM62)
                    } else {
                        None
                    }
                }
                _ => None,
            }
        }
        _ => None,
    }
}
/// The missing score-block name field is supplied explicitly. Custom/file
/// matrices and other compiled-in matrix names are outside this scoped port.
pub fn blast_score_blk_protein_matrix_load(sbp: &mut BlastScoreBlk, matrix_name: &[u8]) -> i16 {
    let Some(psm) = ncbism_get_standard_matrix(matrix_name) else {
        return 1;
    };
    assert_eq!(sbp.alphabet_size as usize, BLASTAA_SIZE);
    assert_eq!(sbp.matrix.len(), BLASTAA_SIZE);
    assert!(sbp.matrix.iter().all(|r| r.len() == BLASTAA_SIZE));
    for i in 0..sbp.alphabet_size as usize {
        for j in 0..sbp.alphabet_size as usize {
            sbp.matrix[i][j] = BLAST_SCORE_MIN;
        }
    }
    for i in 0..sbp.alphabet_size as usize {
        for j in 0..sbp.alphabet_size as usize {
            if i == AMINOACID_TO_NCBISTDAA[b'U' as usize] as usize
                || i == AMINOACID_TO_NCBISTDAA[b'O' as usize] as usize
                || i == AMINOACID_TO_NCBISTDAA[b'-' as usize] as usize
                || j == AMINOACID_TO_NCBISTDAA[b'U' as usize] as usize
                || j == AMINOACID_TO_NCBISTDAA[b'O' as usize] as usize
                || j == AMINOACID_TO_NCBISTDAA[b'-' as usize] as usize
            {
                continue;
            }
            sbp.matrix[i][j] = ncbism_get_score(psm, i as i32, j as i32);
        }
    }
    let x = AMINOACID_TO_NCBISTDAA[b'X' as usize] as usize;
    let u = AMINOACID_TO_NCBISTDAA[b'U' as usize] as usize;
    let o = AMINOACID_TO_NCBISTDAA[b'O' as usize] as usize;
    let c = AMINOACID_TO_NCBISTDAA[b'C' as usize] as usize;
    for i in 0..sbp.alphabet_size as usize {
        sbp.matrix[u][i] = sbp.matrix[c][i];
        sbp.matrix[i][u] = sbp.matrix[i][c];
        sbp.matrix[o][i] = sbp.matrix[x][i];
        sbp.matrix[i][o] = sbp.matrix[i][x];
    }
    0
}

/// Optional output vectors replace the original optional allocated arrays.
/// Unsupported matrix names return zero without changing output arrays.
pub fn blast_get_matrix_values(
    matrix: Option<&[u8]>,
    mut open: Option<&mut Vec<i32>>,
    mut extension: Option<&mut Vec<i32>>,
    mut lambda: Option<&mut Vec<f64>>,
    mut k: Option<&mut Vec<f64>>,
    mut h: Option<&mut Vec<f64>>,
    mut alpha: Option<&mut Vec<f64>>,
    mut beta: Option<&mut Vec<f64>>,
    mut pref_flags: Option<&mut Vec<i32>>,
) -> i16 {
    let Some(matrix) = matrix else {
        return 0;
    };
    let len = matrix.iter().position(|v| *v == 0).unwrap_or(matrix.len());
    if !matrix[..len].eq_ignore_ascii_case(b"BLOSUM62") {
        return 0;
    }
    let n = BLOSUM62_VALUES.len();
    if let Some(out) = open.as_mut() {
        **out = vec![0; n];
    }
    if let Some(out) = extension.as_mut() {
        **out = vec![0; n];
    }
    if let Some(out) = lambda.as_mut() {
        **out = vec![0.0; n];
    }
    if let Some(out) = k.as_mut() {
        **out = vec![0.0; n];
    }
    if let Some(out) = h.as_mut() {
        **out = vec![0.0; n];
    }
    if let Some(out) = alpha.as_mut() {
        **out = vec![0.0; n];
    }
    if let Some(out) = beta.as_mut() {
        **out = vec![0.0; n];
    }
    if let Some(out) = pref_flags.as_mut() {
        **out = vec![0; n];
    }
    for index in 0..n {
        let values = BLOSUM62_VALUES[index];
        if let Some(out) = open.as_mut() {
            out[index] = values[0] as i32;
        }
        if let Some(out) = extension.as_mut() {
            out[index] = values[1] as i32;
        }
        if let Some(out) = lambda.as_mut() {
            out[index] = values[3];
        }
        if let Some(out) = k.as_mut() {
            out[index] = values[4];
        }
        if let Some(out) = h.as_mut() {
            out[index] = values[5];
        }
        if let Some(out) = alpha.as_mut() {
            out[index] = values[6];
        }
        if let Some(out) = beta.as_mut() {
            out[index] = values[7];
        }
        if let Some(out) = pref_flags.as_mut() {
            out[index] = BLOSUM62_PREFS[index];
        }
    }
    n as i16
}
pub fn blast_get_alpha_beta(
    matrix_name: Option<&[u8]>,
    alpha: &mut f64,
    beta: &mut f64,
    gapped: bool,
    gap_open: i32,
    gap_extend: i32,
    kbp_ungapped: &BlastKarlinBlk,
) {
    let (mut open, mut extension, mut prefs, mut alpha_arr, mut beta_arr) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let num_values = blast_get_matrix_values(
        matrix_name,
        Some(&mut open),
        Some(&mut extension),
        None,
        None,
        None,
        Some(&mut alpha_arr),
        Some(&mut beta_arr),
        Some(&mut prefs),
    );
    if gapped {
        if gap_open == 0 && gap_extend == 0 {
            for i in 1..num_values as usize {
                if prefs[i] == BLAST_MATRIX_BEST {
                    *alpha = alpha_arr[i];
                    *beta = beta_arr[i];
                    break;
                }
            }
        } else {
            for i in 1..num_values as usize {
                if open[i] == gap_open && extension[i] == gap_extend {
                    *alpha = alpha_arr[i];
                    *beta = beta_arr[i];
                    break;
                }
            }
        }
    } else if num_values > 0 {
        *alpha = alpha_arr[0];
        *beta = beta_arr[0];
    } else {
        *alpha = kbp_ungapped.lambda / kbp_ungapped.h;
        *beta = 0.0;
    }
}
pub fn blast_get_protein_gap_existence_extend_params(
    matrix_name: Option<&[u8]>,
    gap_existence: &mut i32,
    gap_extension: &mut i32,
) -> i16 {
    let (mut open, mut extension, mut prefs) = (Vec::new(), Vec::new(), Vec::new());
    let num_values = blast_get_matrix_values(
        matrix_name,
        Some(&mut open),
        Some(&mut extension),
        None,
        None,
        None,
        None,
        None,
        Some(&mut prefs),
    );
    if num_values <= 0 {
        return -1;
    }
    for i in 1..num_values as usize {
        if prefs[i] == BLAST_MATRIX_BEST {
            *gap_existence = open[i];
            *gap_extension = extension[i];
            break;
        }
    }
    0
}
pub fn blast_karlin_blk_gapped_load_from_tables(
    kbp: Option<&mut BlastKarlinBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&[u8]>,
    _standard_only: bool,
) -> i16 {
    let Some(matrix) = matrix_name else {
        return -1;
    };
    let len = matrix.iter().position(|v| *v == 0).unwrap_or(matrix.len());
    if !matrix[..len].eq_ignore_ascii_case(b"BLOSUM62") {
        return 1;
    }
    for row in BLOSUM62_VALUES {
        if blast_nint(row[0]) == gap_open as std::os::raw::c_long
            && blast_nint(row[1]) == gap_extend as std::os::raw::c_long
        {
            if let Some(kbp) = kbp {
                kbp.lambda = row[3];
                kbp.k = row[4];
                kbp.log_k = kbp.k.ln();
                kbp.h = row[5];
            }
            return 0;
        }
    }
    2
}
