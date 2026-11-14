// NCBI BLAST+ 2.17.0+, selected Gumbel/Spouge definitions in blast_stat.c/.h.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct BlastGumbelBlk {
    pub lambda: f64,
    pub c: f64,
    pub g: f64,
    pub a: f64,
    pub alpha: f64,
    pub sigma: f64,
    pub a_un: f64,
    pub alpha_un: f64,
    pub b: f64,
    pub beta: f64,
    pub tau: f64,
    pub db_length: i64,
    pub filled: bool,
}

use super::blast_protein_matrix::BLOSUM62_VALUES;
use super::blast_stat::BlastKarlinBlk;
use super::boost_erf::erf_c;
use super::ncbi_math::blast_nint;

/// BLOSUM62 is the supported default protein matrix in this scoped library.
pub fn blast_gumbel_blk_load_from_tables(
    gbp: Option<&mut BlastGumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&[u8]>,
) -> i16 {
    let Some(name) = matrix_name else {
        return -1;
    };
    let name = &name[..name.iter().position(|v| *v == 0).unwrap_or(name.len())];
    if !name.eq_ignore_ascii_case(b"BLOSUM62") {
        return 1;
    }
    for values in &BLOSUM62_VALUES {
        if blast_nint(values[0]) == gap_open as std::os::raw::c_long
            && blast_nint(values[1]) == gap_extend as std::os::raw::c_long
        {
            if let Some(gbp) = gbp {
                gbp.lambda = values[3];
                gbp.c = values[8];
                gbp.g = gap_open.wrapping_add(gap_extend) as f64;
                gbp.a = values[6];
                gbp.alpha = values[9];
                gbp.sigma = values[10];
                gbp.a_un = BLOSUM62_VALUES[0][6];
                gbp.alpha_un = BLOSUM62_VALUES[0][9];
                gbp.b = 2.0 * gbp.g * (gbp.a_un - gbp.a);
                gbp.beta = 2.0 * gbp.g * (gbp.alpha_un - gbp.alpha);
                gbp.tau = 2.0 * gbp.g * (gbp.alpha_un - gbp.sigma);
                gbp.filled = true;
            }
            return 0;
        }
    }
    2
}
pub fn blast_karlin_report_allowed_values(
    matrix_name: Option<&[u8]>,
    error_return: Option<&mut Vec<String>>,
) -> i16 {
    let Some(name) = matrix_name else {
        return 0;
    };
    let name = &name[..name.iter().position(|v| *v == 0).unwrap_or(name.len())];
    if !name.eq_ignore_ascii_case(b"BLOSUM62") {
        return 0;
    }
    if let Some(messages) = error_return {
        for values in &BLOSUM62_VALUES {
            let text = if blast_nint(values[2]) == i16::MAX as std::os::raw::c_long {
                format!(
                    "Gap existence and extension values of {} and {} are supported",
                    blast_nint(values[0]),
                    blast_nint(values[1])
                )
            } else {
                format!("Gap existence, extension and decline-to-align values of {}, {} and {} are supported",blast_nint(values[0]),blast_nint(values[1]),blast_nint(values[2]))
            };
            messages.push(c_error_buffer(text));
        }
    }
    0
}
pub fn blast_gumbel_blk_calc(
    gbp: Option<&mut BlastGumbelBlk>,
    gap_open: i32,
    gap_extend: i32,
    matrix_name: Option<&[u8]>,
    error_return: Option<&mut Vec<String>>,
) -> i16 {
    let status = blast_gumbel_blk_load_from_tables(gbp, gap_open, gap_extend, matrix_name);
    if status != 0 {
        if let Some(messages) = error_return {
            let name = matrix_name.unwrap_or(b"");
            let name = &name[..name.iter().position(|v| *v == 0).unwrap_or(name.len())];
            let name_text = String::from_utf8_lossy(name);
            if status == 1 {
                messages.push(c_error_buffer(format!(
                    "{name_text} is not a supported matrix"
                )));
                // The original loader enumerates these names; only BLOSUM62's
                // numerical data is translated for the supported default.
                for supported in [
                    "BLOSUM80", "BLOSUM62", "BLOSUM50", "BLOSUM45", "PAM250", "BLOSUM90", "PAM30",
                    "PAM70", "IDENTITY",
                ] {
                    messages.push(c_error_buffer(format!("{supported} is a supported matrix")));
                }
            } else if status == 2 {
                messages.push(c_error_buffer(format!("Gap existence and extension values of {gap_open} and {gap_extend} not supported for {name_text}")));
                blast_karlin_report_allowed_values(Some(name), Some(messages));
            }
        }
    }
    status
}
fn c_error_buffer(text: String) -> String {
    let bytes = text.as_bytes();
    String::from_utf8_lossy(&bytes[..bytes.len().min(255)]).into_owned()
}

pub fn blast_spouge_sto_e(
    y: i32,
    kbp: &BlastKarlinBlk,
    gbp: &BlastGumbelBlk,
    m: i32,
    n: i32,
) -> f64 {
    let scale_factor = kbp.lambda / gbp.lambda;
    let db_scale_factor = if gbp.db_length != 0 {
        gbp.db_length as f64 / n as f64
    } else {
        1.0
    };
    let lambda = kbp.lambda;
    let k = kbp.k;
    let ai_hat = gbp.a * scale_factor;
    let bi_hat = gbp.b;
    let alphai_hat = gbp.alpha * scale_factor;
    let betai_hat = gbp.beta;
    let sigma_hat = gbp.sigma * scale_factor;
    let tau_hat = gbp.tau;
    let aj_hat = ai_hat;
    let bj_hat = bi_hat;
    let alphaj_hat = alphai_hat;
    let betaj_hat = betai_hat;
    let const_val = 0.39894228040143267793994605993438;
    let m_li_y = m as f64 - (ai_hat * y as f64 + bi_hat);
    let v1 = 2.0 * alphai_hat / lambda;
    let v2 = alphai_hat * y as f64 + betai_hat;
    let vi_y = if v1 > v2 { v1 } else { v2 };
    let sqrt_vi_y = vi_y.sqrt();
    let m_f = m_li_y / sqrt_vi_y;
    let p_m_f = erf_c(-m_f / 2.0f64.sqrt()) / 2.0;
    let p1 = m_li_y * p_m_f + sqrt_vi_y * const_val * (-0.5 * m_f * m_f).exp();
    let n_lj_y = n as f64 - (aj_hat * y as f64 + bj_hat);
    let v1 = 2.0 * alphaj_hat / lambda;
    let v2 = alphaj_hat * y as f64 + betaj_hat;
    let vj_y = if v1 > v2 { v1 } else { v2 };
    let sqrt_vj_y = vj_y.sqrt();
    let n_f = n_lj_y / sqrt_vj_y;
    let p_n_f = erf_c(-n_f / 2.0f64.sqrt()) / 2.0;
    let p2 = n_lj_y * p_n_f + sqrt_vj_y * const_val * (-0.5 * n_f * n_f).exp();
    let c1 = 2.0 * sigma_hat / lambda;
    let c2 = sigma_hat * y as f64 + tau_hat;
    let c_y = if c1 > c2 { c1 } else { c2 };
    let area = p1 * p2 + c_y * p_m_f * p_n_f;
    let e_value = area * k * (-lambda * y as f64).exp() * db_scale_factor;
    assert!(e_value >= 0.0);
    e_value
}
pub fn blast_spouge_eto_s(
    e0: f64,
    kbp: &BlastKarlinBlk,
    gbp: &BlastGumbelBlk,
    m: i32,
    n: i32,
) -> i32 {
    let mut a = 0;
    let db_scale_factor = if gbp.db_length != 0 {
        gbp.db_length as f64
    } else {
        1.0
    };
    let mut b = ((db_scale_factor / e0).ln() / kbp.lambda) as i32;
    b = b.max(2);
    let mut e = blast_spouge_sto_e(b, kbp, gbp, m, n);
    if e > e0 {
        while e > e0 {
            a = b;
            b = b.wrapping_mul(2);
            e = blast_spouge_sto_e(b, kbp, gbp, m, n);
        }
    } else {
        a = 0;
    }
    while b.wrapping_sub(a) > 1 {
        let c = a.wrapping_add(b) / 2;
        e = blast_spouge_sto_e(c, kbp, gbp, m, n);
        if e > e0 {
            a = c;
        } else {
            b = c;
        }
    }
    a
}
