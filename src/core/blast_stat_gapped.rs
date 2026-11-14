// NCBI BLAST+ 2.17.0+, selected nucleotide gapped-statistic definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_stat::{blast_karlin_blk_copy, BlastKarlinBlk};
use super::ncbi_math::blast_gcd;

pub const BLAST_NUM_STAT_VALUES: usize = 11;
/// Historical C name retained; the current upstream row has eleven columns.
pub type ArrayOf8 = [f64; BLAST_NUM_STAT_VALUES];
pub const BLASTN_VALUES_1_5: [ArrayOf8; 2] = [
    [0.0, 0.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0, 0.0, 0.0, 0.0],
    [3.0, 3.0, 1.39, 0.747, 1.38, 1.00, 0.0, 100.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_1_4: [ArrayOf8; 5] = [
    [
        0.0, 0.0, 1.383, 0.738, 1.36, 1.02, 0.0, 100.0, 0.0, 0.0, 0.0,
    ],
    [1.0, 2.0, 1.36, 0.67, 1.2, 1.1, 0.0, 98.0, 0.0, 0.0, 0.0],
    [0.0, 2.0, 1.26, 0.43, 0.90, 1.4, -1.0, 91.0, 0.0, 0.0, 0.0],
    [2.0, 1.0, 1.35, 0.61, 1.1, 1.2, -1.0, 98.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.22, 0.35, 0.72, 1.7, -3.0, 88.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_2_7: [ArrayOf8; 5] = [
    [0.0, 0.0, 0.69, 0.73, 1.34, 0.515, 0.0, 100.0, 0.0, 0.0, 0.0],
    [2.0, 4.0, 0.68, 0.67, 1.2, 0.55, 0.0, 99.0, 0.0, 0.0, 0.0],
    [0.0, 4.0, 0.63, 0.43, 0.90, 0.7, -1.0, 91.0, 0.0, 0.0, 0.0],
    [4.0, 2.0, 0.675, 0.62, 1.1, 0.6, -1.0, 98.0, 0.0, 0.0, 0.0],
    [2.0, 2.0, 0.61, 0.35, 0.72, 1.7, -3.0, 88.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_1_3: [ArrayOf8; 6] = [
    [
        0.0, 0.0, 1.374, 0.711, 1.31, 1.05, 0.0, 100.0, 0.0, 0.0, 0.0,
    ],
    [2.0, 2.0, 1.37, 0.70, 1.2, 1.1, 0.0, 99.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.35, 0.64, 1.1, 1.2, -1.0, 98.0, 0.0, 0.0, 0.0],
    [0.0, 2.0, 1.25, 0.42, 0.83, 1.5, -2.0, 91.0, 0.0, 0.0, 0.0],
    [2.0, 1.0, 1.34, 0.60, 1.1, 1.2, -1.0, 97.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.21, 0.34, 0.71, 1.7, -2.0, 88.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_2_5: [ArrayOf8; 5] = [
    [0.0, 0.0, 0.675, 0.65, 1.1, 0.6, -1.0, 99.0, 0.0, 0.0, 0.0],
    [2.0, 4.0, 0.67, 0.59, 1.1, 0.6, -1.0, 98.0, 0.0, 0.0, 0.0],
    [0.0, 4.0, 0.62, 0.39, 0.78, 0.8, -2.0, 91.0, 0.0, 0.0, 0.0],
    [4.0, 2.0, 0.67, 0.61, 1.0, 0.65, -2.0, 98.0, 0.0, 0.0, 0.0],
    [2.0, 2.0, 0.56, 0.32, 0.59, 0.95, -4.0, 82.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_1_2: [ArrayOf8; 7] = [
    [0.0, 0.0, 1.28, 0.46, 0.85, 1.5, -2.0, 96.0, 0.0, 0.0, 0.0],
    [2.0, 2.0, 1.33, 0.62, 1.1, 1.2, 0.0, 99.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.30, 0.52, 0.93, 1.4, -2.0, 97.0, 0.0, 0.0, 0.0],
    [0.0, 2.0, 1.19, 0.34, 0.66, 1.8, -3.0, 89.0, 0.0, 0.0, 0.0],
    [3.0, 1.0, 1.32, 0.57, 1.0, 1.3, -1.0, 99.0, 0.0, 0.0, 0.0],
    [2.0, 1.0, 1.29, 0.49, 0.92, 1.4, -1.0, 96.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 1.14, 0.26, 0.52, 2.2, -5.0, 85.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_2_3: [ArrayOf8; 9] = [
    [0.0, 0.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0, 0.0, 0.0, 0.0],
    [4.0, 4.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0, 0.0, 0.0, 0.0],
    [2.0, 4.0, 0.615, 0.37, 0.72, 0.85, -3.0, 97.0, 0.0, 0.0, 0.0],
    [0.0, 4.0, 0.55, 0.21, 0.46, 1.2, -5.0, 87.0, 0.0, 0.0, 0.0],
    [3.0, 3.0, 0.615, 0.37, 0.68, 0.9, -3.0, 97.0, 0.0, 0.0, 0.0],
    [6.0, 2.0, 0.63, 0.42, 0.84, 0.75, -2.0, 99.0, 0.0, 0.0, 0.0],
    [5.0, 2.0, 0.625, 0.41, 0.78, 0.8, -2.0, 99.0, 0.0, 0.0, 0.0],
    [4.0, 2.0, 0.61, 0.35, 0.68, 0.9, -3.0, 96.0, 0.0, 0.0, 0.0],
    [2.0, 2.0, 0.515, 0.14, 0.33, 1.55, -9.0, 81.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_3_4: [ArrayOf8; 6] = [
    [6.0, 3.0, 0.389, 0.25, 0.56, 0.7, -5.0, 95.0, 0.0, 0.0, 0.0],
    [5.0, 3.0, 0.375, 0.21, 0.47, 0.8, -6.0, 92.0, 0.0, 0.0, 0.0],
    [4.0, 3.0, 0.351, 0.14, 0.35, 1.0, -9.0, 86.0, 0.0, 0.0, 0.0],
    [6.0, 2.0, 0.362, 0.16, 0.45, 0.8, -4.0, 88.0, 0.0, 0.0, 0.0],
    [
        5.0, 2.0, 0.330, 0.092, 0.28, 1.2, -13.0, 81.0, 0.0, 0.0, 0.0,
    ],
    [
        4.0, 2.0, 0.281, 0.046, 0.16, 1.8, -23.0, 69.0, 0.0, 0.0, 0.0,
    ],
];
pub const BLASTN_VALUES_4_5: [ArrayOf8; 5] = [
    [0.0, 0.0, 0.22, 0.061, 0.22, 1.0, -15.0, 74.0, 0.0, 0.0, 0.0],
    [6.0, 5.0, 0.28, 0.21, 0.47, 0.6, -7.0, 93.0, 0.0, 0.0, 0.0],
    [5.0, 5.0, 0.27, 0.17, 0.39, 0.7, -9.0, 90.0, 0.0, 0.0, 0.0],
    [4.0, 5.0, 0.25, 0.10, 0.31, 0.8, -10.0, 83.0, 0.0, 0.0, 0.0],
    [3.0, 5.0, 0.23, 0.065, 0.25, 0.9, -11.0, 76.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_1_1: [ArrayOf8; 7] = [
    [3.0, 2.0, 1.09, 0.31, 0.55, 2.0, -2.0, 99.0, 0.0, 0.0, 0.0],
    [2.0, 2.0, 1.07, 0.27, 0.49, 2.2, -3.0, 97.0, 0.0, 0.0, 0.0],
    [1.0, 2.0, 1.02, 0.21, 0.36, 2.8, -6.0, 92.0, 0.0, 0.0, 0.0],
    [0.0, 2.0, 0.80, 0.064, 0.17, 4.8, -16.0, 72.0, 0.0, 0.0, 0.0],
    [4.0, 1.0, 1.08, 0.28, 0.54, 2.0, -2.0, 98.0, 0.0, 0.0, 0.0],
    [3.0, 1.0, 1.06, 0.25, 0.46, 2.3, -4.0, 96.0, 0.0, 0.0, 0.0],
    [2.0, 1.0, 0.99, 0.17, 0.30, 3.3, -10.0, 90.0, 0.0, 0.0, 0.0],
];
pub const BLASTN_VALUES_3_2: [ArrayOf8; 1] = [[
    5.0, 5.0, 0.208, 0.030, 0.072, 2.9, -47.0, 77.0, 0.0, 0.0, 0.0,
]];
pub const BLASTN_VALUES_5_4: [ArrayOf8; 2] = [
    [
        10.0, 6.0, 0.163, 0.068, 0.16, 1.0, -19.0, 85.0, 0.0, 0.0, 0.0,
    ],
    [
        8.0, 6.0, 0.146, 0.039, 0.11, 1.3, -29.0, 76.0, 0.0, 0.0, 0.0,
    ],
];

/// Pointer outputs become offsets into the input rows, without an allocation.
pub fn s_split_array_of8(
    input: Option<&[ArrayOf8]>,
    normal: &mut Option<usize>,
    non_affine: &mut Option<usize>,
    split: &mut bool,
) -> i16 {
    let Some(input) = input else { return -1 };
    assert!(!input.is_empty());
    *normal = None;
    *non_affine = None;
    if input[0][0] == 0.0 && input[0][1] == 0.0 {
        *normal = Some(1);
        *non_affine = Some(0);
        *split = true;
    } else {
        *normal = Some(0);
        *split = false;
    }
    0
}

pub fn s_adjust_gap_parameters_by_gcd(
    normal: Option<&mut [ArrayOf8]>,
    linear: Option<&mut [ArrayOf8]>,
    size: i32,
    gap_existence_max: &mut i32,
    gap_extend_max: &mut i32,
    divisor: i32,
) -> i16 {
    if divisor == 1 {
        return 0;
    }
    if size <= 0 {
        return 1;
    }
    *gap_existence_max = gap_existence_max.wrapping_mul(divisor);
    *gap_extend_max = gap_extend_max.wrapping_mul(divisor);
    if let Some(normal) = normal {
        for row in &mut normal[..size as usize] {
            row[0] *= divisor as f64;
            row[1] *= divisor as f64;
            row[2] /= divisor as f64;
            row[5] /= divisor as f64;
        }
    }
    if let Some(linear) = linear {
        linear[0][0] *= divisor as f64;
        linear[0][1] *= divisor as f64;
        linear[0][2] /= divisor as f64;
        linear[0][5] /= divisor as f64;
    }
    0
}

/// Error messages are returned as owned strings instead of toolkit list nodes.
pub fn s_get_nucl_values_array(
    mut reward: i32,
    mut penalty: i32,
    array_size: &mut i32,
    normal: &mut Option<Vec<ArrayOf8>>,
    non_affine: &mut Option<Vec<ArrayOf8>>,
    gap_open_max: &mut i32,
    gap_extend_max: &mut i32,
    round_down: &mut bool,
    error_return: Option<&mut Vec<String>>,
) -> i16 {
    let divisor = blast_gcd(reward, penalty);
    *round_down = false;
    *array_size = 0;
    *normal = None;
    *non_affine = None;
    // The C implementation divides by zero for (reward,penalty)==(0,0).
    assert_ne!(divisor, 0, "nonzero substitution score divisor required");
    if divisor != 1 {
        reward /= divisor;
        penalty /= divisor;
    }
    let (values, open_max, extend_max, round): (&[ArrayOf8], i32, i32, bool) =
        match (reward, penalty) {
            (1, -5) => (&BLASTN_VALUES_1_5, 3, 3, false),
            (1, -4) => (&BLASTN_VALUES_1_4, 2, 2, false),
            (2, -7) => (&BLASTN_VALUES_2_7, 4, 4, true),
            (1, -3) => (&BLASTN_VALUES_1_3, 2, 2, false),
            (2, -5) => (&BLASTN_VALUES_2_5, 4, 4, true),
            (1, -2) => (&BLASTN_VALUES_1_2, 2, 2, false),
            (2, -3) => (&BLASTN_VALUES_2_3, 6, 4, true),
            (3, -4) => (&BLASTN_VALUES_3_4, 6, 3, true),
            (4, -5) => (&BLASTN_VALUES_4_5, 12, 8, false),
            (1, -1) => (&BLASTN_VALUES_1_1, 4, 2, false),
            (3, -2) => (&BLASTN_VALUES_3_2, 5, 5, false),
            (5, -4) => (&BLASTN_VALUES_5_4, 25, 10, false),
            _ => {
                if let Some(messages) = error_return {
                    messages.push(format!(
                        "Substitution scores {reward} and {penalty} are not supported"
                    ));
                }
                return -1;
            }
        };
    let mut k_values = None;
    let mut k_non_affine = None;
    let mut split = false;
    let status = s_split_array_of8(Some(values), &mut k_values, &mut k_non_affine, &mut split);
    if status != 0 {
        return status;
    }
    *round_down = round;
    *array_size = values.len() as i32;
    *gap_open_max = open_max;
    *gap_extend_max = extend_max;
    if split {
        *array_size -= 1;
    }
    if *array_size > 0 {
        *normal =
            Some(values[k_values.unwrap()..k_values.unwrap() + *array_size as usize].to_vec());
    }
    if let Some(start) = k_non_affine {
        *non_affine = Some(vec![values[start]]);
    }
    s_adjust_gap_parameters_by_gcd(
        normal.as_deref_mut(),
        non_affine.as_deref_mut(),
        *array_size,
        gap_open_max,
        gap_extend_max,
        divisor,
    )
}

pub fn blast_get_nucleotide_gap_existence_extend_params(
    reward: i32,
    penalty: i32,
    gap_existence: &mut i32,
    gap_extension: &mut i32,
) -> i16 {
    let mut array_size = 0;
    let mut normal = None;
    let mut non_affine = None;
    let mut round_down = false;
    let mut gap_existence_max = 0;
    let mut gap_extension_max = 0;
    let status = s_get_nucl_values_array(
        reward,
        penalty,
        &mut array_size,
        &mut normal,
        &mut non_affine,
        &mut gap_existence_max,
        &mut gap_extension_max,
        &mut round_down,
        None,
    );
    if status != 0 {
        return status;
    }
    if *gap_existence == 0 && *gap_extension == 0 && non_affine.is_some() {
        return 0;
    }
    let mut found = false;
    for row in normal
        .as_deref()
        .unwrap_or(&[])
        .iter()
        .take(array_size as usize)
    {
        if *gap_existence as f64 == row[0] && *gap_extension as f64 == row[1] {
            found = true;
            break;
        }
    }
    if !found && (*gap_existence < gap_existence_max || *gap_extension < gap_extension_max) {
        *gap_existence = gap_existence_max;
        *gap_extension = gap_extension_max;
    }
    0
}

pub fn blast_karlin_blk_nucl_gapped_calc(
    kbp: &mut BlastKarlinBlk,
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    kbp_ungap: &BlastKarlinBlk,
    round_down: &mut bool,
    error_return: Option<&mut Vec<String>>,
) -> i16 {
    let mut num_combinations = 0;
    let mut gap_open_max = 0;
    let mut gap_extend_max = 0;
    let mut normal = None;
    let mut linear = None;
    let mut errors = error_return;
    let status = s_get_nucl_values_array(
        reward,
        penalty,
        &mut num_combinations,
        &mut normal,
        &mut linear,
        &mut gap_open_max,
        &mut gap_extend_max,
        round_down,
        errors.as_deref_mut(),
    );
    if status != 0 {
        return status;
    }
    if gap_open == 0 && gap_extend == 0 && linear.is_some() {
        let row = linear.as_ref().unwrap()[0];
        kbp.lambda = row[2];
        kbp.k = row[3];
        kbp.log_k = kbp.k.ln();
        kbp.h = row[4];
    } else {
        let mut index = 0;
        while index < num_combinations {
            let row = normal.as_ref().unwrap()[index as usize];
            if row[0] == gap_open as f64 && row[1] == gap_extend as f64 {
                kbp.lambda = row[2];
                kbp.k = row[3];
                kbp.log_k = kbp.k.ln();
                kbp.h = row[4];
                break;
            }
            index += 1;
        }
        if index == num_combinations {
            if gap_open >= gap_open_max && gap_extend >= gap_extend_max {
                blast_karlin_blk_copy(Some(kbp), Some(kbp_ungap));
            } else if let Some(messages) = errors {
                let mut buffer=format!("Gap existence and extension values {gap_open} and {gap_extend} are not supported for substitution scores {reward} and {penalty}\n");
                for row in normal
                    .as_deref()
                    .unwrap_or(&[])
                    .iter()
                    .take(num_combinations as usize)
                {
                    buffer.push_str(&format!(
                        "{} and {} are supported existence and extension values\n",
                        row[0] as i64, row[1] as i64
                    ));
                }
                buffer.push_str(&format!("{gap_open_max} and {gap_extend_max} are supported existence and extension values\n"));
                buffer.push_str(&format!("Any values more stringent than {gap_open_max} and {gap_extend_max} are supported\n"));
                messages.push(buffer);
                return 1;
            }
        }
    }
    0
}

pub fn s_get_ungapped_beta(reward: i32, penalty: i32) -> f64 {
    let mut beta = 0.0;
    if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        beta = -2.0;
    }
    beta
}

pub fn blast_get_nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    kbp: &BlastKarlinBlk,
    gapped_calculation: bool,
    alpha: &mut f64,
    beta: &mut f64,
) -> i16 {
    let mut num_combinations = 0;
    let mut gap_open_max = 0;
    let mut gap_extend_max = 0;
    let mut normal = None;
    let mut linear = None;
    let mut round_down = false;
    let mut found = false;
    let status = s_get_nucl_values_array(
        reward,
        penalty,
        &mut num_combinations,
        &mut normal,
        &mut linear,
        &mut gap_open_max,
        &mut gap_extend_max,
        &mut round_down,
        None,
    );
    if status != 0 {
        return status;
    }
    if gapped_calculation && normal.is_some() {
        if gap_open == 0 && gap_extend == 0 && linear.is_some() {
            let row = linear.as_ref().unwrap()[0];
            *alpha = row[5];
            *beta = row[6];
            found = true;
        } else {
            for row in normal
                .as_deref()
                .unwrap_or(&[])
                .iter()
                .take(num_combinations as usize)
            {
                if row[0] == gap_open as f64 && row[1] == gap_extend as f64 {
                    *alpha = row[5];
                    *beta = row[6];
                    found = true;
                    break;
                }
            }
        }
    }
    if !found {
        *alpha = kbp.lambda / kbp.h;
        *beta = s_get_ungapped_beta(reward, penalty);
    }
    0
}

pub fn blast_check_reward_penalty_scores(reward: i32, penalty: i32) -> bool {
    let mut array_size = 0;
    let mut normal = None;
    let mut non_affine = None;
    let mut round_down = false;
    let mut gap_existence_max = 0;
    let mut gap_extension_max = 0;
    let status = s_get_nucl_values_array(
        reward,
        penalty,
        &mut array_size,
        &mut normal,
        &mut non_affine,
        &mut gap_existence_max,
        &mut gap_extension_max,
        &mut round_down,
        None,
    );
    status == 0
}
