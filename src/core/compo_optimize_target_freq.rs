// NCBI BLAST+ 2.17.0+, optimize_target_freq.c/h.
// Author: E. Michael Gertz. See LICENSE-NCBI.
use super::compo_linear_algebra::*;

pub fn scaled_symmetric_product_a(w: &mut NlmLtriangMatrix, diagonal: &[f64], alphsize: usize) {
    let m = 2 * alphsize - 1;
    for row in 0..m {
        for col in 0..=row {
            w[row][col] = 0.0;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let dd = diagonal[i * alphsize + j];
            w[j][j] += dd;
            if i > 0 {
                w[i + alphsize - 1][j] += dd;
                w[i + alphsize - 1][i + alphsize - 1] += dd;
            }
        }
    }
}
pub fn multiply_by_a(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for i in 0..2 * alphsize - 1 {
            y[i] = 0.0;
        }
    } else if beta != 1.0 {
        for i in 0..2 * alphsize - 1 {
            y[i] *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            y[j] += alpha * x[i * alphsize + j];
        }
    }
    for i in 1..alphsize {
        for j in 0..alphsize {
            y[i + alphsize - 1] += alpha * x[i * alphsize + j];
        }
    }
}
pub fn multiply_by_atranspose(beta: f64, y: &mut [f64], alphsize: usize, alpha: f64, x: &[f64]) {
    if beta == 0.0 {
        for k in 0..alphsize * alphsize {
            y[k] = 0.0;
        }
    } else if beta != 1.0 {
        for k in 0..alphsize * alphsize {
            y[k] *= beta;
        }
    }
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            y[k] += alpha * x[j];
            if i > 0 {
                y[k] += alpha * x[i + alphsize - 1];
            }
        }
    }
}
pub fn residuals_linear_constraints(
    ra: &mut [f64],
    alphsize: usize,
    x: &[f64],
    row_sums: &[f64],
    col_sums: &[f64],
) {
    for i in 0..alphsize {
        ra[i] = col_sums[i];
    }
    for i in 1..alphsize {
        ra[i + alphsize - 1] = row_sums[i];
    }
    multiply_by_a(1.0, ra, alphsize, -1.0, x);
}
pub fn dual_residuals(
    resids_x: &mut [f64],
    alphsize: usize,
    grads: &NlmDenseMatrix,
    z: &[f64],
    constrain_rel_entropy: bool,
) {
    let n = alphsize * alphsize;
    if constrain_rel_entropy {
        let eta = z[2 * alphsize - 1];
        for i in 0..n {
            resids_x[i] = -grads[0][i] + eta * grads[1][i];
        }
    } else {
        for i in 0..n {
            resids_x[i] = -grads[0][i];
        }
    }
    multiply_by_atranspose(1.0, resids_x, alphsize, 1.0, z);
}
pub fn calculate_residuals(
    rnorm: &mut f64,
    resids_x: &mut [f64],
    alphsize: usize,
    resids_z: &mut [f64],
    values: &[f64],
    grads: &NlmDenseMatrix,
    row_sums: &[f64],
    col_sums: &[f64],
    x: &[f64],
    z: &[f64],
    constrain_rel_entropy: bool,
    relative_entropy: f64,
) {
    dual_residuals(resids_x, alphsize, grads, z, constrain_rel_entropy);
    let norm_x = nlm_euclidean_norm(resids_x, alphsize * alphsize);
    residuals_linear_constraints(resids_z, alphsize, x, row_sums, col_sums);
    let norm_z = if constrain_rel_entropy {
        resids_z[2 * alphsize - 1] = relative_entropy - values[1];
        nlm_euclidean_norm(resids_z, 2 * alphsize)
    } else {
        nlm_euclidean_norm(resids_z, 2 * alphsize - 1)
    };
    *rnorm = (norm_x * norm_x + norm_z * norm_z).sqrt();
}

#[derive(Debug)]
pub struct ReNewtonSystem {
    pub alphsize: usize,
    pub constrain_rel_entropy: bool,
    pub w: NlmLtriangMatrix,
    pub dinv: Vec<f64>,
    pub grad_re: Vec<f64>,
}
pub fn re_newton_system_free(system: &mut Option<ReNewtonSystem>) {
    *system = None;
}
pub fn re_newton_system_new(alphsize: usize) -> Option<ReNewtonSystem> {
    Some(ReNewtonSystem {
        alphsize,
        constrain_rel_entropy: true,
        w: nlm_ltriang_matrix_new(2 * alphsize)?,
        dinv: try_zeroed_f64(alphsize * alphsize)?,
        grad_re: try_zeroed_f64(alphsize * alphsize)?,
    })
}
pub fn factor_re_newton_system(
    system: &mut ReNewtonSystem,
    x: &[f64],
    z: &[f64],
    grads: &NlmDenseMatrix,
    constrain_rel_entropy: bool,
    workspace: &mut [f64],
) {
    let alphsize = system.alphsize;
    let n = alphsize * alphsize;
    let m = if constrain_rel_entropy {
        2 * alphsize
    } else {
        2 * alphsize - 1
    };
    system.constrain_rel_entropy = constrain_rel_entropy;
    if constrain_rel_entropy {
        let eta = z[m - 1];
        for i in 0..n {
            system.dinv[i] = x[i] / (1.0 - eta);
        }
    } else {
        system.dinv[..n].copy_from_slice(&x[..n]);
    }
    scaled_symmetric_product_a(&mut system.w, &system.dinv, alphsize);
    if constrain_rel_entropy {
        system.grad_re[..n].copy_from_slice(&grads[1][..n]);
        system.w[m - 1][m - 1] = 0.0;
        for i in 0..n {
            workspace[i] = system.dinv[i] * system.grad_re[i];
            system.w[m - 1][m - 1] += system.grad_re[i] * workspace[i];
        }
        multiply_by_a(0.0, &mut system.w[m - 1], alphsize, 1.0, workspace);
    }
    nlm_factor_ltriang_pos_def(&mut system.w, m);
}
pub fn solve_re_newton_system(
    x: &mut [f64],
    z: &mut [f64],
    system: &ReNewtonSystem,
    workspace: &mut [f64],
) {
    let alphsize = system.alphsize;
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if system.constrain_rel_entropy {
        ma + 1
    } else {
        ma
    };
    for i in 0..n {
        workspace[i] = x[i] * system.dinv[i];
    }
    multiply_by_a(1.0, z, alphsize, -1.0, workspace);
    if system.constrain_rel_entropy {
        for i in 0..n {
            z[m - 1] -= system.grad_re[i] * workspace[i];
        }
    }
    nlm_solve_ltriang_pos_def(z, m, &system.w);
    if system.constrain_rel_entropy {
        for i in 0..n {
            x[i] += system.grad_re[i] * z[m - 1];
        }
    }
    multiply_by_atranspose(1.0, x, alphsize, 1.0, z);
    for i in 0..n {
        x[i] *= system.dinv[i];
    }
}
pub fn evaluate_re_functions(
    values: &mut [f64],
    grads: &mut NlmDenseMatrix,
    alphsize: usize,
    x: &[f64],
    q: &[f64],
    scores: &[f64],
    constrain_rel_entropy: bool,
) {
    values[0] = 0.0;
    values[1] = 0.0;
    for k in 0..alphsize * alphsize {
        let mut temp = (x[k] / q[k]).ln();
        values[0] += x[k] * temp;
        grads[0][k] = temp + 1.0;
        if constrain_rel_entropy {
            temp += scores[k];
            values[1] += x[k] * temp;
            grads[1][k] = temp + 1.0;
        }
    }
}
pub fn compute_scores_from_probs(
    scores: &mut [f64],
    alphsize: usize,
    target_freqs: &[f64],
    row_freqs: &[f64],
    col_freqs: &[f64],
) {
    for i in 0..alphsize {
        for j in 0..alphsize {
            let k = i * alphsize + j;
            scores[k] = (target_freqs[k] / (row_freqs[i] * col_freqs[j])).ln();
        }
    }
}

pub fn blast_optimize_target_frequencies(
    x: &mut [f64],
    alphsize: usize,
    iteration_count: &mut i32,
    q: &[f64],
    row_probs: &[f64],
    col_probs: &[f64],
    constrain_rel_entropy: bool,
    desired_re: f64,
    tolerance: f64,
    max_iterations: i32,
) -> i32 {
    assert!(alphsize > 0);
    let n = alphsize * alphsize;
    let ma = 2 * alphsize - 1;
    let m = if constrain_rel_entropy { ma + 1 } else { ma };
    assert!(
        x.len() >= n && q.len() >= n && row_probs.len() >= alphsize && col_probs.len() >= alphsize
    );
    let allocated = (|| {
        Some((
            re_newton_system_new(alphsize)?,
            try_zeroed_f64(n)?,
            try_zeroed_f64(ma + 1)?,
            try_zeroed_f64(ma + 1)?,
            try_zeroed_f64(n)?,
            try_zeroed_f64(n)?,
            nlm_dense_matrix_new(2, n)?,
        ))
    })();
    let Some((
        mut system,
        mut resids_x,
        mut resids_z,
        mut z,
        mut old_scores,
        mut workspace,
        mut grads,
    )) = allocated
    else {
        *iteration_count = 0;
        return -1;
    };
    let mut values = [0.0; 2];
    compute_scores_from_probs(&mut old_scores, alphsize, q, row_probs, col_probs);
    x[..n].copy_from_slice(&q[..n]);
    let mut its = 0;
    let mut rnorm = f64::NAN;
    while its <= max_iterations {
        evaluate_re_functions(
            &mut values,
            &mut grads,
            alphsize,
            x,
            q,
            &old_scores,
            constrain_rel_entropy,
        );
        calculate_residuals(
            &mut rnorm,
            &mut resids_x,
            alphsize,
            &mut resids_z,
            &values,
            &grads,
            row_probs,
            col_probs,
            x,
            &z,
            constrain_rel_entropy,
            desired_re,
        );
        if !(rnorm > tolerance) {
            break;
        }
        its += 1;
        if its <= max_iterations {
            factor_re_newton_system(
                &mut system,
                x,
                &z,
                &grads,
                constrain_rel_entropy,
                &mut workspace,
            );
            solve_re_newton_system(&mut resids_x, &mut resids_z, &system, &mut workspace);
            let mut alpha = nlm_step_bound(x, n, &resids_x, 1.0 / 0.95);
            alpha *= 0.95;
            nlm_add_vectors(x, n, alpha, &resids_x);
            nlm_add_vectors(&mut z, m, alpha, &resids_z);
        }
    }
    let mut converged = false;
    if its <= max_iterations && rnorm <= tolerance {
        if !constrain_rel_entropy || z[m - 1] < 1.0 {
            converged = true;
        }
    }
    *iteration_count = its;
    if converged {
        0
    } else {
        1
    }
}
