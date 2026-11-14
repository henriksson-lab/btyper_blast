// NCBI BLAST+ 2.17.0+, selected nlm_linear_algebra.c/h definitions.
// Author: E. Michael Gertz. See LICENSE-NCBI.
use std::ops::{Index, IndexMut};

/// Original contiguous dense storage with explicit row dimensions.
#[derive(Clone, Debug)]
pub struct NlmDenseMatrix {
    pub nrows: usize,
    pub ncols: usize,
    pub data: Vec<f64>,
}
impl Index<usize> for NlmDenseMatrix {
    type Output = [f64];
    fn index(&self, row: usize) -> &Self::Output {
        &self.data[row * self.ncols..(row + 1) * self.ncols]
    }
}
impl IndexMut<usize> for NlmDenseMatrix {
    fn index_mut(&mut self, row: usize) -> &mut Self::Output {
        &mut self.data[row * self.ncols..(row + 1) * self.ncols]
    }
}
/// Original packed lower triangle; row i contains i+1 consecutive values.
#[derive(Clone, Debug)]
pub struct NlmLtriangMatrix {
    pub n: usize,
    pub data: Vec<f64>,
}
impl Index<usize> for NlmLtriangMatrix {
    type Output = [f64];
    fn index(&self, row: usize) -> &Self::Output {
        let start = row * (row + 1) / 2;
        &self.data[start..start + row + 1]
    }
}
impl IndexMut<usize> for NlmLtriangMatrix {
    fn index_mut(&mut self, row: usize) -> &mut Self::Output {
        let start = row * (row + 1) / 2;
        &mut self.data[start..start + row + 1]
    }
}

/// Rust allocation binding: retain C's recoverable allocation-failure path.
pub(crate) fn try_zeroed_f64(length: usize) -> Option<Vec<f64>> {
    let mut data = Vec::new();
    data.try_reserve_exact(length).ok()?;
    data.resize(length, 0.0);
    Some(data)
}
pub fn nlm_dense_matrix_new(nrows: usize, ncols: usize) -> Option<NlmDenseMatrix> {
    assert!(nrows > 0);
    Some(NlmDenseMatrix {
        nrows,
        ncols,
        data: try_zeroed_f64(nrows.checked_mul(ncols)?)?,
    })
}
pub fn nlm_ltriang_matrix_new(n: usize) -> Option<NlmLtriangMatrix> {
    assert!(n > 0);
    Some(NlmLtriangMatrix {
        n,
        data: try_zeroed_f64(n.checked_mul(n.checked_add(1)?)? / 2)?,
    })
}
pub fn nlm_dense_matrix_free<T>(matrix: &mut Option<T>) {
    *matrix = None;
}

pub fn nlm_factor_ltriang_pos_def(a: &mut NlmLtriangMatrix, n: usize) {
    for i in 0..n {
        for j in 0..i {
            let mut temp = a[i][j];
            for k in 0..j {
                temp -= a[i][k] * a[j][k];
            }
            a[i][j] = temp / a[j][j];
        }
        let mut temp = a[i][i];
        for k in 0..i {
            temp -= a[i][k] * a[i][k];
        }
        a[i][i] = temp.sqrt();
    }
}
pub fn nlm_solve_ltriang_pos_def(x: &mut [f64], n: usize, l: &NlmLtriangMatrix) {
    for i in 0..n {
        let mut temp = x[i];
        for j in 0..i {
            temp -= l[i][j] * x[j];
        }
        x[i] = temp / l[i][i];
    }
    for j in (0..n).rev() {
        x[j] /= l[j][j];
        for i in 0..j {
            x[i] -= l[j][i] * x[j];
        }
    }
}
pub fn nlm_euclidean_norm(v: &[f64], n: usize) -> f64 {
    let mut sum = 1.0;
    let mut scale = 0.0;
    for i in 0..n {
        if v[i] != 0.0 {
            let absvi = v[i].abs();
            if scale < absvi {
                sum = 1.0 + sum * (scale / absvi) * (scale / absvi);
                scale = absvi;
            } else {
                sum += (absvi / scale) * (absvi / scale);
            }
        }
    }
    scale * sum.sqrt()
}
pub fn nlm_add_vectors(y: &mut [f64], n: usize, alpha: f64, x: &[f64]) {
    for i in 0..n {
        y[i] += alpha * x[i];
    }
}
pub fn nlm_step_bound(x: &[f64], n: usize, step_x: &[f64], max: f64) -> f64 {
    let mut alpha = max;
    for i in 0..n {
        let alpha_i = -x[i] / step_x[i];
        if alpha_i >= 0.0 && alpha_i < alpha {
            alpha = alpha_i;
        }
    }
    alpha
}
