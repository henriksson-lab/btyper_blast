// NCBI BLAST+ 2.17.0+, blast_stat.c sum-HSP statistics and ncbi_math.c closure.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::ncbi_math::{blast_expm1, blast_log1p, blast_powi};
#[path = "hsp_sumstat_tables.rs"]
mod hsp_sumstat_tables;
use hsp_sumstat_tables::*;
pub const LOGDERIV_ORDER_MAX: i32 = 4;
pub const POLYGAMMA_ORDER_MAX: i32 = LOGDERIV_ORDER_MAX;
pub const NCBIMATH_PI: f64 = 3.1415926535897932384626433832795;
pub const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;
pub const NCBIMATH_LNPI: f64 = 1.1447298858494001741434273513531;
pub const MAX_DIAGS: usize = 20;
pub fn s_log_derivative(order: i32, u: &[f64]) -> f64 {
    if !(0..=LOGDERIV_ORDER_MAX).contains(&order) {
        return f64::INFINITY;
    }
    if order > 0 && u[0] == 0.0 {
        return f64::INFINITY;
    }
    let mut y = [0.0; 5];
    for i in 1..=order as usize {
        y[i] = u[i] / u[0];
    }
    match order {
        0 => {
            if u[0] > 0.0 {
                u[0].ln()
            } else {
                f64::INFINITY
            }
        }
        1 => y[1],
        2 => y[2] - y[1] * y[1],
        3 => y[3] - 3.0 * y[2] * y[1] + 2.0 * y[1] * y[1] * y[1],
        4 => {
            let tmp = y[1] * y[1];
            let mut value = y[4] - 4.0 * y[3] * y[1] - 3.0 * y[2] * y[2] + 12.0 * y[2] * tmp;
            value -= 6.0 * tmp * tmp;
            value
        }
        _ => f64::INFINITY,
    }
}
pub fn s_general_ln_gamma(x: f64, order: i32) -> f64 {
    let xx = x - 1.0;
    let tx = xx + DEFAULT_GAMMA_COEF.len() as f64;
    let mut y = [0.0; 5];
    for i in 0..=order {
        let mut tmp = tx;
        let mut coef = DEFAULT_GAMMA_COEF.len() - 1;
        let mut value;
        if i == 0 {
            value = DEFAULT_GAMMA_COEF[coef] / tmp;
            while coef > 0 {
                coef -= 1;
                tmp -= 1.0;
                value += DEFAULT_GAMMA_COEF[coef] / tmp;
            }
        } else {
            value = DEFAULT_GAMMA_COEF[coef] / blast_powi(tmp, i + 1);
            while coef > 0 {
                coef -= 1;
                tmp -= 1.0;
                value += DEFAULT_GAMMA_COEF[coef] / blast_powi(tmp, i + 1);
            }
            let tmp = blast_factorial(i);
            value *= if i % 2 == 0 { tmp } else { -tmp };
        }
        y[i as usize] = value;
    }
    y[0] += 1.0;
    let mut value = s_log_derivative(order, &y);
    let tmp = tx + 0.5;
    let dim = DEFAULT_GAMMA_COEF.len() as f64;
    match order {
        0 => value += ((NCBIMATH_LNPI + NCBIMATH_LN2) / 2.0) + (xx + 0.5) * tmp.ln() - tmp,
        1 => value += tmp.ln() - dim / tmp,
        2 => value += (tmp + dim) / (tmp * tmp),
        3 => value -= (1.0 + 2.0 * dim / tmp) / (tmp * tmp),
        4 => value += 2.0 * (1.0 + 3.0 * dim / tmp) / (tmp * tmp * tmp),
        _ => {
            let t = blast_factorial(order - 2)
                * blast_powi(tmp, 1 - order)
                * (1.0 + (order - 1) as f64 * dim / tmp);
            if order % 2 == 0 {
                value += t;
            } else {
                value -= t;
            }
        }
    }
    value
}
pub fn s_poly_gamma(mut x: f64, order: i32) -> f64 {
    if !(0..=POLYGAMMA_ORDER_MAX).contains(&order) {
        return f64::INFINITY;
    }
    if x >= 1.0 {
        return s_general_ln_gamma(x, order);
    }
    let mut value;
    if x < 0.0 {
        value = s_general_ln_gamma(1.0 - x, order);
        value = if (order - 1) % 2 == 0 { value } else { -value };
        if order == 0 {
            let mut sx = (NCBIMATH_PI * x).sin();
            sx = sx.abs();
            if (x < -0.1 && (x.ceil() == x || sx < 2.0 * f64::EPSILON)) || sx == 0.0 {
                return f64::INFINITY;
            }
            value += NCBIMATH_LNPI - sx.ln();
        } else {
            let mut y = [0.0; 5];
            x *= NCBIMATH_PI;
            y[0] = x.sin();
            let mut tmp = 1.0;
            for k in 1..=order as usize {
                tmp *= NCBIMATH_PI;
                x += NCBIMATH_PI / 2.0;
                y[k] = tmp * x.sin();
            }
            value -= s_log_derivative(order, &y);
        }
    } else {
        value = s_general_ln_gamma(1.0 + x, order);
        if order == 0 {
            if x == 0.0 {
                return f64::INFINITY;
            }
            value -= x.ln();
        } else {
            let tmp = blast_factorial(order - 1) * blast_powi(x, -order);
            value += if order % 2 == 0 { tmp } else { -tmp };
        }
    }
    value
}
pub fn s_ln_gamma(x: f64) -> f64 {
    s_poly_gamma(x, 0)
}
pub fn blast_factorial(n: i32) -> f64 {
    if n < 0 {
        return 0.0;
    }
    if (n as usize) < K_PRECOMPUTED_FACTORIAL.len() {
        return K_PRECOMPUTED_FACTORIAL[n as usize];
    }
    s_ln_gamma((n + 1) as f64).exp()
}
pub fn blast_ln_gamma_int(n: i32) -> f64 {
    if n > 1 && (n as usize) < K_PRECOMPUTED_FACTORIAL.len() {
        return K_PRECOMPUTED_FACTORIAL[n as usize - 1].ln();
    }
    s_ln_gamma(n as f64)
}
pub fn blast_ln_factorial(x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else {
        s_ln_gamma(x + 1.0)
    }
}
/// Generic callback/borrow replaces C callback plus opaque mutable argument.
pub fn blast_romberg_integrate<A>(
    f: fn(f64, &mut A) -> f64,
    fargs: &mut A,
    p: f64,
    q: f64,
    eps: f64,
    mut epsit: i32,
    mut itmin: i32,
) -> f64 {
    let mut romb = [0.0; MAX_DIAGS];
    let mut epsit_cnt = 0;
    itmin = itmin.max(1);
    itmin = itmin.min(MAX_DIAGS as i32 - 1);
    epsit = epsit.max(1);
    epsit = epsit.min(3);
    let epsck = itmin - epsit;
    let mut npts = 1;
    let mut h = q - p;
    let mut x = f(p, fargs);
    if x.abs() == f64::INFINITY {
        return x;
    }
    let mut y = f(q, fargs);
    if y.abs() == f64::INFINITY {
        return y;
    }
    romb[0] = 0.5 * h * (x + y);
    for i in 1..MAX_DIAGS {
        let mut sum = 0.0;
        x = p + 0.5 * h;
        for _ in 0..npts {
            y = f(x, fargs);
            if y.abs() == f64::INFINITY {
                return y;
            }
            sum += y;
            x += h;
        }
        romb[i] = 0.5 * (romb[i - 1] + h * sum);
        let mut n: std::os::raw::c_long = 4;
        for j in (0..i).rev() {
            romb[j] = (n as f64 * romb[j + 1] - romb[j]) / (n - 1) as f64;
            n = n.wrapping_mul(4);
        }
        if i as i32 > epsck {
            if (romb[1] - romb[0]).abs() > eps * romb[0].abs() {
                epsit_cnt = 0;
            } else {
                epsit_cnt += 1;
                if i as i32 >= itmin && epsit_cnt >= epsit {
                    return romb[0];
                }
            }
        }
        npts *= 2;
        h *= 0.5;
    }
    f64::INFINITY
}
#[derive(Clone, Debug, Default)]
pub struct SRombergCbackArgs {
    pub num_hsps: i32,
    pub num_hsps_minus_2: i32,
    pub adj1: f64,
    pub adj2: f64,
    pub sdvir: f64,
    pub epsilon: f64,
}
pub fn s_outer_integral_cback(x: f64, args: &mut SRombergCbackArgs) -> f64 {
    let y = (x - args.sdvir).exp();
    if y == f64::INFINITY {
        return 0.0;
    }
    if args.num_hsps_minus_2 == 0 {
        return (args.adj2 - y).exp();
    }
    if x == 0.0 {
        return 0.0;
    }
    (args.num_hsps_minus_2 as f64 * x.ln() + args.adj2 - y).exp()
}
pub fn s_inner_integral_cback(s: f64, args: &mut SRombergCbackArgs) -> f64 {
    args.adj2 = args.adj1 - s;
    args.sdvir = s / args.num_hsps as f64;
    let mx = if s > 0.0 { args.sdvir + 3.0 } else { 3.0 };
    blast_romberg_integrate(s_outer_integral_cback, args, 0.0, mx, args.epsilon, 0, 1)
}
pub fn s_blast_sum_p_calc(r: i32, s: f64) -> f64 {
    const K_SUMP_EPSILON: f64 = 0.002;
    if r == 1 {
        if s > 8.0 {
            return (-s).exp();
        }
        return -blast_expm1(-(-s).exp());
    }
    if r < 1 {
        return 0.0;
    }
    if r < 8 {
        if s <= -2.3 * r as f64 {
            return 1.0;
        }
    } else if r < 15 {
        if s <= -2.5 * r as f64 {
            return 1.0;
        }
    } else if r < 27 {
        if s <= -3.0 * r as f64 {
            return 1.0;
        }
    } else if r < 51 {
        if s <= -3.4 * r as f64 {
            return 1.0;
        }
    } else if r < 101 && s <= -4.0 * r as f64 {
        return 1.0;
    }
    let stddev = (r as f64).sqrt();
    let stddev4 = 4.0 * stddev;
    let r1 = r - 1;
    if r > 100 {
        let est_mean = (-r * r1) as f64;
        if s <= est_mean - stddev4 {
            return 1.0;
        }
    }
    let logr = (r as f64).ln();
    let mean = r as f64 * (1.0 - logr) - 0.5;
    if s <= mean - stddev4 {
        return 1.0;
    }
    let (t, mut itmin) = if s >= mean {
        (s + 6.0 * stddev, 1)
    } else {
        (mean + 6.0 * stddev, 2)
    };
    let mut args = SRombergCbackArgs {
        num_hsps: r,
        num_hsps_minus_2: r - 2,
        adj1: (r - 2) as f64 * logr - blast_ln_gamma_int(r1) - blast_ln_gamma_int(r),
        epsilon: K_SUMP_EPSILON,
        ..Default::default()
    };
    loop {
        let d = blast_romberg_integrate(
            s_inner_integral_cback,
            &mut args,
            s,
            t,
            K_SUMP_EPSILON,
            0,
            itmin,
        );
        if d == f64::INFINITY {
            return d;
        }
        let retry = if s < mean && d < 0.4 {
            let old = itmin;
            itmin += 1;
            old < 4
        } else {
            false
        };
        if !retry {
            return if d < 1.0 { d } else { 1.0 };
        }
    }
}
pub fn s_blast_sum_p(r: i32, s: f64) -> f64 {
    let table: [&[f64]; 3] = [&K_TAB2, &K_TAB3, &K_TAB4];
    let tabsize = [K_TAB2.len() - 1, K_TAB3.len() - 1, K_TAB4.len() - 1];
    if r == 1 {
        return -blast_expm1(-(-s).exp());
    }
    if r <= 4 {
        if r < 1 {
            return 0.0;
        }
        let r1 = r - 1;
        if s >= (r * r + r1) as f64 {
            let a = blast_ln_gamma_int(r + 1);
            return r as f64 * (r1 as f64 * s.ln() - s - a - a).exp();
        }
        if s > (-2 * r) as f64 {
            let mut a = s + s + (4 * r) as f64;
            let mut i = a as i32;
            a -= i as f64;
            let r2 = (r - 2) as usize;
            i = tabsize[r2] as i32 - i;
            return a * table[r2][i as usize - 1] + (1.0 - a) * table[r2][i as usize];
        }
        return 1.0;
    }
    s_blast_sum_p_calc(r, s)
}
pub fn blast_karlin_pto_e(p: f64) -> f64 {
    if p < 0.0 || p > 1.0 {
        return i32::MIN as f64;
    }
    if p == 1.0 {
        return i32::MAX as f64;
    }
    -blast_log1p(-p)
}
pub fn blast_uneven_gap_sum_e(
    query_start_points: i32,
    subject_start_points: i32,
    num: i16,
    mut xsum: f64,
    query_length: i32,
    subject_length: i32,
    searchsp_eff: i64,
    weight_divisor: f64,
) -> f64 {
    let mut sum_e;
    if num == 1 {
        sum_e = searchsp_eff as f64 * (-xsum).exp();
    } else {
        let pair_search_space = subject_length as f64 * query_length as f64;
        xsum -= pair_search_space.ln()
            + (num - 1) as f64
                * ((query_start_points as f64).ln() + (subject_start_points as f64).ln());
        xsum -= blast_ln_factorial(num as f64);
        let sum_p = s_blast_sum_p(num as i32, xsum);
        sum_e = blast_karlin_pto_e(sum_p) * (searchsp_eff as f64 / pair_search_space);
    }
    if weight_divisor == 0.0 {
        sum_e = i32::MAX as f64;
    } else {
        sum_e /= weight_divisor;
        if sum_e > i32::MAX as f64 {
            sum_e = i32::MAX as f64;
        }
    }
    sum_e
}
