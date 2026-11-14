// NCBI BLAST+ 2.17.0+, selected functions from ncbi_math.c. See LICENSE-NCBI.

/// Matches C `BLAST_Expm1`; retain the upstream polynomial and cutoffs.
// Retain polynomial nesting verbatim; rustfmt is extremely slow on this expression.
#[rustfmt::skip]
pub fn blast_expm1(x: f64) -> f64 {
    let absx = if x >= 0.0 { x } else { -x };
    if absx > 0.33 {
        return x.exp() - 1.0;
    }
    if absx < 1.0e-16 {
        return x;
    }
    x * (1.0
        + x * (1.0 / 2.0
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0 + x / 6227020800.0))))))))))))
}

/// Matches C `BLAST_Log1p`; do not replace with the standard library algorithm.
pub fn blast_log1p(x: f64) -> f64 {
    let absx = if x >= 0.0 { x } else { -x };
    if absx >= 0.2 {
        return (x + 1.0).ln();
    }
    let mut i = 0;
    let mut sum = 0.0;
    let mut y = x;
    while i < 500 {
        i += 1;
        sum += y / f64::from(i);
        if (if y >= 0.0 { y } else { -y }) < f64::EPSILON {
            break;
        }
        y *= x;
        i += 1;
        sum -= y / f64::from(i);
        if y < f64::EPSILON {
            break;
        }
        y *= x;
    }
    sum
}

/// Matches C `BLAST_Gcd`; callers must avoid C's undefined ABS(i32::MIN).
pub fn blast_gcd(mut a: i32, mut b: i32) -> i32 {
    b = if b >= 0 { b } else { -b };
    if b > a {
        let c = a;
        a = b;
        b = c;
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// Matches C `BLAST_Gdb3` (the original symbol spells Gdb, not Gcd).
pub fn blast_gdb3(a: &mut i32, b: &mut i32, c: &mut i32) -> i32 {
    let g = if *b == 0 {
        blast_gcd(*a, *c)
    } else {
        blast_gcd(*a, blast_gcd(*b, *c))
    };
    if g > 1 {
        *a /= g;
        *b /= g;
        *c /= g;
    }
    g
}

/// Matches C `BLAST_Nint` with the target ABI's C long width.
/// Inputs must be finite and round to a representable C long.
pub fn blast_nint(mut x: f64) -> std::os::raw::c_long {
    x += if x >= 0.0 { 0.5 } else { -0.5 };
    x as std::os::raw::c_long
}

/// Matches C `BLAST_Powi`; n must not be i32::MIN (undefined negation in C).
pub fn blast_powi(mut x: f64, mut n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }
    if x == 0.0 {
        if n < 0 {
            return f64::INFINITY;
        }
        return 0.0;
    }
    if n < 0 {
        x = 1.0 / x;
        n = -n;
    }
    let mut y = 1.0;
    while n > 0 {
        if n & 1 != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}
