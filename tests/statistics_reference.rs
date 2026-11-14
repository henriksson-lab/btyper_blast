use btyper_blast::core::blast_options::*;
use btyper_blast::core::blast_program::EBlastProgramType;
use btyper_blast::core::blast_stat::*;

fn n<'a>(t: &mut impl Iterator<Item = &'a str>) -> i32 {
    t.next().unwrap().parse().unwrap()
}
fn f<'a>(t: &mut impl Iterator<Item = &'a str>) -> f64 {
    f64::from_bits(u64::from_str_radix(t.next().unwrap(), 16).unwrap())
}
fn close(actual: f64, expected: f64, line: usize) {
    if expected.is_nan() {
        assert!(actual.is_nan(), "line {line}: {actual:?} != NaN");
    } else if expected.is_infinite() || expected == 0.0 {
        assert_eq!(actual.to_bits(), expected.to_bits(), "line {line}");
    } else {
        assert!(
            actual.is_finite()
                && actual.is_sign_negative() == expected.is_sign_negative()
                && actual.to_bits().abs_diff(expected.to_bits()) <= 8,
            "line {line}: {actual:.17e} != {expected:.17e}, {} ULP",
            actual.to_bits().abs_diff(expected.to_bits())
        );
    }
}
fn block<'a>(actual: BlastKarlinBlk, t: &mut impl Iterator<Item = &'a str>, line: usize) {
    for value in [
        actual.lambda,
        actual.k,
        actual.log_k,
        actual.h,
        actual.param_c,
    ] {
        close(value, f(t), line);
    }
}
fn verify(text: &str) {
    let mut freq = None;
    let mut to = BlastKarlinBlk {
        lambda: 9.0,
        k: 9.0,
        log_k: 9.0,
        h: 9.0,
        param_c: 9.0,
    };
    for (line, row) in text.lines().enumerate() {
        let mut t = row.split_whitespace();
        match t.next().unwrap() {
            "options" => {
                let program = match n(&mut t) {
                    12 => EBlastProgramType::Blastn,
                    41 => EBlastProgramType::Tblastn,
                    p => panic!("unknown program {p}"),
                };
                let opt = blast_initial_word_options_new(program);
                assert_eq!([opt.window_size, opt.scan_range], [n(&mut t), n(&mut t)]);
                close(opt.gap_trigger, f(&mut t), line);
                close(opt.x_dropoff, f(&mut t), line);
                assert_eq!(opt.program_number as i32, n(&mut t));
            }
            "new" => {
                let low = n(&mut t);
                let high = n(&mut t);
                assert_eq!(blast_score_chk(low, high) as i32, n(&mut t));
                let actual = blast_score_freq_new(low, high);
                assert_eq!(actual.is_some() as i32, n(&mut t));
                if let Some(s) = actual {
                    assert_eq!(
                        [s.score_min, s.score_max, s.obs_min, s.obs_max],
                        [n(&mut t), n(&mut t), n(&mut t), n(&mut t)]
                    );
                    close(s.score_avg, f(&mut t), line);
                    for p in s.sprob {
                        close(p, f(&mut t), line);
                    }
                }
            }
            "zero" => block(blast_karlin_blk_new(), &mut t, line),
            "copy" => {
                let mode = n(&mut t);
                let from = BlastKarlinBlk {
                    lambda: 1.0,
                    k: 2.0,
                    log_k: 3.0,
                    h: 4.0,
                    param_c: 5.0,
                };
                let status = blast_karlin_blk_copy(
                    if mode & 1 != 0 { Some(&mut to) } else { None },
                    if mode & 2 != 0 { Some(&from) } else { None },
                );
                assert_eq!(status as i32, n(&mut t));
                block(to, &mut t, line);
            }
            "frequency" => {
                let lo = n(&mut t);
                let hi = n(&mut t);
                let mut s = blast_score_freq_new(lo, hi).unwrap();
                s.obs_min = lo;
                s.obs_max = hi;
                s.score_avg = f(&mut t);
                for p in &mut s.sprob {
                    *p = f(&mut t);
                }
                freq = Some(s);
            }
            "lambda" => {
                let guess = f(&mut t);
                close(
                    blast_karlin_lambda_nr(freq.as_ref().unwrap(), guess),
                    f(&mut t),
                    line,
                );
            }
            "h" => {
                let lambda = f(&mut t);
                close(
                    blast_karlin_lto_h(freq.as_ref().unwrap(), lambda),
                    f(&mut t),
                    line,
                );
            }
            "k" => {
                let lambda = f(&mut t);
                let h = f(&mut t);
                close(
                    blast_karlin_l_hto_k(freq.as_ref().unwrap(), lambda, h),
                    f(&mut t),
                    line,
                );
            }
            "nr" => {
                let d = n(&mut t);
                let itmax = n(&mut t);
                let max_newton = n(&mut t);
                let expected_itn = n(&mut t);
                let s = freq.as_ref().unwrap();
                let mut itn = -101;
                let value = nlm_karlin_lambda_nr(
                    s, d, s.obs_min, s.obs_max, 0.5, 1e-5, itmax, max_newton, &mut itn,
                );
                assert_eq!(itn, expected_itn, "line {line}");
                close(value, f(&mut t), line);
            }
            "block" => {
                let mut k = BlastKarlinBlk {
                    lambda: 91.0,
                    k: 92.0,
                    log_k: 93.0,
                    h: 94.0,
                    param_c: 95.0,
                };
                let status = blast_karlin_blk_ungapped_calc(Some(&mut k), freq.as_ref());
                assert_eq!(status as i32, n(&mut t));
                block(k, &mut t, line);
            }
            "null" => {
                let mut k = BlastKarlinBlk::default();
                assert_eq!(
                    blast_karlin_blk_ungapped_calc(None, freq.as_ref()) as i32,
                    n(&mut t)
                );
                assert_eq!(
                    blast_karlin_blk_ungapped_calc(Some(&mut k), None) as i32,
                    n(&mut t)
                );
                assert_eq!(blast_karlin_blk_ungapped_calc(None, None) as i32, n(&mut t));
            }
            _ => panic!("unknown record {row}"),
        }
        assert_eq!(t.next(), None, "extra data line {line}");
    }
}
#[test]
fn original_ungapped_statistics() {
    verify(include_str!("fixtures/statistics-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_STATISTICS_FIXTURE"]
fn freshly_compiled_ungapped_statistics() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_STATISTICS_FIXTURE").unwrap()).unwrap());
}
