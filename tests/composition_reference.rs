use btyper_blast::core::blast_stat::*;
fn n<'a>(t: &mut impl Iterator<Item = &'a str>) -> i32 {
    t.next().unwrap().parse().unwrap()
}
fn f<'a>(t: &mut impl Iterator<Item = &'a str>) -> f64 {
    f64::from_bits(u64::from_str_radix(t.next().unwrap(), 16).unwrap())
}
fn bytes(s: &str) -> Vec<u8> {
    if s == "-" {
        return vec![];
    }
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).unwrap())
        .collect()
}
fn close(a: f64, b: f64, line: usize, ulps: u64) {
    if b.is_nan() {
        assert!(a.is_nan(), "line {line}");
    } else if b.is_infinite() || b == 0.0 {
        assert_eq!(a.to_bits(), b.to_bits(), "line {line}");
    } else {
        assert!(
            a.is_finite()
                && a.is_sign_negative() == b.is_sign_negative()
                && a.to_bits().abs_diff(b.to_bits()) <= ulps,
            "line {line}: {a:.17e} != {b:.17e}"
        );
    }
}
fn probs<'a>(p: &BlastResFreq, t: &mut impl Iterator<Item = &'a str>, line: usize) {
    for v in &p.prob {
        close(*v, f(t), line, 0);
    }
}
fn block<'a>(b: BlastKarlinBlk, t: &mut impl Iterator<Item = &'a str>, line: usize) {
    for v in [b.lambda, b.k, b.log_k, b.h, b.param_c] {
        close(v, f(t), line, 8);
    }
}
fn verify(text: &str) {
    let mut sbp: Option<BlastScoreBlk> = None;
    let mut freq: Option<BlastResFreq> = None;
    let mut standard: Option<BlastResFreq> = None;
    let mut comp: Option<BlastResComp> = None;
    let mut sequence = Vec::new();
    let mut score: Option<BlastScoreFreq> = None;
    let mut background = 0;
    let mut nt_background = 0;
    for (line, row) in text.lines().enumerate() {
        let mut t = row.split_whitespace();
        match t.next().unwrap() {
            "null_frequency" => assert_eq!(blast_res_freq_new(None).is_some() as i32, n(&mut t)),
            "ideal_null" => assert_eq!(blast_score_blk_kbp_ideal_calc(None) as i32, n(&mut t)),
            "std_alphabet" => {
                let code = n(&mut t) as u8;
                let size = n(&mut t) as usize;
                let expected = n(&mut t);
                let mut out = [173; 21];
                assert_eq!(
                    blast_get_std_alphabet(code, &mut out[..size]) as i32,
                    expected
                );
                assert_eq!(out.as_slice(), bytes(t.next().unwrap()));
            }
            "background" => {
                assert_eq!(ROBINSON_PROB[background].0 as i32, n(&mut t));
                close(ROBINSON_PROB[background].1, f(&mut t), line, 0);
                background += 1;
            }
            "nt_background" => {
                assert_eq!(NT_PROB[nt_background].0 as i32, n(&mut t));
                close(NT_PROB[nt_background].1, f(&mut t), line, 0);
                nt_background += 1;
            }
            "matrix_bounds" => {
                let matrix = (0..2)
                    .map(|_| (0..2).map(|_| n(&mut t)).collect())
                    .collect();
                let mut sb = BlastScoreBlk {
                    alphabet_code: 0,
                    alphabet_size: 2,
                    alphabet_start: 0,
                    protein_alphabet: false,
                    ambiguous_res: vec![],
                    loscore: 0,
                    hiscore: 0,
                    matrix,
                    kbp_ideal: None,
                    ..Default::default()
                };
                assert_eq!(blast_score_blk_max_score_set(&mut sb), 0);
                assert_eq!([sb.loscore, sb.hiscore], [n(&mut t), n(&mut t)]);
            }
            "scaled_score" | "score_error" => {
                let sb = sbp.as_ref().unwrap();
                // Read the record prefix without changing the current sequence frequency.
                let op = row.split_whitespace().next().unwrap();
                let (lo, missing, p) = if op == "scaled_score" {
                    let scale = f(&mut t);
                    let mut p = standard.as_ref().unwrap().clone();
                    for v in &mut p.prob {
                        *v *= scale;
                    }
                    (sb.loscore, 0, p)
                } else {
                    (sb.loscore + 1, n(&mut t), freq.as_ref().unwrap().clone())
                };
                let mut out = blast_score_freq_new(lo, sb.hiscore).unwrap();
                if op == "score_error" {
                    out.obs_min = -111;
                    out.obs_max = 222;
                    out.score_avg = 333.0;
                    out.sprob.fill(444.0);
                }
                let status = blast_score_freq_calc(
                    if missing & 1 != 0 { None } else { Some(sb) },
                    if missing & 2 != 0 {
                        None
                    } else {
                        Some(&mut out)
                    },
                    &p,
                    if op == "scaled_score" {
                        &p
                    } else {
                        standard.as_ref().unwrap()
                    },
                );
                assert_eq!(status as i32, n(&mut t));
                assert_eq!([out.obs_min, out.obs_max], [n(&mut t), n(&mut t)]);
                close(out.score_avg, f(&mut t), line, 0);
                for v in out.sprob {
                    close(v, f(&mut t), line, 0);
                }
            }
            "alphabet" => {
                let code = n(&mut t) as u8;
                let size = n(&mut t) as usize;
                let count = n(&mut t);
                let ambiguous_res = (0..count).map(|_| n(&mut t) as u8).collect();
                let matrix = (0..size)
                    .map(|_| (0..size).map(|_| n(&mut t)).collect())
                    .collect();
                let sb = BlastScoreBlk {
                    alphabet_code: code,
                    alphabet_size: size as i16,
                    alphabet_start: 0,
                    protein_alphabet: size == 28,
                    ambiguous_res,
                    loscore: 0,
                    hiscore: 0,
                    matrix,
                    kbp_ideal: None,
                    ..Default::default()
                };
                freq = blast_res_freq_new(Some(&sb));
                standard = blast_res_freq_new(Some(&sb));
                comp = Some(blast_res_comp_new(&sb));
                sbp = Some(sb);
            }
            "bounds" => {
                let sb = sbp.as_mut().unwrap();
                assert_eq!(blast_score_blk_max_score_set(sb), 0);
                assert_eq!([sb.loscore, sb.hiscore], [n(&mut t), n(&mut t)]);
            }
            "new_frequency" => {
                let p = freq.as_ref().unwrap();
                assert_eq!(p.alphabet_code as i32, n(&mut t));
                probs(p, &mut t, line);
            }
            "new_comp" => {
                let c = comp.as_ref().unwrap();
                assert_eq!(c.alphabet_code as i32, n(&mut t));
                for v in &c.comp {
                    assert_eq!(*v, n(&mut t));
                }
            }
            "standard" => {
                let sb = sbp.as_ref().unwrap();
                let p = standard.as_mut().unwrap();
                assert_eq!(blast_res_freq_std_comp(sb, p), 0);
                probs(p, &mut t, line);
            }
            "standard_dirty" => {
                let dirty = n(&mut t);
                let p = freq.as_mut().unwrap();
                for (i, v) in p.prob.iter_mut().enumerate() {
                    *v = if dirty != 0 {
                        0.1 * (i + 1) as f64
                    } else {
                        0.0
                    };
                }
                assert_eq!(
                    blast_res_freq_std_comp(sbp.as_ref().unwrap(), p) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            "ideal" => {
                let sb = sbp.as_mut().unwrap();
                assert_eq!(blast_score_blk_kbp_ideal_calc(Some(sb)), 0);
                block(sb.kbp_ideal.unwrap(), &mut t, line);
            }
            "ideal_invalid" => {
                let sb = sbp.as_mut().unwrap();
                let lo = sb.loscore;
                let hi = sb.hiscore;
                sb.loscore = 0;
                sb.hiscore = 0;
                assert_eq!(blast_score_blk_kbp_ideal_calc(Some(sb)) as i32, n(&mut t));
                block(sb.kbp_ideal.unwrap(), &mut t, line);
                sb.loscore = lo;
                sb.hiscore = hi;
            }
            "sequence" => {
                sequence = bytes(t.next().unwrap());
                let c = comp.as_mut().unwrap();
                c.comp.fill(991);
                assert_eq!(
                    blast_res_comp_str(sbp.as_ref(), Some(c), Some(&sequence)) as i32,
                    n(&mut t)
                );
                for v in &c.comp {
                    assert_eq!(*v, n(&mut t), "line {line}");
                }
            }
            "composition_freq" => {
                let p = freq.as_mut().unwrap();
                p.prob.fill(991.0);
                assert_eq!(
                    blast_res_freq_res_comp(sbp.as_ref().unwrap(), Some(p), comp.as_ref()) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            "string_freq" => {
                let p = freq.as_mut().unwrap();
                p.prob.fill(991.0);
                assert_eq!(
                    blast_res_freq_string(sbp.as_ref().unwrap(), p, &sequence) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            "score_freq" => {
                let sb = sbp.as_ref().unwrap();
                let mut s = blast_score_freq_new(sb.loscore, sb.hiscore).unwrap();
                assert_eq!(
                    blast_score_freq_calc(
                        Some(sb),
                        Some(&mut s),
                        freq.as_ref().unwrap(),
                        standard.as_ref().unwrap()
                    ) as i32,
                    n(&mut t)
                );
                assert_eq!(
                    [s.obs_min, s.obs_max],
                    [n(&mut t), n(&mut t)],
                    "line {line}"
                );
                close(s.score_avg, f(&mut t), line, 0);
                for v in &s.sprob {
                    close(*v, f(&mut t), line, 0);
                }
                score = Some(s);
            }
            "karlin" => {
                let mut out = BlastKarlinBlk {
                    lambda: 91.0,
                    k: 92.0,
                    log_k: 93.0,
                    h: 94.0,
                    param_c: 95.0,
                };
                assert_eq!(
                    blast_karlin_blk_ungapped_calc(Some(&mut out), score.as_ref()) as i32,
                    n(&mut t)
                );
                block(out, &mut t, line);
            }
            "normalize" => {
                let mode = n(&mut t);
                let norm = f(&mut t);
                let p = freq.as_mut().unwrap();
                for (i, v) in p.prob.iter_mut().enumerate() {
                    *v = match mode {
                        0 => 0.0,
                        1 => (i + 1) as f64,
                        2 => {
                            if i == 3 {
                                -1.0
                            } else {
                                i as f64
                            }
                        }
                        3 => {
                            if i == 3 {
                                f64::NAN
                            } else {
                                0.1
                            }
                        }
                        _ => 0.1,
                    };
                }
                assert_eq!(
                    blast_res_freq_normalize(sbp.as_ref().unwrap(), p, norm) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            "comp_null" => {
                let missing = n(&mut t);
                let c = comp.as_mut().unwrap();
                c.comp.fill(991);
                assert_eq!(
                    blast_res_comp_str(
                        if missing & 1 != 0 { None } else { sbp.as_ref() },
                        if missing & 2 != 0 { None } else { Some(c) },
                        if missing & 4 != 0 { None } else { Some(&[0]) }
                    ) as i32,
                    n(&mut t)
                );
                for v in &c.comp {
                    assert_eq!(*v, n(&mut t));
                }
            }
            "comp_mismatch" => {
                let c = comp.as_mut().unwrap();
                c.comp.fill(991);
                c.alphabet_code = 0;
                assert_eq!(
                    blast_res_comp_str(sbp.as_ref(), Some(c), Some(&[0])) as i32,
                    n(&mut t)
                );
                for v in &c.comp {
                    assert_eq!(*v, n(&mut t));
                }
            }
            "freq_mismatch" => {
                let p = freq.as_mut().unwrap();
                p.prob.fill(991.0);
                assert_eq!(
                    blast_res_freq_res_comp(sbp.as_ref().unwrap(), Some(p), comp.as_ref()) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            "freq_null" => {
                let missing = n(&mut t);
                let p = freq.as_mut().unwrap();
                assert_eq!(
                    blast_res_freq_res_comp(
                        sbp.as_ref().unwrap(),
                        if missing & 1 != 0 { None } else { Some(p) },
                        if missing & 2 != 0 {
                            None
                        } else {
                            comp.as_ref()
                        }
                    ) as i32,
                    n(&mut t)
                );
                probs(p, &mut t, line);
            }
            _ => panic!("unknown row {line}"),
        }
        assert_eq!(t.next(), None, "extra data line {line}");
    }
}
#[test]
fn original_sequence_composition() {
    verify(include_str!("fixtures/composition-2.17.0.txt"));
}
#[test]
#[ignore = "requires original C output via BLAST_COMPOSITION_FIXTURE"]
fn freshly_compiled_sequence_composition() {
    verify(&std::fs::read_to_string(std::env::var("BLAST_COMPOSITION_FIXTURE").unwrap()).unwrap());
}
