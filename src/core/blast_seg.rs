// NCBI BLAST+ 2.17.0+, selected blast_seg.c/h definitions.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
use super::blast_lookup::SSeqRange;
use std::collections::VecDeque;
#[path = "seg_tables.rs"]
mod tables;
use tables::{LNFACT, LOG_WIN10};

pub const K_SEG_WINDOW: i32 = 12;
pub const K_SEG_LOCUT: f64 = 2.2;
pub const K_SEG_HICUT: f64 = 2.5;
pub const K_PROT_ALPHABET: i32 = 2;
const LN_2: f64 = 0.69314718055994530941723212145818;
const FENCE_SENTRY: u8 = 201;

#[derive(Clone, Debug)]
pub struct SegParameters {
    pub window: i32,
    pub locut: f64,
    pub hicut: f64,
    pub period: i32,
    pub hilenmin: i32,
    pub overlaps: bool,
    pub maxtrim: i32,
    pub maxbogus: i32,
}
#[derive(Debug)]
pub struct Alpha {
    pub alphabet: i32,
    pub alphasize: i32,
    pub lnalphasize: f64,
    pub alphaindex: [i32; 128],
    pub alphaflag: [bool; 128],
}

/// The immutable root buffer and its current offset replace the mutable C
/// sequence pointer. Parent length is the only parent field SEG reads.
#[derive(Debug)]
pub struct SSequence<'a> {
    pub buffer: &'a [u8],
    pub seq_offset: usize,
    pub parent_length: i32,
    pub palpha: Option<&'a Alpha>,
    pub start: i32,
    pub length: i32,
    pub bogus: i32,
    pub punctuation: bool,
    pub composition: Option<Vec<i32>>,
    pub state: Option<Vec<i32>>,
    pub entropy: f64,
}
#[derive(Debug)]
pub struct SSeg {
    pub begin: i32,
    pub end: i32,
    pub next: Option<Box<SSeg>>,
}

// C s_SegFree unlinks iteratively; retain that teardown behavior rather than
// recursively dropping a long owned next chain.
impl Drop for SSeg {
    fn drop(&mut self) {
        let mut next = self.next.take();
        while let Some(mut node) = next {
            next = node.next.take();
        }
    }
}

pub fn s_s_sequence_new<'a>() -> SSequence<'a> {
    SSequence {
        buffer: &[],
        seq_offset: 0,
        parent_length: 0,
        palpha: None,
        start: 0,
        length: 0,
        bogus: 0,
        punctuation: false,
        composition: None,
        state: None,
        entropy: 0.0,
    }
}
pub fn s_has_dash(win: &SSequence<'_>) -> bool {
    let mut seq = win.seq_offset;
    let max = seq + win.length as usize;
    while seq < max {
        let residue = win.buffer[seq];
        seq += 1;
        if residue == b'-' {
            return true;
        }
    }
    false
}
pub fn s_state_cmp(s1: i32, s2: i32) -> i32 {
    s2.wrapping_sub(s1)
}
pub fn s_comp_on(win: &mut SSequence<'_>) {
    let alpha = win.palpha.unwrap();
    let mut comp = vec![0i32; alpha.alphasize as usize];
    let mut seq = win.seq_offset;
    let max = seq + win.length as usize;
    while seq < max {
        let letter = win.buffer[seq] as usize;
        seq += 1;
        if !alpha.alphaflag[letter] {
            let slot = alpha.alphaindex[letter] as usize;
            comp[slot] = comp[slot].wrapping_add(1);
        } else {
            win.bogus = win.bogus.wrapping_add(1);
        }
    }
    win.composition = Some(comp);
}
pub fn s_state_on(win: &mut SSequence<'_>) {
    let alphasize = win.palpha.unwrap().alphasize as usize;
    if win.composition.is_none() {
        s_comp_on(win);
    }
    let composition = win.composition.as_ref().unwrap();
    let mut state = vec![0; alphasize + 1];
    let mut nel = 0;
    for &count in &composition[..alphasize] {
        if count == 0 {
            continue;
        }
        state[nel] = count;
        nel += 1;
    }
    for count in &mut state[nel..] {
        *count = 0;
    }
    state[..nel].sort_by(|a, b| s_state_cmp(*a, *b).cmp(&0));
    win.state = Some(state);
}
pub fn s_open_win<'a>(parent: &SSequence<'a>, start: i32, length: i32) -> Option<SSequence<'a>> {
    if start < 0 || length < 0 || start + length > parent.length {
        return None;
    }
    let mut win = SSequence {
        buffer: parent.buffer,
        seq_offset: parent.seq_offset + start as usize,
        parent_length: parent.length,
        palpha: parent.palpha,
        start,
        length,
        bogus: 0,
        punctuation: false,
        entropy: -2.0,
        state: None,
        composition: None,
    };
    s_state_on(&mut win);
    Some(win)
}
pub fn s_entropy(sv: &[i32]) -> f64 {
    let mut total = 0i32;
    let mut i = 0;
    while sv[i] != 0 {
        total = total.wrapping_add(sv[i]);
        i += 1;
    }
    if total == 0 {
        return 0.0;
    }
    let mut ent = 0.0;
    i = 0;
    if total == 10 {
        while sv[i] != 0 {
            ent += (sv[i] as f64) * LOG_WIN10[sv[i] as usize] / LN_2;
            i += 1;
        }
    } else {
        while sv[i] != 0 {
            ent += (sv[i] as f64) * ((sv[i] as f64) / (total as f64)).ln() / LN_2;
            i += 1;
        }
    }
    (ent / (total as f64)).abs()
}
pub fn s_decrement_sv(sv: &mut [i32], class: i32) {
    let mut i = 0;
    loop {
        let svi = sv[i];
        i += 1;
        if svi == 0 {
            break;
        }
        if svi == class && sv[i] < class {
            sv[i - 1] = svi - 1;
            break;
        }
    }
}
pub fn s_increment_sv(sv: &mut [i32], class: i32) {
    let mut i = 0;
    loop {
        let value = sv[i];
        i += 1;
        if value == class {
            sv[i - 1] = sv[i - 1].wrapping_add(1);
            break;
        }
    }
}
pub fn s_shift_win1(win: &mut SSequence<'_>) -> bool {
    let length = win.length as usize;
    if win.buffer[win.seq_offset + length] == FENCE_SENTRY {
        return false;
    }
    win.start += 1;
    if (win.start as u32).wrapping_add(length as u32) > win.parent_length as u32 {
        win.start -= 1;
        return false;
    }
    let alpha = win.palpha.unwrap();
    let comp = win.composition.as_mut().unwrap();
    let state = win.state.as_mut().unwrap();
    let mut j = win.buffer[win.seq_offset] as usize;
    if !alpha.alphaflag[j] {
        let slot = alpha.alphaindex[j] as usize;
        let class = comp[slot];
        comp[slot] -= 1;
        s_decrement_sv(state, class);
    } else {
        win.bogus -= 1;
    }
    j = win.buffer[win.seq_offset + length] as usize;
    win.seq_offset += 1;
    if !alpha.alphaflag[j] {
        let slot = alpha.alphaindex[j] as usize;
        let class = comp[slot];
        comp[slot] = comp[slot].wrapping_add(1);
        s_increment_sv(state, class);
    } else {
        win.bogus += 1;
    }
    if win.entropy > -2.0 {
        win.entropy = s_entropy(win.state.as_ref().unwrap());
    }
    true
}
pub fn s_entropy_on(win: &mut SSequence<'_>) {
    if win.state.is_none() {
        s_state_on(win);
    }
    win.entropy = s_entropy(win.state.as_ref().unwrap());
}
pub fn s_seq_entropy(seq: &SSequence<'_>, window: i32, maxbogus: i32) -> Option<Vec<f64>> {
    let downset = (window + 1) / 2 - 1;
    let upset = window - downset;
    if window > seq.length {
        return None;
    }
    let mut h = vec![-1.0; seq.length as usize];
    let mut win = s_open_win(seq, 0, window).unwrap();
    s_entropy_on(&mut win);
    let first = downset;
    let last = seq.length - upset;
    for i in first..=last {
        if seq.punctuation && s_has_dash(&win) {
            h[i as usize] = -1.0;
            s_shift_win1(&mut win);
            continue;
        }
        if win.bogus > maxbogus {
            h[i as usize] = -1.0;
            s_shift_win1(&mut win);
            continue;
        }
        h[i as usize] = win.entropy;
        s_shift_win1(&mut win);
    }
    Some(h)
}
pub fn s_find_low(i: i32, limit: i32, hicut: f64, h: &[f64]) -> i32 {
    let mut j = i;
    while j >= limit {
        if h[j as usize] == -1.0 {
            break;
        }
        if h[j as usize] > hicut {
            break;
        }
        j -= 1;
    }
    j + 1
}
pub fn s_find_high(i: i32, limit: i32, hicut: f64, h: &[f64]) -> i32 {
    let mut j = i;
    while j <= limit {
        if h[j as usize] == -1.0 {
            break;
        }
        if h[j as usize] > hicut {
            break;
        }
        j += 1;
    }
    j - 1
}
pub fn s_lnfact(n: i32) -> f64 {
    if (n as usize) < LNFACT.len() {
        LNFACT[n as usize]
    } else {
        (n as f64 + 0.5) * (n as f64).ln() - n as f64 + 0.9189385332
    }
}
pub fn s_ln_perm(sv: &[i32], window_length: i32) -> f64 {
    let mut ans = s_lnfact(window_length);
    let mut i = 0;
    while sv[i] != 0 {
        ans -= s_lnfact(sv[i]);
        i += 1;
    }
    ans
}
pub fn s_ln_ass(sv: &[i32], alphasize: i32) -> f64 {
    let mut ans = LNFACT[alphasize as usize];
    if sv[0] == 0 {
        return ans;
    }
    let mut total = alphasize;
    let mut class = 1;
    let mut svim1 = sv[0];
    let mut i = 0;
    let mut ptr = 0;
    loop {
        i += 1;
        if i == alphasize {
            ans -= s_lnfact(class);
            break;
        }
        ptr += 1;
        let svi = sv[ptr];
        if svi == svim1 {
            class += 1;
            svim1 = svi;
            continue;
        }
        total -= class;
        ans -= s_lnfact(class);
        if svi == 0 {
            ans -= s_lnfact(total);
            break;
        } else {
            class = 1;
            svim1 = svi;
            continue;
        }
    }
    ans
}
pub fn s_get_prob(sv: &[i32], total: i32, alpha: &Alpha) -> f64 {
    let totseq = (total as f64) * alpha.lnalphasize;
    let ans1 = s_ln_ass(sv, alpha.alphasize);
    let mut ans2 = 0.0;
    if ans1 > -100000.0 && sv[0] != i32::MIN {
        ans2 = s_ln_perm(sv, total);
    }
    ans1 + ans2 - totseq
}
pub fn s_trim(
    seq: SSequence<'_>,
    leftend: &mut i32,
    rightend: &mut i32,
    params: &SegParameters,
) -> i16 {
    let mut lend = 0;
    let mut rend = seq.length - 1;
    let mut minlen = 1;
    let maxtrim = params.maxtrim;
    if seq.length - maxtrim > minlen {
        minlen = seq.length - maxtrim;
    }
    let mut minprob = 1.0;
    let mut len = seq.length;
    while len > minlen {
        let mut shift = true;
        let mut i = 0;
        let mut win = s_open_win(&seq, 0, len).unwrap();
        while shift {
            let prob = s_get_prob(win.state.as_ref().unwrap(), len, win.palpha.unwrap());
            if prob < minprob {
                minprob = prob;
                lend = i;
                rend = len + i - 1;
            }
            shift = s_shift_win1(&mut win);
            i += 1;
        }
        len -= 1;
    }
    *leftend += lend;
    *rightend -= seq.length - rend - 1;
    0
}
pub fn s_seg_seq(
    seq: &SSequence<'_>,
    params: &mut SegParameters,
    segs: &mut Option<Box<SSeg>>,
    offset: i32,
) -> i16 {
    if params.window <= 0 {
        return 0;
    }
    if params.locut <= 0.0 {
        params.locut = 0.0;
    }
    if params.hicut <= 0.0 {
        params.hicut = 0.0;
    }
    let window = params.window;
    let locut = params.locut;
    let hicut = params.hicut;
    let downset = (window + 1) / 2 - 1;
    let upset = window - downset;
    let Some(h) = s_seq_entropy(seq, window, params.maxbogus) else {
        return 0;
    };
    let first = downset;
    let last = seq.length - upset;
    let mut lowlim = first;
    let mut i = first;
    while i <= last {
        if h[i as usize] <= locut && h[i as usize] != -1.0 {
            let loi = s_find_low(i, lowlim, hicut, &h);
            let hii = s_find_high(i, last, hicut, &h);
            let mut leftend = loi - downset;
            let mut rightend = hii + upset - 1;
            let temp_seq = s_open_win(seq, leftend, rightend - leftend + 1).unwrap();
            let status = s_trim(temp_seq, &mut leftend, &mut rightend, params);
            if status < 0 {
                return status;
            }
            if i + upset - 1 < leftend {
                let lend = loi - downset;
                let rend = leftend - 1;
                let leftseq = s_open_win(seq, lend, rend - lend + 1).unwrap();
                let mut leftsegs = None;
                let status = s_seg_seq(&leftseq, params, &mut leftsegs, offset + lend);
                if status < 0 {
                    return status;
                }
                // The C source replaces the recursive head's next link here.
                if let Some(mut head) = leftsegs {
                    head.next = segs.take();
                    *segs = Some(head);
                }
            }
            *segs = Some(Box::new(SSeg {
                begin: leftend + offset,
                end: rightend + offset,
                next: segs.take(),
            }));
            i = hii.min(rightend + downset);
            lowlim = i + 1;
        }
        i += 1;
    }
    0
}
pub fn s_merge_segs(seq: &SSequence<'_>, segs: &mut Option<Box<SSeg>>) {
    let hilenmin = 0;
    let Some(mut head) = segs.take() else {
        return;
    };
    if seq.length - 1 - head.end < hilenmin {
        head.end = seq.length - 1;
    }
    let mut seg = &mut head;
    while seg.next.is_some() {
        if seg.begin - seg.next.as_ref().unwrap().end - 1 < hilenmin {
            let mut next = seg.next.take().unwrap();
            if seg.end < next.end {
                seg.end = next.end;
            }
            if seg.begin > next.begin {
                seg.begin = next.begin;
            }
            seg.next = next.next.take();
        } else {
            seg = seg.next.as_mut().unwrap();
        }
    }
    if seg.begin < hilenmin {
        seg.begin = 0;
    }
    *segs = Some(head);
}
pub fn s_segs_to_blast_seq_loc(
    segs: Option<&SSeg>,
    offset: i32,
    seg_locs: &mut VecDeque<SSeqRange>,
) -> i16 {
    let mut current = segs;
    while let Some(seg) = current {
        seg_locs.push_front(SSeqRange {
            left: seg.begin + offset,
            right: seg.end + offset,
        });
        current = seg.next.as_deref();
    }
    0
}
pub fn s_aa20alpha_std() -> Alpha {
    let mut alpha = Alpha {
        alphabet: K_PROT_ALPHABET,
        alphasize: 20,
        lnalphasize: 2.9957322735539909,
        alphaindex: [0; 128],
        alphaflag: [false; 128],
    };
    let mut i = 0;
    for c in 0..128 {
        if c == 1 || (3..=20).contains(&c) || c == 22 {
            alpha.alphaflag[c] = false;
            alpha.alphaindex[c] = i;
            i += 1;
        } else {
            alpha.alphaflag[c] = true;
            alpha.alphaindex[c] = 20;
        }
    }
    alpha
}
pub fn seg_parameters_new_aa() -> SegParameters {
    SegParameters {
        window: K_SEG_WINDOW,
        locut: K_SEG_LOCUT,
        hicut: K_SEG_HICUT,
        period: 1,
        hilenmin: 0,
        overlaps: false,
        maxtrim: 50,
        maxbogus: 2,
    }
}
pub fn s_seg_parameters_check(params: Option<&mut SegParameters>) {
    let Some(params) = params else {
        return;
    };
    if params.window <= 0 {
        params.window = 12;
    }
    if params.locut < 0.0 {
        params.locut = 0.0;
    }
    if params.hicut < 0.0 {
        params.hicut = 0.0;
    }
    if params.locut > params.hicut {
        params.hicut = params.locut;
    }
    if params.maxbogus < 0 {
        params.maxbogus = 0;
    }
    if params.maxbogus > params.window {
        params.maxbogus = params.window;
    }
    if params.period <= 0 {
        params.period = 1;
    }
    if params.maxtrim < 0 {
        params.maxtrim = 0;
    }
}
/// Run original protein SEG over NCBIstdaa data. The input must include a
/// readable trailing byte after `length`, as required by s_ShiftWin1 in C.
/// Locations are inclusive and emitted in the original linked-list order.
pub fn seq_buffer_seg(
    sequence: &[u8],
    length: i32,
    offset: i32,
    params: Option<&mut SegParameters>,
    seg_locs: &mut VecDeque<SSeqRange>,
) -> i16 {
    assert!(
        length >= 0 && sequence.len() > length as usize,
        "SEG needs its original readable trailing byte"
    );
    assert!(
        sequence[..length as usize].iter().all(|r| *r < 128),
        "SEG protein residue exceeds original alphabet table"
    );
    let mut defaults;
    let params = if let Some(params) = params {
        s_seg_parameters_check(Some(params));
        params
    } else {
        defaults = seg_parameters_new_aa();
        s_seg_parameters_check(Some(&mut defaults));
        &mut defaults
    };
    let alpha = s_aa20alpha_std();
    let mut seq = s_s_sequence_new();
    seq.buffer = sequence;
    seq.length = length;
    seq.palpha = Some(&alpha);
    seg_locs.clear();
    let mut segs = None;
    let status = s_seg_seq(&seq, params, &mut segs, 0);
    if status < 0 {
        return status;
    }
    if params.overlaps {
        s_merge_segs(&seq, &mut segs);
    }
    s_segs_to_blast_seq_loc(segs.as_deref(), offset, seg_locs);
    0
}
