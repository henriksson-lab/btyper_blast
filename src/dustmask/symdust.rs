// NCBI BLAST+ 2.17.0+, symdust.cpp and symdust.hpp.
// Author: Aleksandr Morgulis. See LICENSE-NCBI.

use crate::util::random_gen::CRandom;
use std::collections::VecDeque;

pub const DEFAULT_LEVEL: u32 = 20;
pub const DEFAULT_WINDOW: u32 = 64;
pub const DEFAULT_LINKER: u32 = 1;
const TRIPLET_MASK: u8 = 0x3f;

/// Flattened C++ CSymDustMasker::CIupac2Ncbi2na_converter.
#[derive(Debug)]
pub struct CIupac2Ncbi2naConverter {
    m_random: CRandom,
}

impl CIupac2Ncbi2naConverter {
    /// Matches C++ conversion operator. Input is already IUPAC-coded; uppercase
    /// N consumes the persistent default RNG, other ambiguities convert to A.
    pub fn convert(&mut self, r: u8) -> u8 {
        match r {
            67 => 1,
            71 => 2,
            84 => 3,
            78 => (self.m_random.get_rand() & 0x3) as u8,
            _ => 0,
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Perfect {
    bounds_: (u32, u32),
    score_: u32,
    len_: u32,
}

impl Perfect {
    /// Matches C++ CSymDustMasker::perfect constructor; bounds are inclusive.
    pub fn new(start: u32, stop: u32, score: u32, len: u32) -> Self {
        Self {
            bounds_: (start, stop),
            score_: score,
            len_: len,
        }
    }
}

/// Flattened C++ CSymDustMasker::triplets. Perfect-list and threshold aliases are
/// passed to methods instead of stored, so Rust can borrow the masker's fields
/// separately without shared mutable references or per-base runtime locks.
#[derive(Debug)]
pub struct Triplets {
    triplet_list_: VecDeque<u8>,
    start_: u32,
    stop_: u32,
    max_size_: u32,
    low_k_: u8,
    l: u32,
    c_w: [u8; 64],
    c_v: [u8; 64],
    r_w: u32,
    r_v: u32,
    num_diff: u32,
}

impl Triplets {
    /// Matches C++ triplets constructor.
    pub fn new(window: u32, low_k: u8) -> Self {
        Self {
            triplet_list_: VecDeque::new(),
            start_: 0,
            stop_: 0,
            max_size_: window - 2,
            low_k_: low_k,
            l: 0,
            c_w: [0; 64],
            c_v: [0; 64],
            r_w: 0,
            r_v: 0,
            num_diff: 0,
        }
    }

    pub fn start(&self) -> u32 {
        self.start_
    }
    pub fn stop(&self) -> u32 {
        self.stop_
    }

    /// Matches C++ triplets::add_triplet_info; implicit unused this is removed.
    pub fn add_triplet_info(r: &mut u32, c: &mut [u8; 64], t: u8) {
        *r = r.wrapping_add(u32::from(c[t as usize]));
        c[t as usize] = c[t as usize].wrapping_add(1);
    }

    /// Matches C++ triplets::rem_triplet_info.
    pub fn rem_triplet_info(r: &mut u32, c: &mut [u8; 64], t: u8) {
        c[t as usize] = c[t as usize].wrapping_sub(1);
        *r = r.wrapping_sub(u32::from(c[t as usize]));
    }

    /// Matches C++ triplets::needs_processing.
    pub fn needs_processing(&self, thresholds: &[u32]) -> bool {
        let count = self.stop_.wrapping_sub(self.l);
        (count as usize) < self.triplet_list_.len()
            && 10u32.wrapping_mul(self.r_w) > thresholds[count as usize]
    }

    /// Matches C++ triplets::shift_high, including its unchanged suffix state.
    pub fn shift_high(&mut self, t: u8, p: &mut VecDeque<Perfect>) -> bool {
        let s = *self.triplet_list_.back().unwrap();
        self.triplet_list_.pop_back();
        Self::rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
        if self.c_w[s as usize] == 0 {
            self.num_diff = self.num_diff.wrapping_sub(1);
        }
        self.start_ += 1;
        self.triplet_list_.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        self.stop_ += 1;
        if self.num_diff <= 1 {
            p.insert(0, Perfect::new(self.start_, self.stop_ + 1, 0, 0));
            false
        } else {
            true
        }
    }

    /// Matches C++ triplets::shift_window.
    pub fn shift_window(&mut self, t: u8, p: &mut VecDeque<Perfect>) -> bool {
        if self.triplet_list_.len() >= self.max_size_ as usize {
            if self.num_diff <= 1 {
                return self.shift_high(t, p);
            }
            let s = *self.triplet_list_.back().unwrap();
            self.triplet_list_.pop_back();
            Self::rem_triplet_info(&mut self.r_w, &mut self.c_w, s);
            if self.c_w[s as usize] == 0 {
                self.num_diff = self.num_diff.wrapping_sub(1);
            }
            if self.l == self.start_ {
                self.l += 1;
                Self::rem_triplet_info(&mut self.r_v, &mut self.c_v, s);
            }
            self.start_ += 1;
        }
        self.triplet_list_.push_front(t);
        if self.c_w[t as usize] == 0 {
            self.num_diff += 1;
        }
        Self::add_triplet_info(&mut self.r_w, &mut self.c_w, t);
        Self::add_triplet_info(&mut self.r_v, &mut self.c_v, t);
        if self.c_v[t as usize] > self.low_k_ {
            let mut off = (self.triplet_list_.len() as u32)
                .wrapping_sub(self.l.wrapping_sub(self.start_))
                .wrapping_sub(1);
            loop {
                Self::rem_triplet_info(
                    &mut self.r_v,
                    &mut self.c_v,
                    self.triplet_list_[off as usize],
                );
                self.l += 1;
                let last = self.triplet_list_[off as usize] == t;
                off = off.wrapping_sub(1);
                if last {
                    break;
                }
            }
        }
        self.stop_ += 1;
        if self.triplet_list_.len() >= self.max_size_ as usize && self.num_diff <= 1 {
            p.clear();
            p.insert(0, Perfect::new(self.start_, self.stop_ + 1, 0, 0));
            false
        } else {
            true
        }
    }

    /// Matches C++ triplets::find_perfect. A deque insertion index preserves the
    /// original list iterator's position, including returning the inserted node.
    pub fn find_perfect(&self, p: &mut VecDeque<Perfect>, thresholds: &[u32]) {
        let mut count = self.stop_.wrapping_sub(self.l);
        let mut counts = self.c_v;
        let mut score = self.r_v;
        let mut perfect_iter = 0;
        let mut max_perfect_score: u32 = 0;
        let mut max_len: u32 = 0;
        let mut pos = self.l.wrapping_sub(1);
        let mut it = count as usize;
        while it != self.triplet_list_.len() {
            let cnt = counts[self.triplet_list_[it] as usize];
            Self::add_triplet_info(&mut score, &mut counts, self.triplet_list_[it]);
            if cnt > 0 && score.wrapping_mul(10) > thresholds[count as usize] {
                while perfect_iter != p.len() && pos <= p[perfect_iter].bounds_.0 {
                    if max_perfect_score == 0
                        || max_len.wrapping_mul(p[perfect_iter].score_)
                            > max_perfect_score.wrapping_mul(p[perfect_iter].len_)
                    {
                        max_perfect_score = p[perfect_iter].score_;
                        max_len = p[perfect_iter].len_;
                    }
                    perfect_iter += 1;
                }
                if max_perfect_score == 0
                    || score.wrapping_mul(max_len) >= max_perfect_score.wrapping_mul(count)
                {
                    max_perfect_score = score;
                    max_len = count;
                    p.insert(
                        perfect_iter,
                        Perfect::new(pos, self.stop_ + 1, max_perfect_score, count),
                    );
                }
            }
            it += 1;
            count += 1;
            pos = pos.wrapping_sub(1);
        }
    }
}

#[derive(Debug)]
pub struct CSymDustMasker {
    level_: u32,
    window_: u32,
    linker_: u32,
    low_k_: u8,
    p: VecDeque<Perfect>,
    thresholds_: Vec<u32>,
    converter_: CIupac2Ncbi2naConverter,
}

impl CSymDustMasker {
    /// Matches C++ CSymDustMasker constructor, including parameter clamping.
    pub fn new(level: u32, window: u32, linker: u32) -> Self {
        let level_ = if (2..=64).contains(&level) {
            level
        } else {
            DEFAULT_LEVEL
        };
        let window_ = if (8..=64).contains(&window) {
            window
        } else {
            DEFAULT_WINDOW
        };
        let linker_ = if (1..=32).contains(&linker) {
            linker
        } else {
            DEFAULT_LINKER
        };
        let mut result = Self {
            level_,
            window_,
            linker_,
            low_k_: (level_ / 5) as u8,
            p: VecDeque::new(),
            thresholds_: Vec::with_capacity((window_ - 2) as usize),
            converter_: CIupac2Ncbi2naConverter {
                m_random: CRandom::new_default(),
            },
        };
        result.thresholds_.push(1);
        for i in 1..window_ - 2 {
            result.thresholds_.push(i * result.level_);
        }
        result
    }

    /// Matches C++ save_masked_regions; result and stored bounds are inclusive.
    pub fn save_masked_regions(&mut self, res: &mut Vec<(u32, u32)>, wstart: u32, start: u32) {
        if !self.p.is_empty() {
            let b = self.p.back().unwrap().bounds_;
            if b.0 < wstart {
                let b1 = (b.0 + start, b.1 + start);
                if !res.is_empty() {
                    let s = res.last().unwrap().1;
                    if s + self.linker_ >= b1.0 {
                        res.last_mut().unwrap().1 = s.max(b1.1);
                    } else {
                        res.push(b1);
                    }
                } else {
                    res.push(b1);
                }
                while !self.p.is_empty() && self.p.back().unwrap().bounds_.0 < wstart {
                    self.p.pop_back();
                }
            }
        }
    }

    /// Matches the C++ ranged operator(). Sequence bytes are IUPAC, already
    /// normalized by the query adapter. An owned Vec replaces unique_ptr<vector>.
    pub fn mask_range(&mut self, seq: &[u8], mut start: u32, mut stop: u32) -> Vec<(u32, u32)> {
        let mut res = Vec::new();
        if seq.is_empty() {
            return res;
        }
        assert!(seq.len() <= i32::MAX as usize);
        if stop >= seq.len() as u32 {
            stop = seq.len() as u32 - 1;
        }
        if start > stop {
            start = stop;
        }
        while stop > 2 + start {
            self.p.clear();
            let mut w = Triplets::new(self.window_, self.low_k_);
            let mut it = start;
            let c1 = seq[it as usize];
            it += 1;
            let c2 = seq[it as usize];
            let mut t = (self.converter_.convert(c1) << 2) + self.converter_.convert(c2);
            it = start + w.stop() + 2;
            let mut done = false;
            while !done && it <= stop {
                self.save_masked_regions(&mut res, w.start(), start);
                t = ((t << 2) & TRIPLET_MASK) + (self.converter_.convert(seq[it as usize]) & 0x3);
                it += 1;
                if w.shift_window(t, &mut self.p) {
                    if w.needs_processing(&self.thresholds_) {
                        w.find_perfect(&mut self.p, &self.thresholds_);
                    }
                } else {
                    while it <= stop {
                        self.save_masked_regions(&mut res, w.start(), start);
                        t = ((t << 2) & TRIPLET_MASK)
                            + (self.converter_.convert(seq[it as usize]) & 0x3);
                        if w.shift_window(t, &mut self.p) {
                            done = true;
                            break;
                        }
                        it += 1;
                    }
                }
            }
            let mut wstart = w.start();
            while !self.p.is_empty() {
                self.save_masked_regions(&mut res, wstart, start);
                wstart += 1;
            }
            if w.start() > 0 {
                start += w.start();
            } else {
                break;
            }
        }
        res
    }

    /// Matches the whole-sequence C++ operator(), preserving default RNG state
    /// across calls on this same masker. The query adapter creates one per query.
    pub fn mask(&mut self, seq: &[u8]) -> Vec<(u32, u32)> {
        self.mask_range(seq, 0, (seq.len() as u32).wrapping_sub(1))
    }
}
