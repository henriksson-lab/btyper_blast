// NCBI BLAST+ 2.17.0+, blast_app_util.hpp/.cpp CBatchSizeMixer.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.
#[derive(Clone, Debug)]
pub struct CBatchSizeMixer {
    k_mix_in: f64,
    m_target_hits: i32,
    m_ratio: f64,
    m_batch_size: i32,
    k_min_batch_size: i32,
    k_max_batch_size: i32,
    k_min_target_hits: i32,
}
impl CBatchSizeMixer {
    pub fn c_batch_size_mixer(max_batch_size: i32) -> Self {
        Self {
            k_mix_in: 0.3,
            m_target_hits: 1_000_000,
            m_ratio: -1.0,
            m_batch_size: 5000,
            k_min_batch_size: 100,
            k_max_batch_size: max_batch_size,
            k_min_target_hits: 1_000_000,
        }
    }
    pub fn set_target_hits(&mut self, target: i32) {
        self.m_target_hits = if target < self.k_min_target_hits {
            self.k_min_target_hits
        } else {
            target
        };
        self.m_batch_size = self.m_target_hits / 200;
    }
    /// Pass -1 for the original default argument (no new extension count).
    pub fn get_batch_size(&mut self, hits: i32) -> i32 {
        if hits >= 0 {
            let ratio = 1.0 * (hits + 1) as f64 / self.m_batch_size as f64;
            self.m_ratio = if self.m_ratio < 0.0 {
                ratio
            } else {
                self.k_mix_in * ratio + (1.0 - self.k_mix_in) * self.m_ratio
            };
            let next_size = 1.0 * self.m_target_hits as f64 / self.m_ratio;
            // C++'s out-of-range Int4 conversion is undefined. The pinned x86
            // oracle yields INT_MIN, then the source lower clamp applies.
            // Define that result portably instead of Rust's saturating cast.
            self.m_batch_size = if next_size >= i32::MIN as f64 && next_size < 2_147_483_648.0 {
                next_size as i32
            } else {
                i32::MIN
            };
        }
        if self.m_batch_size > self.k_max_batch_size {
            self.m_batch_size = self.k_max_batch_size;
            self.m_ratio = -1.0;
        } else if self.m_batch_size < self.k_min_batch_size {
            self.m_batch_size = self.k_min_batch_size;
            self.m_ratio = -1.0;
        }
        self.m_batch_size
    }
    pub fn get_target_hits(&self) -> i32 {
        self.m_target_hits
    }
    pub fn get_ratio(&self) -> f64 {
        self.m_ratio
    }
}
