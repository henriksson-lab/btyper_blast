// NCBI BLAST+ 2.17.0+, random_gen.cpp and random_gen.hpp.
// Authors: Clifford Clausen, Denis Vakatov, Jim Ostell, Jonathan Kans,
// Greg Schuler, Eugene Vasilchenko. See LICENSE-NCBI.

const K_STATE_SIZE: usize = 33;
const K_STATE_OFFSET: i32 = 12;

/// LFG modes required by the database writer and symmetric DUST converter.
#[derive(Debug)]
pub struct CRandom {
    m_state: [u32; K_STATE_SIZE],
    m_rj: i32,
    m_rk: i32,
    m_seed: u32,
}

impl CRandom {
    /// Matches C++ `CRandom::CRandom(EGetRandMethod)` for default LFG mode.
    /// System randomness is outside the supported algorithm paths.
    pub fn new_default() -> Self {
        let mut result = Self {
            m_state: [0; K_STATE_SIZE],
            m_rj: 0,
            m_rk: 0,
            m_seed: 0,
        };
        result.reset();
        result
    }

    /// Matches C++ `CRandom::Reset(void)` for LFG mode. This uses the original
    /// startup state directly; it is NOT equivalent to SetSeed(0).
    pub fn reset(&mut self) {
        const SM_STATE: [u32; K_STATE_SIZE] = [
            0xd53f1852, 0xdfc78b83, 0x4f256096, 0xe643df7, 0x82c359bf, 0xc7794dfa, 0xd5e9ffaa,
            0x2c8cb64a, 0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589, 0xa5853646, 0xe71576e2,
            0xdae30df, 0xb09ce711, 0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8, 0xd788d765,
            0x67498087, 0x9d7aba26, 0x261351d4, 0x411ee7ea, 0x393a263, 0x2c5a5835, 0xc115fcd8,
            0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98, 0x6e37bd55,
        ];
        assert!((K_STATE_OFFSET as usize) < K_STATE_SIZE);
        for i in 0..K_STATE_SIZE {
            self.m_state[i] = SM_STATE[i];
        }
        self.m_rj = K_STATE_OFFSET;
        self.m_rk = K_STATE_SIZE as i32 - 1;
    }

    /// Matches C++ `CRandom::CRandom(TValue seed)`.
    pub fn new(seed: u32) -> Self {
        let mut result = Self {
            m_state: [0; K_STATE_SIZE],
            m_rj: 0,
            m_rk: 0,
            m_seed: 0,
        };
        result.set_seed(seed);
        result
    }

    /// Matches C++ `CRandom::SetSeed` for LFG mode.
    pub fn set_seed(&mut self, seed: u32) {
        assert!((K_STATE_OFFSET as usize) < K_STATE_SIZE);
        self.m_state[0] = seed;
        self.m_seed = seed;
        for i in 1..K_STATE_SIZE {
            self.m_state[i] = 1103515245u32
                .wrapping_mul(self.m_state[i - 1])
                .wrapping_add(12345);
        }
        self.m_rj = K_STATE_OFFSET;
        self.m_rk = K_STATE_SIZE as i32 - 1;
        for _ in 0..10 * K_STATE_SIZE {
            self.get_rand();
        }
    }

    /// Matches C++ `CRandom::x_GetRand32Bits` for LFG mode.
    pub fn x_get_rand32_bits(&mut self) -> u32 {
        let r = self.m_state[self.m_rk as usize].wrapping_add(self.m_state[self.m_rj as usize]);
        self.m_rj -= 1;
        self.m_state[self.m_rk as usize] = r;
        self.m_rk -= 1;
        if self.m_rk < 0 {
            self.m_rk = K_STATE_SIZE as i32 - 1;
        } else if self.m_rj < 0 {
            self.m_rj = K_STATE_SIZE as i32 - 1;
        }
        r
    }

    /// Matches C++ `CRandom::GetRand(void)`.
    pub fn get_rand(&mut self) -> u32 {
        self.x_get_rand32_bits() >> 1
    }
}
