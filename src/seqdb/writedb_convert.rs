// NCBI BLAST+ 2.17.0+, selected writedb_convert.cpp and writedb_convert.hpp.
// Author: Kevin Bealer. See LICENSE-NCBI.

use super::writedb_general::s_divide_round_up;
use crate::util::random_gen::CRandom;
use std::sync::OnceLock;

#[derive(Debug)]
pub struct CAmbiguousRegion {
    m_start: i32,
    m_end: i32,
    m_value: i32,
}

impl CAmbiguousRegion {
    /// Matches C++ `CAmbiguousRegion::CAmbiguousRegion(int value, int offset)`.
    pub fn new(value: i32, offset: i32) -> Self {
        Self {
            m_start: offset,
            m_end: offset + 1,
            m_value: value,
        }
    }

    /// Matches C++ `CAmbiguousRegion::Append`.
    pub fn append(&mut self, value: i32, offset: i32) -> bool {
        assert!(self.m_end != 0 && self.m_value != 0);
        if value == self.m_value && offset == self.m_end {
            if self.length() < 0xfff {
                self.m_end += 1;
                return true;
            }
        }
        false
    }

    /// Matches C++ `CAmbiguousRegion::Value`.
    pub fn value(&self) -> i32 {
        self.m_value
    }

    /// Matches C++ `CAmbiguousRegion::Length`.
    pub fn length(&self) -> i32 {
        self.m_end - self.m_start
    }

    /// Matches C++ `CAmbiguousRegion::Offset`.
    pub fn offset(&self) -> i32 {
        self.m_start
    }
}

#[derive(Debug)]
pub struct CAmbigDataBuilder {
    m_log2: [i32; 16],
    m_size: i32,
    m_regions: Vec<CAmbiguousRegion>,
    m_random: CRandom,
}

impl CAmbigDataBuilder {
    /// Matches C++ `CAmbigDataBuilder::CAmbigDataBuilder(int sz)`.
    pub fn new(sz: i32) -> Self {
        let mut result = Self {
            m_log2: [0; 16],
            m_size: sz,
            m_regions: Vec::new(),
            m_random: CRandom::new(sz as u32),
        };
        for i in 0..16 {
            result.m_log2[i] = -1;
        }
        result.m_log2[1] = 0;
        result.m_log2[2] = 1;
        result.m_log2[4] = 2;
        result.m_log2[8] = 3;
        result
    }

    /// Matches C++ `CAmbigDataBuilder::Check`.
    pub fn check(&mut self, data: i32, offset: i32) -> i32 {
        if offset >= self.m_size {
            return 0;
        }
        assert_ne!(data, 0);
        let rv = self.m_log2[(data & 0xf) as usize];
        if rv != -1 {
            return rv;
        }
        self.x_add_ambig(data, offset);
        self.x_random(data)
    }

    /// Matches C++ `CAmbigDataBuilder::GetAmbig` (appends to the output).
    pub fn get_ambig(&self, amb: &mut Vec<u8>) {
        let mut new_format = false;
        if self.m_size > 0xffffff {
            new_format = true;
        } else {
            for i in 0..self.m_regions.len() {
                if self.m_regions[i].length() > 0xf {
                    new_format = true;
                    break;
                }
            }
        }
        let num_amb = self.m_regions.len() as i32;
        let amb_words = if new_format {
            0x80000000u32 | (num_amb * 2) as u32
        } else {
            num_amb as u32
        };
        amb.reserve((1 + self.m_regions.len()) * 8);
        s_append_int4(amb, amb_words as i32);
        for i in 0..num_amb as usize {
            if new_format {
                self.x_pack_new_ambig(amb, &self.m_regions[i]);
            } else {
                self.x_pack_old_ambig(amb, &self.m_regions[i]);
            }
        }
    }

    /// Matches C++ `CAmbigDataBuilder::x_PackNewAmbig`.
    pub fn x_pack_new_ambig(&self, amb: &mut Vec<u8>, r: &CAmbiguousRegion) {
        let length_m1 = r.length() - 1;
        assert!(r.value() <= 15);
        assert_eq!(length_m1 >> 12, 0);
        let ch0 = (r.value() << 4) | (length_m1 >> 8);
        let ch1 = length_m1 & 0xff;
        let a1 = [ch0 as u8, ch1 as u8, 0, 0];
        amb.extend_from_slice(&a1);
        s_append_int4(amb, r.offset());
    }

    /// Matches C++ `CAmbigDataBuilder::x_PackOldAmbig`.
    pub fn x_pack_old_ambig(&self, amb: &mut Vec<u8>, r: &CAmbiguousRegion) {
        let length_m1 = r.length() - 1;
        let off = r.offset();
        assert!(r.value() <= 15);
        assert_eq!(length_m1 >> 4, 0);
        assert!(off <= 0xffffff);
        let ch0 = (r.value() << 4) | length_m1;
        let a1 = [
            ch0 as u8,
            ((off >> 16) & 0xff) as u8,
            ((off >> 8) & 0xff) as u8,
            (off & 0xff) as u8,
        ];
        amb.extend_from_slice(&a1);
    }

    /// Matches C++ `CAmbigDataBuilder::x_AddAmbig`.
    pub fn x_add_ambig(&mut self, value: i32, offset: i32) {
        if !self.m_regions.is_empty() {
            if self.m_regions.last_mut().unwrap().append(value, offset) {
                return;
            }
        }
        let r = CAmbiguousRegion::new(value, offset);
        self.m_regions.push(r);
    }

    /// Matches C++ `CAmbigDataBuilder::x_Random`.
    pub fn x_random(&mut self, value: i32) -> i32 {
        if value == 0xf {
            return (self.m_random.get_rand() & 0x3) as i32;
        }
        if value == 0 {
            eprintln!("Error: '0' ambiguity code found, changing to 15.");
            return (self.m_random.get_rand() & 0x3) as i32;
        }
        let bitcount = (value & 1) + ((value >> 1) & 1) + ((value >> 2) & 1) + ((value >> 3) & 1);
        assert!(bitcount >= 2);
        assert!(bitcount <= 3);
        let mut pick = self.m_random.get_rand() % bitcount as u32;
        for i in 0..4 {
            if value & (1 << i) == 0 {
                continue;
            }
            if pick == 0 {
                return i;
            }
            pick -= 1;
        }
        unreachable!("NCBI ambiguity selection exhausted its bits");
    }
}

/// Matches C++ `s_AppendInt4`.
pub fn s_append_int4(outp: &mut Vec<u8>, x: i32) {
    let buf = [
        ((x >> 24) & 0xff) as u8,
        ((x >> 16) & 0xff) as u8,
        ((x >> 8) & 0xff) as u8,
        (x & 0xff) as u8,
    ];
    outp.extend_from_slice(&buf);
}

/// Matches C++ `s_BuildNa4ToNa2Table`.
pub fn s_build_na4_to_na2_table() -> Vec<u8> {
    let mut ctable = vec![0xff; 16];
    for i in 0..4 {
        ctable[1 << i] = i as u8;
    }
    ctable
}

/// Matches the raw-buffer overload of C++ `WriteDB_Ncbi4naToBinary`.
/// Input contains two NCBI4na nibbles per byte. Output buffers must be empty or
/// cleared by the caller: the upstream function ORs into resized sequence bytes
/// and appends ambiguity data. Length and nibble validation belong to the adapter.
pub fn write_db_ncbi4na_to_binary(
    ncbi4na: &[u8],
    byte_length: i32,
    base_length: i32,
    seq: &mut Vec<u8>,
    amb: &mut Vec<u8>,
) {
    static CTABLE: OnceLock<Vec<u8>> = OnceLock::new();
    let ctable = CTABLE.get_or_init(s_build_na4_to_na2_table);
    let inp_bytes = s_divide_round_up(base_length, 2);
    let blast_bytes = base_length / 4 + 1;
    let remainder = base_length & 3;
    let last_byte = blast_bytes - 1;
    let mut ambiguities = CAmbigDataBuilder::new(base_length);
    if inp_bytes != byte_length {
        println!("ib={inp_bytes},n4sz={byte_length}");
    }
    assert_eq!(inp_bytes, byte_length);
    assert!(base_length >= 0 && ncbi4na.len() >= byte_length as usize);
    seq.resize(blast_bytes as usize, 0);
    for i in 0..inp_bytes {
        let inp = ncbi4na[i as usize];
        let mut b1 = inp >> 4;
        let mut b2 = inp & 0xf;
        let c1 = ctable[b1 as usize];
        let c2 = ctable[b2 as usize];
        let mut half: u8 = 0;
        if ((c1 | c2) & 0x80) == 0 {
            half = (c1 << 2) | c2;
        } else {
            if b1 == 0 {
                b1 = 0xf;
            }
            if b2 == 0 && i * 2 + 1 < base_length {
                b2 = 0xf;
            }
            half |= (ambiguities.check(i32::from(b1), i * 2) << 2) as u8;
            half |= ambiguities.check(i32::from(b2), i * 2 + 1) as u8;
        }
        seq[i as usize / 2] |= if i & 1 != 0 { half } else { half << 4 };
    }
    seq[last_byte as usize] &= 255 - 3;
    seq[last_byte as usize] |= remainder as u8;
    ambiguities.get_ambig(amb);
}
