// NCBI BLAST+ 2.17.0+, selected blast_setup_cxx.cpp definition.
// Author: Christiam Camacho. See LICENSE-NCBI.
use super::blast_encoding::{EBlastEncoding, K_NUCL_SENTINEL, K_PROT_SENTINEL, NCBI4NA_TO_BLASTNA};
use super::blast_query_info::BlastQueryInfo;

/// ASN.1 Na-strand discriminants, without generated object-manager classes.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum ENaStrand {
    Unknown = 0,
    Plus = 1,
    Minus = 2,
    Both = 3,
    BothRev = 4,
    Other = 255,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(i32)]
pub enum ESentinelType {
    Sentinels = 0,
    NoSentinels = 1,
}

/// Matches C++ `CalculateSeqBufferLength` for valid encoding/strand pairs.
pub fn calculate_seq_buffer_length(
    sequence_length: u32,
    encoding: EBlastEncoding,
    strand: ENaStrand,
    sentinel: ESentinelType,
) -> u32 {
    if sequence_length == 0 {
        return 0;
    }
    match encoding {
        EBlastEncoding::Ncbi2na => {
            assert_eq!(sentinel, ESentinelType::NoSentinels);
            assert_eq!(strand, ENaStrand::Plus);
            sequence_length / 4 + 1
        }
        EBlastEncoding::Ncbi4na | EBlastEncoding::Nucleotide => {
            if sentinel == ESentinelType::Sentinels {
                if strand == ENaStrand::Both {
                    sequence_length.wrapping_mul(2).wrapping_add(3)
                } else {
                    sequence_length.wrapping_add(2)
                }
            } else if strand == ENaStrand::Both {
                sequence_length.wrapping_mul(2)
            } else {
                sequence_length
            }
        }
        EBlastEncoding::Protein => {
            assert_eq!(sentinel, ESentinelType::Sentinels);
            assert_eq!(strand, ENaStrand::Unknown);
            sequence_length.wrapping_add(2)
        }
        EBlastEncoding::Error => panic!("unsupported encoding"),
    }
}

/// Matches C++ `GetSentinelByte`. Unsupported encoding throws in C++;
/// this low-level Rust binding panics for that unsupported input.
pub fn get_sentinel_byte(encoding: EBlastEncoding) -> u8 {
    match encoding {
        EBlastEncoding::Protein => K_PROT_SENTINEL,
        EBlastEncoding::Ncbi4na | EBlastEncoding::Nucleotide => K_NUCL_SENTINEL,
        EBlastEncoding::Ncbi2na | EBlastEncoding::Error => panic!("unsupported sentinel encoding"),
    }
}

/// Matches C++ `GetSequenceSingleNucleotideStrand`'s defined nucleotide paths.
/// Expanded NCBI4na strand slices replace the sequence-vector input interface;
/// the caller supplies both strands, each of the same nonzero size.
pub fn get_sequence_single_nucleotide_strand(
    plus: &[u8],
    minus: &[u8],
    encoding: EBlastEncoding,
    strand: ENaStrand,
    sentinel: ESentinelType,
) -> Vec<u8> {
    assert!(matches!(strand, ENaStrand::Plus | ENaStrand::Minus));
    assert!(matches!(
        encoding,
        EBlastEncoding::Ncbi4na | EBlastEncoding::Nucleotide
    ));
    assert_eq!(plus.len(), minus.len());
    let size = plus.len();
    let buflen = calculate_seq_buffer_length(size as u32, encoding, strand, sentinel);
    assert_ne!(buflen, 0);
    let mut buffer = vec![0; buflen as usize];
    let mut pos = 0;
    if sentinel == ESentinelType::Sentinels {
        buffer[pos] = get_sentinel_byte(encoding);
        pos += 1;
    }
    let data = if strand == ENaStrand::Plus {
        plus
    } else {
        minus
    };
    buffer[pos..pos + size].copy_from_slice(data);
    if encoding == EBlastEncoding::Nucleotide {
        for i in 0..size {
            assert!(plus[i] < 16);
            buffer[pos + i] = NCBI4NA_TO_BLASTNA[buffer[pos + i] as usize];
        }
    }
    pos += size;
    if sentinel == ESentinelType::Sentinels {
        buffer[pos] = get_sentinel_byte(encoding);
    }
    buffer
}

/// Matches C++ `GetSequenceNucleotideBothStrands`, including its three
/// sentinels. Sequence-vector acquisition becomes explicit input slices.
pub fn get_sequence_nucleotide_both_strands(
    plus: &[u8],
    minus: &[u8],
    encoding: EBlastEncoding,
    sentinel: ESentinelType,
) -> Vec<u8> {
    let plus_data = get_sequence_single_nucleotide_strand(
        plus,
        minus,
        encoding,
        ENaStrand::Plus,
        ESentinelType::NoSentinels,
    );
    let minus_data = get_sequence_single_nucleotide_strand(
        plus,
        minus,
        encoding,
        ENaStrand::Minus,
        ESentinelType::NoSentinels,
    );
    let buflen =
        calculate_seq_buffer_length(plus.len() as u32, encoding, ENaStrand::Both, sentinel);
    let mut retval = vec![0; buflen as usize];
    let mut pos = 0;
    if sentinel == ESentinelType::Sentinels {
        retval[pos] = get_sentinel_byte(encoding);
        pos += 1;
    }
    retval[pos..pos + plus_data.len()].copy_from_slice(&plus_data);
    pos += plus_data.len();
    if sentinel == ESentinelType::Sentinels {
        retval[pos] = get_sentinel_byte(encoding);
        pos += 1;
    }
    retval[pos..pos + minus_data.len()].copy_from_slice(&minus_data);
    pos += minus_data.len();
    if sentinel == ESentinelType::Sentinels {
        retval[pos] = get_sentinel_byte(encoding);
    }
    retval
}

/// Matches C++ `s_QueryInfo_SetContext`. Zero lengths invalidate a context;
/// subsequent nonzero assignments deliberately do not restore is_valid.
pub fn s_query_info_set_context(qinfo: &mut BlastQueryInfo, index: u32, length: u32) {
    assert!(index <= qinfo.last_context as u32);
    let index = index as usize;
    if index != 0 {
        let prev_loc = qinfo.contexts[index - 1].query_offset as u32;
        let prev_len = qinfo.contexts[index - 1].query_length as u32;
        let shift = if prev_len != 0 {
            prev_len.wrapping_add(1)
        } else {
            0
        };
        qinfo.contexts[index].query_offset = prev_loc.wrapping_add(shift) as i32;
        qinfo.contexts[index].query_length = length as i32;
        if length == 0 {
            qinfo.contexts[index].is_valid = false;
        }
    } else {
        qinfo.contexts[0].query_offset = 0;
        qinfo.contexts[0].query_length = length as i32;
        if length == 0 {
            qinfo.contexts[0].is_valid = false;
        }
    }
}
