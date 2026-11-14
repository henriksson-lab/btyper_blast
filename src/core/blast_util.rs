// NCBI BLAST+ 2.17.0+, selected functions from blast_util.c.
// Authors: NCBI BLAST developers. See LICENSE-NCBI.

use super::blast_encoding::{EBlastEncoding, AMINOACID_TO_NCBISTDAA, NCBI4NA_TO_BLASTNA};
use super::blast_program::*;

pub const FENCE_SENTRY: u8 = 201;
pub const CODON_LENGTH: i32 = 3;
pub const COMPRESSION_RATIO: i32 = 4;
pub const NUM_FRAMES: u32 = 6;
pub const NUM_STRANDS: u32 = 2;
const NCBI2NA_MASK: u8 = 0x03;

/// Matches C `BlastCompressBlastnaSequence`. Returns its whole allocation;
/// the rolling sequence starts at index 3, with three prefix bytes and three
/// suffix bytes. This differs from the four-residues-per-byte database format.
pub fn blast_compress_blastna_sequence(sequence: &[u8]) -> Vec<u8> {
    let len = sequence.len();
    let mut new_seq = vec![0; len + 3];
    let max_start = 3.min(len);
    let mut curr_letter = 0_i32;
    let mut i = 0;
    while i < max_start {
        curr_letter = curr_letter.wrapping_shl(2) | i32::from(sequence[i] & 3);
        new_seq[3 + i - max_start] = curr_letter as u8;
        i += 1;
    }
    while i < len {
        curr_letter = curr_letter.wrapping_shl(2) | i32::from(sequence[i] & 3);
        new_seq[3 + i - max_start] = curr_letter as u8;
        i += 1;
    }
    for i in 0..max_start {
        curr_letter = curr_letter.wrapping_shl(2);
        new_seq[3 + len - (max_start - i)] = curr_letter as u8;
    }
    new_seq
}

/// Matches C `s_CodonToAA`. Codon values are NCBI4na; codes are NCBIstdaa.
pub fn s_codon_to_aa(codon: &[u8; 3], codes: &[u8; 64]) -> u8 {
    let k_x_residue = AMINOACID_TO_NCBISTDAA[b'X' as usize];
    let mut aa = 0;
    let mapping: [u8; 4] = [8, 2, 1, 4];

    if (codon[0] | codon[1] | codon[2]) > 15 {
        if codon[0] == FENCE_SENTRY || codon[1] == FENCE_SENTRY || codon[2] == FENCE_SENTRY {
            return FENCE_SENTRY;
        }
    }
    for i in 0..4 {
        if codon[0] & mapping[i] != 0 {
            let index0 = i * 16;
            for j in 0..4 {
                if codon[1] & mapping[j] != 0 {
                    let index1 = index0 + j * 4;
                    for k in 0..4 {
                        if codon[2] & mapping[k] != 0 {
                            let index2 = index1 + k;
                            let taa = codes[index2];
                            if aa == 0 {
                                aa = taa;
                            } else if taa != aa {
                                aa = k_x_residue;
                                break;
                            }
                        }
                        if aa == k_x_residue {
                            break;
                        }
                    }
                }
                if aa == k_x_residue {
                    break;
                }
            }
        }
        if aa == k_x_residue {
            break;
        }
    }
    aa
}

/// Matches C `BLAST_GetTranslation`.
///
/// Positive frames use `query_seq` without a leading sentinel. Negative frames
/// use `query_seq_rev` WITH its leading sentinel. `prot_seq` needs room for both
/// protein sentinels. The return value excludes both sentinels.
pub fn blast_get_translation(
    query_seq: &[u8],
    query_seq_rev: &[u8],
    nt_length: i32,
    frame: i16,
    prot_seq: &mut [u8],
    genetic_code: &[u8; 64],
) -> i32 {
    assert!((-3..=3).contains(&frame) && frame != 0);
    assert!(nt_length >= 0);
    let nucl_seq = if frame >= 0 {
        query_seq
    } else {
        &query_seq_rev[1..]
    };
    prot_seq[0] = 0;
    let mut index_prot = 1;
    let mut index = i32::from(frame.abs()) - 1;
    while index < nt_length - 2 {
        let codon = [
            nucl_seq[index as usize],
            nucl_seq[index as usize + 1],
            nucl_seq[index as usize + 2],
        ];
        let residue = s_codon_to_aa(&codon, genetic_code);
        if residue <= 250 || residue == FENCE_SENTRY {
            prot_seq[index_prot as usize] = residue;
            index_prot += 1;
        }
        index += CODON_LENGTH;
    }
    prot_seq[index_prot as usize] = 0;
    index_prot - 1
}

/// Matches C `GetReverseNuclSequence`; the output includes two zero sentinels.
/// Rust returns the owned output buffer in place of the C output pointer.
pub fn get_reverse_nucl_sequence(sequence: &[u8], length: i32) -> Vec<u8> {
    assert!(length >= 0 && length as usize <= sequence.len());
    let conversion_table: [u8; 16] = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];
    let mut rev_sequence = vec![0; length as usize + 2];
    for index in 0..length as usize {
        rev_sequence[length as usize - index] = if sequence[index] == FENCE_SENTRY {
            FENCE_SENTRY
        } else {
            conversion_table[sequence[index] as usize]
        };
    }
    rev_sequence
}

/// Matches C `BLAST_ContextToFrame`, retaining error value `i8::MAX`.
pub fn blast_context_to_frame(prog_number: EBlastProgramType, mut context_number: u32) -> i8 {
    let mut frame = i8::MAX;
    if prog_number == EBlastProgramType::Blastn || prog_number == EBlastProgramType::Mapping {
        if context_number % NUM_STRANDS == 0 {
            frame = 1;
        } else {
            frame = -1;
        }
    } else if blast_query_is_protein(prog_number) || prog_number == EBlastProgramType::PhiBlastn {
        frame = 0;
    } else if prog_number == EBlastProgramType::Blastx
        || prog_number == EBlastProgramType::Tblastx
        || prog_number == EBlastProgramType::RpsTblastn
    {
        context_number %= NUM_FRAMES;
        frame = match context_number {
            0 => 1,
            1 => 2,
            2 => 3,
            3 => -1,
            4 => -2,
            5 => -3,
            _ => unreachable!(),
        };
    }
    frame
}

/// Matches C `BLAST_PackDNA`; includes the final length-remainder byte.
/// Input encoding and ambiguity mapping are intentionally identical to C.
pub fn blast_pack_dna(buffer: &[u8], length: i32, encoding: EBlastEncoding) -> Vec<u8> {
    assert!(length >= 0 && length as usize <= buffer.len());
    let new_length = length / COMPRESSION_RATIO + 1;
    let mut new_buffer = vec![0; new_length as usize];
    let mut index: i32 = 0;
    let mut new_index: i32 = 0;
    while new_index < new_length - 1 {
        let i = index as usize;
        if encoding == EBlastEncoding::Nucleotide {
            new_buffer[new_index as usize] = ((buffer[i] & NCBI2NA_MASK) << 6)
                | ((buffer[i + 1] & NCBI2NA_MASK) << 4)
                | ((buffer[i + 2] & NCBI2NA_MASK) << 2)
                | (buffer[i + 3] & NCBI2NA_MASK);
        } else {
            new_buffer[new_index as usize] =
                ((NCBI4NA_TO_BLASTNA[buffer[i] as usize] & NCBI2NA_MASK) << 6)
                    | ((NCBI4NA_TO_BLASTNA[buffer[i + 1] as usize] & NCBI2NA_MASK) << 4)
                    | ((NCBI4NA_TO_BLASTNA[buffer[i + 2] as usize] & NCBI2NA_MASK) << 2)
                    | (NCBI4NA_TO_BLASTNA[buffer[i + 3] as usize] & NCBI2NA_MASK);
        }
        new_index += 1;
        index += COMPRESSION_RATIO;
    }
    new_buffer[new_index as usize] = (length % COMPRESSION_RATIO) as u8;
    while index < length {
        let shift = match index % COMPRESSION_RATIO {
            0 => 6,
            1 => 4,
            2 => 2,
            _ => unreachable!(),
        };
        if encoding == EBlastEncoding::Nucleotide {
            new_buffer[new_index as usize] |= (buffer[index as usize] & NCBI2NA_MASK) << shift;
        } else {
            new_buffer[new_index as usize] |=
                (NCBI4NA_TO_BLASTNA[buffer[index as usize] as usize] & NCBI2NA_MASK) << shift;
        }
        index += 1;
    }
    new_buffer
}

/// Matches C `BLAST_FrameToContext` (including subject frames for tblastn).
pub fn blast_frame_to_context(frame: i16, program: EBlastProgramType) -> i32 {
    if blast_query_is_translated(program) || blast_subject_is_translated(program) {
        assert!((-3..=3).contains(&frame) && frame != 0);
        if frame > 0 {
            i32::from(frame) - 1
        } else {
            2 - i32::from(frame)
        }
    } else if blast_query_is_nucleotide(program) || blast_subject_is_nucleotide(program) {
        assert!(frame == 1 || frame == -1);
        if frame == 1 {
            0
        } else {
            1
        }
    } else {
        assert_eq!(frame, 0);
        0
    }
}

/// Matches C `BLAST_GetNumberOfContexts` (query contexts, not subject frames).
pub fn blast_get_number_of_contexts(p: EBlastProgramType) -> u32 {
    if blast_query_is_translated(p) {
        NUM_FRAMES
    } else if blast_query_is_nucleotide(p) {
        NUM_STRANDS
    } else if blast_program_is_valid(p) {
        1
    } else {
        0
    }
}

/// Matches C `BLAST_GetTranslatedProteinLength`.
pub fn blast_get_translated_protein_length(nucleotide_length: usize, context: u32) -> usize {
    if nucleotide_length == 0 || nucleotide_length <= (context % CODON_LENGTH as u32) as usize {
        return 0;
    }
    (nucleotide_length - (context % CODON_LENGTH as u32) as usize) / CODON_LENGTH as usize
}

/// Matches C `s_BlastGetTranslationTable` in NCBI2na codon order.
/// An owned fixed-size array replaces the C allocated array and null output.
pub fn s_blast_get_translation_table(
    genetic_code: Option<&[u8; 64]>,
    reverse_complement: bool,
) -> Option<[u8; 64]> {
    let mapping: [i16; 4] = [2, 1, 3, 0];
    let genetic_code = genetic_code?;
    let mut translation = [0; 64];
    for index1 in 0..4 {
        for index2 in 0..4 {
            for index3 in 0..4 {
                if reverse_complement {
                    let bp1 = 3 - index1;
                    let bp2 = 3 - index2;
                    let bp3 = 3 - index3;
                    let codon = (mapping[bp1] << 4) + (mapping[bp2] << 2) + mapping[bp3];
                    translation[(index3 << 4) + (index2 << 2) + index1] =
                        genetic_code[codon as usize];
                } else {
                    let codon = (mapping[index1] << 4) + (mapping[index2] << 2) + mapping[index3];
                    translation[(index1 << 4) + (index2 << 2) + index3] =
                        genetic_code[codon as usize];
                }
            }
        }
    }
    Some(translation)
}

/// Matches C `BLAST_TranslateCompressedSequence`, retaining the unrolled loops.
///
/// `nt_seq` uses NCBI2na packing. Each strand needs the corresponding table from
/// `s_blast_get_translation_table`. Output includes two protein sentinels, except
/// that the original early return leaves the output entirely untouched.
///
/// Inputs entering the state machine must have at least four bases. The upstream
/// three-base case can read before `nt_seq` and is outside this kernel's supported
/// domain; short subjects must follow the search engine's short-subject handling.
pub fn blast_translate_compressed_sequence(
    translation: &[u8; 64],
    length: i32,
    nt_seq: Option<&[u8]>,
    frame: i16,
    prot_seq: Option<&mut [u8]>,
) -> i32 {
    let (nt_seq, prot_seq) = match (nt_seq, prot_seq) {
        (Some(nt_seq), Some(prot_seq)) if length - i32::from(frame).abs() + 1 >= CODON_LENGTH => {
            (nt_seq, prot_seq)
        }
        _ => return 0,
    };
    assert!((-3..=3).contains(&frame) && frame != 0);
    assert!(length >= 4);
    assert!(nt_seq.len() > length as usize / 4);
    let mut nt: isize = 0;
    prot_seq[0] = 0;
    let mut prot: usize = 1;
    let prot_seq_start = prot;
    let remainder = (length % 4) as u8;
    let mut state: i32;
    let mut codon: i32 = -1;
    let mut byte_value: i32;

    if frame > 0 {
        let nt_seq_end = (length / 4 - 1) as isize;
        let last_remainder = ((4 * (length / 4) - i32::from(frame) + 1) % CODON_LENGTH) as u8;
        let total_remainder = i16::from(last_remainder) + i16::from(remainder);
        state = i32::from(frame) - 1;
        byte_value = i32::from(nt_seq[nt as usize]);

        while nt < nt_seq_end {
            // C switch fallthrough is made explicit, without moving work between cases.
            let mut dispatch = state;
            loop {
                match dispatch {
                    0 => {
                        codon = byte_value >> 2;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        dispatch = 3;
                    }
                    3 => {
                        codon = (byte_value & 3) << 4;
                        nt += 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        codon += byte_value >> 4;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        if nt >= nt_seq_end {
                            state = 2;
                            break;
                        }
                        dispatch = 2;
                    }
                    2 => {
                        codon = (byte_value & 15) << 2;
                        nt += 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        codon += byte_value >> 6;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        if nt >= nt_seq_end {
                            state = 1;
                            break;
                        }
                        dispatch = 1;
                    }
                    1 => {
                        codon = byte_value & 63;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        nt += 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        state = 0;
                        break;
                    }
                    _ => unreachable!(),
                }
            }

            // Original optimized loop: eight codons from six packed bytes.
            while nt < nt_seq_end - 10 {
                nt += 1;
                let byte_value1 = i32::from(nt_seq[nt as usize]);
                nt += 1;
                let byte_value2 = i32::from(nt_seq[nt as usize]);
                nt += 1;
                let byte_value3 = i32::from(nt_seq[nt as usize]);
                codon = byte_value >> 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = (byte_value & 3) << 4;
                codon += byte_value1 >> 4;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                nt += 1;
                let byte_value4 = i32::from(nt_seq[nt as usize]);
                codon = (byte_value1 & 15) << 2;
                codon += byte_value2 >> 6;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value2 & 63;
                nt += 1;
                let byte_value5 = i32::from(nt_seq[nt as usize]);
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value3 >> 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                nt += 1;
                byte_value = i32::from(nt_seq[nt as usize]);
                codon = (byte_value3 & 3) << 4;
                codon += byte_value4 >> 4;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = (byte_value4 & 15) << 2;
                codon += byte_value5 >> 6;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value5 & 63;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                state = 0;
            }
        }
        if state == 1 {
            byte_value = i32::from(nt_seq[nt as usize]);
            codon = byte_value & 63;
            state = 0;
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        } else if state == 0 {
            byte_value = i32::from(nt_seq[nt as usize]);
            codon = byte_value >> 2;
            state = 3;
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        }
        if total_remainder >= CODON_LENGTH as i16 {
            byte_value = i32::from(nt_seq[nt_seq_end as usize]);
            let last_byte = i32::from(nt_seq[(nt_seq_end + 1) as usize]);
            if state == 0 {
                codon = last_byte >> 2;
            } else if state == 2 {
                codon = (byte_value & 15) << 2;
                codon += last_byte >> 6;
            } else if state == 3 {
                codon = (byte_value & 3) << 4;
                codon += last_byte >> 4;
            }
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        }
    } else {
        let nt_seq_start = nt;
        nt += (length / 4) as isize;
        state = i32::from(remainder) + i32::from(frame);
        if state >= 0 {
            let last_byte = i32::from(nt_seq[nt as usize]);
            nt -= 1;
            if state == 0 {
                codon = last_byte >> 6;
                byte_value = i32::from(nt_seq[nt as usize]);
                codon += (byte_value & 15) << 2;
                state = 1;
            } else if state == 1 {
                codon = last_byte >> 4;
                byte_value = i32::from(nt_seq[nt as usize]);
                codon += (byte_value & 3) << 4;
                state = 2;
            } else if state == 2 {
                codon = last_byte >> 2;
                state = 3;
            }
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        } else {
            state = 3 + (i32::from(remainder) + i32::from(frame) + 1);
            nt -= 1;
        }
        byte_value = i32::from(nt_seq[nt as usize]);
        while nt > nt_seq_start {
            let mut dispatch = state;
            loop {
                match dispatch {
                    3 => {
                        codon = byte_value & 63;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        dispatch = 0;
                    }
                    0 => {
                        codon = byte_value >> 6;
                        nt -= 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        codon += (byte_value & 15) << 2;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        if nt <= nt_seq_start {
                            state = 1;
                            break;
                        }
                        dispatch = 1;
                    }
                    1 => {
                        codon = byte_value >> 4;
                        nt -= 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        codon += (byte_value & 3) << 4;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        if nt <= nt_seq_start {
                            state = 2;
                            break;
                        }
                        dispatch = 2;
                    }
                    2 => {
                        codon = byte_value >> 2;
                        prot_seq[prot] = translation[codon as usize];
                        prot += 1;
                        nt -= 1;
                        byte_value = i32::from(nt_seq[nt as usize]);
                        state = 3;
                        break;
                    }
                    _ => unreachable!(),
                }
            }
            while nt > nt_seq_start + 10 {
                nt -= 1;
                let byte_value1 = i32::from(nt_seq[nt as usize]);
                nt -= 1;
                let byte_value2 = i32::from(nt_seq[nt as usize]);
                nt -= 1;
                let byte_value3 = i32::from(nt_seq[nt as usize]);
                codon = byte_value & 63;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value >> 6;
                codon += (byte_value1 & 15) << 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                nt -= 1;
                let byte_value4 = i32::from(nt_seq[nt as usize]);
                codon = byte_value1 >> 4;
                codon += (byte_value2 & 3) << 4;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value2 >> 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                nt -= 1;
                let byte_value5 = i32::from(nt_seq[nt as usize]);
                codon = byte_value3 & 63;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                nt -= 1;
                byte_value = i32::from(nt_seq[nt as usize]);
                codon = byte_value3 >> 6;
                codon += (byte_value4 & 15) << 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value4 >> 4;
                codon += (byte_value5 & 3) << 4;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
                codon = byte_value5 >> 2;
                prot_seq[prot] = translation[codon as usize];
                prot += 1;
            }
        }
        byte_value = i32::from(nt_seq[nt as usize]);
        if state == 3 {
            codon = byte_value & 63;
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        } else if state == 2 {
            codon = byte_value >> 2;
            prot_seq[prot] = translation[codon as usize];
            prot += 1;
        }
    }
    prot_seq[prot] = 0;
    (prot - prot_seq_start) as i32
}

/// Matches C `BLAST_GetAllTranslations` for NCBI2na and NCBI4na input.
/// Owned Vec outputs replace C allocated output pointers. Unused allocation tails
/// are zeroed; frame offsets delimit meaningful bytes exactly as in upstream.
#[allow(clippy::too_many_arguments)]
pub fn blast_get_all_translations(
    nucl_seq: &[u8],
    encoding: EBlastEncoding,
    nucl_length: i32,
    genetic_code: &[u8; 64],
    translation_buffer_ptr: Option<&mut Vec<u8>>,
    frame_offsets_ptr: Option<&mut [u32; 7]>,
    mixed_seq_ptr: Option<&mut Vec<u8>>,
) -> i16 {
    if encoding != EBlastEncoding::Ncbi2na && encoding != EBlastEncoding::Ncbi4na {
        return -1;
    }
    // Very short subjects must use upstream search-engine filtering before this
    // kernel; malloc-sized buffers in the original are not valid for every case.
    assert!(nucl_length >= 4);
    let buffer_length = (2 * (nucl_length + 1) + 2) as u32;
    let mut translation_buffer = vec![0; buffer_length as usize];
    let mut translation_table = None;
    let mut translation_table_rc = None;
    let mut nucl_seq_rev = Vec::new();
    if encoding == EBlastEncoding::Ncbi4na {
        nucl_seq_rev = get_reverse_nucl_sequence(nucl_seq, nucl_length);
    } else {
        translation_table = s_blast_get_translation_table(Some(genetic_code), false);
        translation_table_rc = s_blast_get_translation_table(Some(genetic_code), true);
    }
    let mut frame_offsets = [0; 7];
    let mut offset: u32 = 0;
    for context in 0..NUM_FRAMES {
        let frame = i16::from(blast_context_to_frame(EBlastProgramType::Blastx, context));
        let length: u32;
        if encoding == EBlastEncoding::Ncbi2na {
            if frame > 0 {
                length = blast_translate_compressed_sequence(
                    translation_table.as_ref().unwrap(),
                    nucl_length,
                    Some(nucl_seq),
                    frame,
                    Some(&mut translation_buffer[offset as usize..]),
                ) as u32;
            } else {
                length = blast_translate_compressed_sequence(
                    translation_table_rc.as_ref().unwrap(),
                    nucl_length,
                    Some(nucl_seq),
                    frame,
                    Some(&mut translation_buffer[offset as usize..]),
                ) as u32;
            }
        } else {
            length = blast_get_translation(
                nucl_seq,
                &nucl_seq_rev,
                nucl_length,
                frame,
                &mut translation_buffer[offset as usize..],
                genetic_code,
            ) as u32;
        }
        offset += length + 1;
        frame_offsets[context as usize + 1] = offset;
    }
    if let Some(mixed_seq_ptr) = mixed_seq_ptr {
        let mut mixed_seq = vec![0; (2 * nucl_length + 3) as usize];
        let mut seq = 0;
        let mut index = 0;
        while index < NUM_FRAMES {
            for i in 0..=nucl_length {
                let context = i % CODON_LENGTH;
                let offset = i / CODON_LENGTH;
                mixed_seq[seq] = translation_buffer
                    [frame_offsets[(index + context as u32) as usize] as usize + offset as usize];
                seq += 1;
            }
            index += CODON_LENGTH as u32;
        }
        mixed_seq[seq] = 0;
        *mixed_seq_ptr = mixed_seq;
    }
    if let Some(translation_buffer_ptr) = translation_buffer_ptr {
        *translation_buffer_ptr = translation_buffer;
    }
    if let Some(frame_offsets_ptr) = frame_offsets_ptr {
        *frame_offsets_ptr = frame_offsets;
    }
    0
}

/// Matches C `Blast_GetPartialTranslation`, used for selected subject ranges
/// during traceback and composition-adjusted redo alignment.
///
/// C writes only `nucl_length + 1` bytes to mixed output. The final allocation
/// byte is not part of the defined output.
pub fn blast_get_partial_translation(
    nucl_seq: &[u8],
    nucl_length: i32,
    frame: i16,
    genetic_code: &[u8; 64],
    translation_buffer_ptr: Option<&mut Vec<u8>>,
    protein_length: Option<&mut i32>,
    mixed_seq_ptr: Option<&mut Vec<u8>>,
) -> i32 {
    assert!(nucl_length >= 0);
    let mut nucl_seq_rev = Vec::new();
    if frame < 0 {
        nucl_seq_rev = get_reverse_nucl_sequence(nucl_seq, nucl_length);
    }
    let mut translation_buffer;
    if mixed_seq_ptr.is_none() {
        translation_buffer = vec![0; (nucl_length / CODON_LENGTH + 2) as usize];
        let length = blast_get_translation(
            nucl_seq,
            &nucl_seq_rev,
            nucl_length,
            frame,
            &mut translation_buffer,
            genetic_code,
        );
        if let Some(protein_length) = protein_length {
            *protein_length = length;
        }
    } else {
        // Original allocation cannot hold all sentinels for lengths 0 or 1.
        assert!(nucl_length >= 2);
        let frame_sign: i16 = if frame < 0 { -1 } else { 1 };
        let mut offset: u32 = 0;
        let mut frame_offsets = [0; 3];
        translation_buffer = vec![0; (nucl_length + 2) as usize];
        for index in 1..=CODON_LENGTH as i16 {
            let length = blast_get_translation(
                nucl_seq,
                &nucl_seq_rev,
                nucl_length,
                frame_sign * index,
                &mut translation_buffer[offset as usize..],
                genetic_code,
            );
            frame_offsets[index as usize - 1] = offset;
            offset += length as u32 + 1;
        }
        let mut mixed_seq = vec![0; (nucl_length + 2) as usize];
        if let Some(protein_length) = protein_length {
            *protein_length = nucl_length;
        }
        for index in 0..=nucl_length {
            mixed_seq[index as usize] =
                translation_buffer[frame_offsets[(index % CODON_LENGTH) as usize] as usize
                    + (index / CODON_LENGTH) as usize];
        }
        *mixed_seq_ptr.unwrap() = mixed_seq;
    }
    if let Some(translation_buffer_ptr) = translation_buffer_ptr {
        *translation_buffer_ptr = translation_buffer;
    }
    0
}
