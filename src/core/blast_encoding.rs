// NCBI BLAST+ 2.17.0+, blast_encoding.c and blast_encoding.h.
// Author: Christiam Camacho. Public domain notice: ../../LICENSE-NCBI.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u8)]
pub enum EBlastEncoding {
    Protein = 0,
    Nucleotide = 1,
    Ncbi4na = 2,
    Ncbi2na = 3,
    Error = 255,
}

pub const BLAST2NA_SIZE: usize = 4;
pub const BLASTNA_SIZE: usize = 16;
pub const BLASTAA_SIZE: usize = 28;
pub const BLASTNA_SEQ_CODE: usize = 99;
pub const BLASTAA_SEQ_CODE: usize = 11;
pub const NCBI4NA_SEQ_CODE: usize = 4;

pub const NCBI4NA_TO_BLASTNA: [u8; 16] = [15, 0, 1, 6, 2, 4, 9, 13, 3, 8, 5, 12, 7, 11, 10, 14];

pub const BLASTNA_TO_NCBI4NA: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];

pub const BLASTNA_TO_IUPACNA: [u8; 16] = [
    b'A', b'C', b'G', b'T', b'R', b'Y', b'M', b'K', b'W', b'S', b'B', b'D', b'H', b'V', b'N', b'-',
];

pub const NCBI4NA_TO_IUPACNA: [u8; 16] = [
    b'-', b'A', b'C', b'M', b'G', b'R', b'S', b'V', b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N',
];

pub const IUPACNA_TO_BLASTNA: [u8; 128] = [
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0, 10, 1, 11, 15, 15, 2,
    12, 15, 15, 7, 15, 6, 14, 15, 15, 15, 4, 9, 3, 15, 13, 8, 15, 5, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 15, 15, 15,
];

pub const IUPACNA_TO_NCBI4NA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 14, 2, 13, 0, 0, 4, 11, 0, 0, 12, 0, 3, 15, 0, 0, 0, 5, 6, 8, 0, 7, 9, 0, 10, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0,
];

pub const AMINOACID_TO_NCBISTDAA: [u8; 128] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 27, 10, 11, 12, 13, 26, 14, 15, 16, 17, 18, 24, 19, 20, 21,
    22, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
];

pub const NCBISTDAA_TO_AMINOACID: [u8; 28] = [
    b'-', b'A', b'B', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q',
    b'R', b'S', b'T', b'V', b'W', b'X', b'Y', b'Z', b'U', b'*', b'O', b'J',
];

pub const K_PROT_SENTINEL: u8 = 0;
pub const K_NUCL_SENTINEL: u8 = 0xF;
