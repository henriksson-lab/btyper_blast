//! Owned in-memory nucleotide storage for the BLAST sequence-source boundary.
//! These are library input/storage adapters, separate from translated kernels.
//!
//! ```
//! use btyper_blast::DatabaseBuilder;
//!
//! let mut builder = DatabaseBuilder::new();
//! builder.add_sequence("contig_1", b"ACGTRYSWKMN").unwrap();
//! let database = builder.build().unwrap();
//! assert_eq!(database.total_length(), 11);
//! assert_eq!(database.sequence(0).unwrap().id(), "contig_1");
//! ```

use crate::core::blast_encoding::{EBlastEncoding, IUPACNA_TO_NCBI4NA, NCBI4NA_TO_BLASTNA};
use crate::core::blast_util::blast_pack_dna;
use crate::seqdb::writedb_convert::write_db_ncbi4na_to_binary;
use std::fmt;
use std::io::BufRead;

#[derive(Debug)]
pub enum DatabaseError {
    EmptyDatabase,
    EmptyIdentifier,
    EmptySequence,
    SequenceTooLong,
    TooManySequences,
    TotalLengthOverflow,
    InvalidResidue { offset: usize, byte: u8 },
    InvalidNcbi4na { offset: usize, value: u8 },
    MissingFastaHeader { line: usize },
    Io(std::io::Error),
}

impl fmt::Display for DatabaseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDatabase => write!(f, "database contains no sequences"),
            Self::EmptyIdentifier => write!(f, "sequence identifier is empty"),
            Self::EmptySequence => write!(f, "sequence contains no residues"),
            Self::SequenceTooLong => write!(f, "sequence exceeds the BLAST i32 length domain"),
            Self::TooManySequences => {
                write!(f, "sequence count exceeds the BLAST i32 ordinal domain")
            }
            Self::TotalLengthOverflow => {
                write!(f, "total sequence length exceeds the BLAST i64 domain")
            }
            Self::InvalidResidue { offset, byte } => {
                write!(f, "invalid nucleotide byte {byte:#04x} at offset {offset}")
            }
            Self::InvalidNcbi4na { offset, value } => {
                write!(f, "invalid NCBI4na value {value} at offset {offset}")
            }
            Self::MissingFastaHeader { line } => write!(
                f,
                "sequence data before the first FASTA header on line {line}"
            ),
            Self::Io(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for DatabaseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<std::io::Error> for DatabaseError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// An invalid FASTA residue omitted from the sequence, with its input location.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ResidueWarning {
    pub sequence_ordinal: usize,
    pub line: usize,
    pub column: usize,
    pub byte: u8,
}

/// A nucleotide subject. Both encodings exclude leading sentinels; packed data
/// includes the original final length-remainder byte. Ambiguities are restored
/// in the expanded NCBI4na data and retained in encoded database form separately.
#[derive(Debug)]
pub struct NucleotideSequence {
    id: String,
    ncbi4na: Vec<u8>,
    packed_ncbi2na: Vec<u8>,
    ambiguity_data: Vec<u8>,
}

impl NucleotideSequence {
    pub fn id(&self) -> &str {
        &self.id
    }
    pub fn len(&self) -> usize {
        self.ncbi4na.len()
    }
    pub fn is_empty(&self) -> bool {
        self.ncbi4na.is_empty()
    }
    pub fn ncbi4na(&self) -> &[u8] {
        &self.ncbi4na
    }
    pub fn packed_ncbi2na(&self) -> &[u8] {
        &self.packed_ncbi2na[..self.len() / 4 + 1]
    }
    pub(crate) fn packed_for_scanning(&self) -> &[u8] {
        &self.packed_ncbi2na
    }
    pub fn ambiguity_data(&self) -> &[u8] {
        &self.ambiguity_data
    }

    /// BLASTNA residues without leading or trailing sentinels.
    pub fn blastna(&self) -> Vec<u8> {
        self.ncbi4na
            .iter()
            .map(|base| NCBI4NA_TO_BLASTNA[*base as usize])
            .collect()
    }
}

#[derive(Debug)]
pub struct NucleotideDatabase {
    sequences: Vec<NucleotideSequence>,
    total_length: u64,
    max_length: usize,
    warnings: Vec<ResidueWarning>,
}

impl NucleotideDatabase {
    pub fn len(&self) -> usize {
        self.sequences.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
    pub fn total_length(&self) -> u64 {
        self.total_length
    }
    pub fn max_length(&self) -> usize {
        self.max_length
    }
    pub fn sequences(&self) -> &[NucleotideSequence] {
        &self.sequences
    }
    pub fn sequence(&self, ordinal: usize) -> Option<&NucleotideSequence> {
        self.sequences.get(ordinal)
    }
    pub fn warnings(&self) -> &[ResidueWarning] {
        &self.warnings
    }

    /// Read ordinary nucleotide FASTA. Identifiers are the first header token.
    /// U becomes T and '-' becomes N, matching the default NCBI nucleotide reader.
    /// Invalid sequence bytes are omitted with warnings returned in the database;
    /// nothing is printed. Empty records are rejected rather than silently lost.
    pub fn from_fasta<R: BufRead>(reader: R) -> Result<Self, DatabaseError> {
        let mut builder = DatabaseBuilder::new();
        let mut id = None;
        let mut sequence = Vec::new();
        let mut warnings = Vec::new();
        for (line_index, line) in reader.split(b'\n').enumerate() {
            let line = line?;
            if line.first() == Some(&b'>') {
                if let Some(previous_id) = id.take() {
                    builder.add_ncbi4na(previous_id, &sequence)?;
                    sequence.clear();
                }
                let header = std::str::from_utf8(&line[1..])
                    .map_err(|error| std::io::Error::new(std::io::ErrorKind::InvalidData, error))?;
                id = Some(
                    header
                        .split_whitespace()
                        .next()
                        .ok_or(DatabaseError::EmptyIdentifier)?
                        .to_owned(),
                );
                continue;
            }
            if line.first() == Some(&b';') || line.iter().all(|byte| byte.is_ascii_whitespace()) {
                continue;
            }
            if id.is_none() {
                return Err(DatabaseError::MissingFastaHeader {
                    line: line_index + 1,
                });
            }
            for (column, byte) in line.into_iter().enumerate() {
                if byte.is_ascii_whitespace() {
                    continue;
                }
                if let Some(base) = nucleotide_value(byte) {
                    sequence.push(base);
                } else {
                    warnings.push(ResidueWarning {
                        sequence_ordinal: builder.sequences.len(),
                        line: line_index + 1,
                        column: column + 1,
                        byte,
                    });
                }
            }
        }
        if let Some(id) = id {
            builder.add_ncbi4na(id, &sequence)?;
        }
        let mut database = builder.build()?;
        database.warnings = warnings;
        Ok(database)
    }
}

#[derive(Debug, Default)]
pub struct DatabaseBuilder {
    sequences: Vec<NucleotideSequence>,
    total_length: u64,
    max_length: usize,
}

impl DatabaseBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add ASCII IUPAC nucleotide residues. This direct API rejects invalid
    /// bytes; use from_fasta for whitespace and invalid-residue diagnostics.
    pub fn add_sequence(
        &mut self,
        id: impl Into<String>,
        sequence: &[u8],
    ) -> Result<(), DatabaseError> {
        let mut encoded = Vec::with_capacity(sequence.len());
        for (offset, byte) in sequence.iter().enumerate() {
            encoded.push(
                nucleotide_value(*byte).ok_or(DatabaseError::InvalidResidue {
                    offset,
                    byte: *byte,
                })?,
            );
        }
        self.add_ncbi4na(id, &encoded)
    }

    /// Add one expanded NCBI4na value per base (not packed nibbles). Raw gap
    /// code zero becomes N, as in the upstream database conversion. Ordinals and
    /// duplicate identifiers are retained in insertion order.
    pub fn add_ncbi4na(
        &mut self,
        id: impl Into<String>,
        sequence: &[u8],
    ) -> Result<(), DatabaseError> {
        let id = id.into();
        if id.trim().is_empty() {
            return Err(DatabaseError::EmptyIdentifier);
        }
        if sequence.is_empty() {
            return Err(DatabaseError::EmptySequence);
        }
        if sequence.len() > i32::MAX as usize - 1 {
            return Err(DatabaseError::SequenceTooLong);
        }
        if self.sequences.len() >= i32::MAX as usize {
            return Err(DatabaseError::TooManySequences);
        }
        let total_length = self
            .total_length
            .checked_add(sequence.len() as u64)
            .filter(|value| *value <= i64::MAX as u64)
            .ok_or(DatabaseError::TotalLengthOverflow)?;
        let mut ncbi4na = Vec::with_capacity(sequence.len());
        let mut unambiguous = true;
        for (offset, value) in sequence.iter().copied().enumerate() {
            if value > 15 {
                return Err(DatabaseError::InvalidNcbi4na { offset, value });
            }
            let value = if value == 0 { 15 } else { value };
            if value != 1 && value != 2 && value != 4 && value != 8 {
                unambiguous = false;
            }
            ncbi4na.push(value);
        }
        let (mut packed_ncbi2na, mut ambiguity_data) = (Vec::new(), Vec::new());
        // CSeqConvert_imp::CPacker chooses the wide (4na) encoding on a byte-cost
        // tie. For unambiguous records this occurs at lengths 1 and 2. Preserve
        // that choice because it changes whether an empty ambiguity header is
        // written, even though the biological sequence is the same.
        let prefer_ncbi2na = (ncbi4na.len() + 3) / 4 < (ncbi4na.len() + 1) / 2;
        if unambiguous && prefer_ncbi2na {
            // NCBI's input conversion stores unambiguous sequences as NCBI2na,
            // bypassing the ambiguity builder and its empty four-byte header.
            packed_ncbi2na =
                blast_pack_dna(&ncbi4na, ncbi4na.len() as i32, EBlastEncoding::Ncbi4na);
        } else {
            let mut packed_ncbi4na = vec![0; (ncbi4na.len() + 1) / 2];
            for (i, value) in ncbi4na.iter().enumerate() {
                packed_ncbi4na[i / 2] |= value << if i % 2 == 0 { 4 } else { 0 };
            }
            write_db_ncbi4na_to_binary(
                &packed_ncbi4na,
                packed_ncbi4na.len() as i32,
                ncbi4na.len() as i32,
                &mut packed_ncbi2na,
                &mut ambiguity_data,
            );
        }
        // The original word scanners fetch packed windows beyond the final
        // byte containing residues. Keep readable tail bytes in memory while
        // exposing the exact original database payload through its getter.
        packed_ncbi2na.resize(packed_ncbi2na.len() + 3, 0);
        self.max_length = self.max_length.max(ncbi4na.len());
        self.sequences.push(NucleotideSequence {
            id,
            ncbi4na,
            packed_ncbi2na,
            ambiguity_data,
        });
        self.total_length = total_length;
        Ok(())
    }

    pub fn build(self) -> Result<NucleotideDatabase, DatabaseError> {
        if self.sequences.is_empty() {
            return Err(DatabaseError::EmptyDatabase);
        }
        Ok(NucleotideDatabase {
            sequences: self.sequences,
            total_length: self.total_length,
            max_length: self.max_length,
            warnings: Vec::new(),
        })
    }
}

pub(crate) fn nucleotide_value(byte: u8) -> Option<u8> {
    match byte.to_ascii_uppercase() {
        b'U' => Some(8),
        b'-' => Some(15),
        upper => IUPACNA_TO_NCBI4NA
            .get(upper as usize)
            .copied()
            .filter(|value| *value != 0),
    }
}
