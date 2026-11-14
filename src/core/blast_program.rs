// NCBI BLAST+ 2.17.0+, blast_program.c and blast_program.h.
// Authors: Christiam Camacho / Ilya Dondoshansky. See LICENSE-NCBI.

const PROTEIN_QUERY_MASK: u32 = 1 << 0;
const PROTEIN_SUBJECT_MASK: u32 = 1 << 1;
const NUCLEOTIDE_QUERY_MASK: u32 = 1 << 2;
const NUCLEOTIDE_SUBJECT_MASK: u32 = 1 << 3;
const TRANSLATED_QUERY_MASK: u32 = 1 << 4;
const TRANSLATED_SUBJECT_MASK: u32 = 1 << 5;
const MAPPING_MASK: u32 = 1 << 9;

/// Original discriminants are retained even for unsupported search programs.
/// A discriminant here does not imply a public search implementation.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
#[repr(u32)]
pub enum EBlastProgramType {
    Undefined = 0,
    Blastp = 3,
    Blastn = 12,
    Blastx = 22,
    Tblastn = 41,
    Tblastx = 60,
    PsiBlast = 67,
    PsiTblastn = 105,
    RpsBlast = 131,
    RpsTblastn = 150,
    PhiBlastp = 259,
    PhiBlastn = 268,
    Mapping = 524,
}

/// Matches C `Blast_QueryIsProtein`.
pub fn blast_query_is_protein(p: EBlastProgramType) -> bool {
    p as u32 & PROTEIN_QUERY_MASK != 0
}

/// Matches C `Blast_QueryIsNucleotide`.
pub fn blast_query_is_nucleotide(p: EBlastProgramType) -> bool {
    p as u32 & NUCLEOTIDE_QUERY_MASK != 0
}

/// Matches C `Blast_SubjectIsProtein`.
pub fn blast_subject_is_protein(p: EBlastProgramType) -> bool {
    p as u32 & PROTEIN_SUBJECT_MASK != 0
}

/// Matches C `Blast_SubjectIsNucleotide`.
pub fn blast_subject_is_nucleotide(p: EBlastProgramType) -> bool {
    p as u32 & NUCLEOTIDE_SUBJECT_MASK != 0
}

/// Matches C `Blast_QueryIsTranslated`.
pub fn blast_query_is_translated(p: EBlastProgramType) -> bool {
    p as u32 & TRANSLATED_QUERY_MASK != 0
}

/// Matches C `Blast_SubjectIsTranslated`.
pub fn blast_subject_is_translated(p: EBlastProgramType) -> bool {
    p as u32 & TRANSLATED_SUBJECT_MASK != 0
}

/// Matches C `Blast_ProgramIsMapping`.
pub fn blast_program_is_mapping(p: EBlastProgramType) -> bool {
    p as u32 & MAPPING_MASK != 0
}

/// Matches C `Blast_ProgramIsNucleotide`.
pub fn blast_program_is_nucleotide(p: EBlastProgramType) -> bool {
    blast_query_is_nucleotide(p)
        && blast_subject_is_nucleotide(p)
        && !blast_query_is_translated(p)
        && !blast_subject_is_translated(p)
}

/// Matches C `Blast_ProgramIsValid`.
pub fn blast_program_is_valid(p: EBlastProgramType) -> bool {
    match p {
        EBlastProgramType::Blastp
        | EBlastProgramType::Blastn
        | EBlastProgramType::Blastx
        | EBlastProgramType::Tblastn
        | EBlastProgramType::Tblastx
        | EBlastProgramType::PsiBlast
        | EBlastProgramType::PsiTblastn
        | EBlastProgramType::RpsBlast
        | EBlastProgramType::RpsTblastn
        | EBlastProgramType::PhiBlastp
        | EBlastProgramType::PhiBlastn
        | EBlastProgramType::Mapping => true,
        EBlastProgramType::Undefined => false,
    }
}

pub fn blast_query_is_pssm(p: EBlastProgramType) -> bool {
    p as u32 & (1 << 6) != 0
}
pub fn blast_program_is_phi_blast(p: EBlastProgramType) -> bool {
    p as u32 & (1 << 8) != 0
}
pub fn blast_query_is_pattern(p: EBlastProgramType) -> bool {
    p as u32 & (1 << 8) != 0
}

pub fn blast_program_is_rps_blast(p: EBlastProgramType) -> bool {
    p as u32 & (1 << 7) != 0
}
