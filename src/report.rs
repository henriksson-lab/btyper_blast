//! In-memory display projection of completed native HSPs. This adapter does not
//! alter search algorithms or scores and does not construct ASN.1 Seq-aligns.
use crate::core::{
    blast_encoding::{BLASTNA_TO_IUPACNA, NCBISTDAA_TO_AMINOACID},
    blast_hits::BlastHSP,
    gapinfo::EGapAlignOpType,
};
use std::fmt;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReportError {
    InvalidCoordinates,
    MissingTraceback,
    InvalidEditScript,
    UnsupportedFrame,
    InvalidResidue,
}
impl fmt::Display for ReportError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "cannot report alignment: {self:?}")
    }
}
impl std::error::Error for ReportError {}

/// Endpoints use BLAST tabular's one-based, inclusive coordinates. Nucleotide
/// query alignments are oriented to the forward query; tblastn retains the
/// subject translation frame. Raw scores and floating values are unchanged.
#[derive(Debug)]
pub struct AlignmentReport {
    pub qstart: i32,
    pub qend: i32,
    pub sstart: i32,
    pub send: i32,
    pub qframe: i16,
    pub sframe: i16,
    pub length: usize,
    pub nident: i32,
    pub mismatch: usize,
    pub gaps: usize,
    pub gapopen: usize,
    pub score: i32,
    pub bitscore: f64,
    pub evalue: f64,
    pub btop: String,
    pub qseq: String,
    pub sseq: String,
}
impl AlignmentReport {
    /// `query_context` is the original context selected by hsp.context, without
    /// a leading sentinel. `subject` contains restored, uncompressed BLASTna.
    pub fn nucleotide(
        hsp: &BlastHSP,
        query_context: &[u8],
        subject: &[u8],
    ) -> Result<Self, ReportError> {
        Self::nucleotide_range(hsp, query_context, subject, 0)
    }

    /// Display only the restored subject interval covered by this HSP. The
    /// subject origin is absolute; all stored and reported coordinates retain
    /// their original full-subject values.
    pub fn nucleotide_range(
        hsp: &BlastHSP,
        query_context: &[u8],
        subject: &[u8],
        subject_origin: i32,
    ) -> Result<Self, ReportError> {
        let (mut qseq, mut sseq, gapopen) = columns(
            hsp,
            query_context,
            subject,
            subject_origin,
            &BLASTNA_TO_IUPACNA,
        )?;
        let reverse = hsp.query.frame < 0;
        let (qstart, qend, sstart, send) = if reverse {
            reverse_complement_columns(&mut qseq);
            reverse_complement_columns(&mut sseq);
            (
                query_context.len() as i32 - hsp.query.end + 1,
                query_context.len() as i32 - hsp.query.offset,
                hsp.subject.end,
                hsp.subject.offset + 1,
            )
        } else {
            (
                hsp.query.offset + 1,
                hsp.query.end,
                hsp.subject.offset + 1,
                hsp.subject.end,
            )
        };
        Ok(Self::assemble(
            hsp,
            qseq,
            sseq,
            gapopen,
            qstart,
            qend,
            sstart,
            send,
            1,
            if reverse { -1 } else { 1 },
        ))
    }

    /// Both input buffers are unmasked NCBIstdaa; `subject_frame` is the full
    /// translated frame selected by the HSP, and `nucleotide_length` is the
    /// original subject length. Translation is shared with the search core.
    pub fn translated_subject(
        hsp: &BlastHSP,
        query: &[u8],
        subject_frame: &[u8],
        nucleotide_length: i32,
    ) -> Result<Self, ReportError> {
        Self::translated_subject_range(hsp, query, subject_frame, 0, nucleotide_length)
    }

    /// Project a cached partial target translation. `subject_origin` is its
    /// first residue's offset in the complete translated frame. HSP coordinates
    /// retain their original frame offsets throughout display construction.
    pub fn translated_subject_range(
        hsp: &BlastHSP,
        query: &[u8],
        subject_frame: &[u8],
        subject_origin: i32,
        nucleotide_length: i32,
    ) -> Result<Self, ReportError> {
        let frame = hsp.subject.frame as i32;
        if frame == 0 || frame.abs() > 3 || hsp.query.frame != 0 {
            return Err(ReportError::UnsupportedFrame);
        }
        let (qseq, sseq, gapopen) = columns(
            hsp,
            query,
            subject_frame,
            subject_origin,
            &NCBISTDAA_TO_AMINOACID,
        )?;
        let (sstart, send) = if frame > 0 {
            (
                3 * hsp.subject.offset + frame,
                3 * hsp.subject.end + frame - 1,
            )
        } else {
            (
                nucleotide_length - 3 * hsp.subject.offset + frame + 1,
                nucleotide_length - 3 * hsp.subject.end + frame + 2,
            )
        };
        Ok(Self::assemble(
            hsp,
            qseq,
            sseq,
            gapopen,
            hsp.query.offset + 1,
            hsp.query.end,
            sstart,
            send,
            0,
            hsp.subject.frame,
        ))
    }

    #[allow(clippy::too_many_arguments)]
    fn assemble(
        hsp: &BlastHSP,
        qseq: Vec<u8>,
        sseq: Vec<u8>,
        gapopen: usize,
        qstart: i32,
        qend: i32,
        sstart: i32,
        send: i32,
        qframe: i16,
        sframe: i16,
    ) -> Self {
        let mut btop = String::new();
        let (mut matches, mut gaps, mut mismatch, mut nident) = (0usize, 0usize, 0usize, 0i32);
        for (&q, &s) in qseq.iter().zip(&sseq) {
            if q == s {
                matches += 1;
                nident += 1;
            } else {
                if matches != 0 {
                    btop.push_str(&matches.to_string());
                    matches = 0;
                }
                btop.push(q as char);
                btop.push(s as char);
                if q == b'-' || s == b'-' {
                    gaps += 1;
                } else {
                    mismatch += 1;
                }
            }
        }
        if matches != 0 {
            btop.push_str(&matches.to_string());
        }
        Self {
            qstart,
            qend,
            sstart,
            send,
            qframe,
            sframe,
            length: qseq.len(),
            nident,
            mismatch,
            gaps,
            gapopen,
            score: hsp.score,
            bitscore: hsp.bit_score,
            evalue: hsp.evalue,
            btop,
            qseq: String::from_utf8(qseq).expect("ASCII residue table"),
            sseq: String::from_utf8(sseq).expect("ASCII residue table"),
        }
    }
}

fn columns(
    hsp: &BlastHSP,
    query: &[u8],
    subject: &[u8],
    subject_origin: i32,
    alphabet: &[u8],
) -> Result<(Vec<u8>, Vec<u8>, usize), ReportError> {
    let script = hsp.gap_info.as_ref().ok_or(ReportError::MissingTraceback)?;
    if hsp.query.offset < 0
        || subject_origin < 0
        || hsp.subject.offset < subject_origin
        || hsp.query.end < hsp.query.offset
        || hsp.subject.end < hsp.subject.offset
        || hsp.query.end as usize > query.len()
        || (hsp.subject.end - subject_origin) as usize > subject.len()
    {
        return Err(ReportError::InvalidCoordinates);
    }
    let (mut q, mut s) = (
        hsp.query.offset as usize,
        (hsp.subject.offset - subject_origin) as usize,
    );
    let (mut qseq, mut sseq) = (Vec::new(), Vec::new());
    let mut gapopen = 0;
    if script.size < 0
        || script.size as usize > script.num.len()
        || script.size as usize > script.op_type.len()
    {
        return Err(ReportError::InvalidEditScript);
    }
    for i in 0..script.size as usize {
        let (advance_q, advance_s) = match script.op_type[i] {
            EGapAlignOpType::Sub => (true, true),
            EGapAlignOpType::Ins => {
                gapopen += 1;
                (true, false)
            }
            EGapAlignOpType::Del => {
                gapopen += 1;
                (false, true)
            }
            _ => return Err(ReportError::InvalidEditScript),
        };
        if script.num[i] <= 0 {
            return Err(ReportError::InvalidEditScript);
        }
        for _ in 0..script.num[i] {
            let qc = if advance_q {
                let code = *query.get(q).ok_or(ReportError::InvalidCoordinates)?;
                q += 1;
                *alphabet
                    .get(code as usize)
                    .ok_or(ReportError::InvalidResidue)?
            } else {
                b'-'
            };
            let sc = if advance_s {
                let code = *subject.get(s).ok_or(ReportError::InvalidCoordinates)?;
                s += 1;
                *alphabet
                    .get(code as usize)
                    .ok_or(ReportError::InvalidResidue)?
            } else {
                b'-'
            };
            qseq.push(qc);
            sseq.push(sc);
        }
    }
    if q != hsp.query.end as usize || s != (hsp.subject.end - subject_origin) as usize {
        return Err(ReportError::InvalidEditScript);
    }
    Ok((qseq, sseq, gapopen))
}
fn reverse_complement_columns(columns: &mut [u8]) {
    columns.reverse();
    for base in columns {
        *base = match *base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'R' => b'Y',
            b'Y' => b'R',
            b'M' => b'K',
            b'K' => b'M',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            other => other,
        };
    }
}
