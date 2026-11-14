//! A scoped translation of NCBI BLAST+ for nucleotide databases, blastn, and tblastn.
//!
//! The source investigation and function mapping live in `upstream/`.
//! Standard nucleotide and translated-subject search APIs assemble the pinned
//! kernels over in-memory databases. Results remain owned Rust HSP values;
//! `report` provides optional display projection. Versioned differential and
//! workflow benchmark evidence is recorded in `upstream/`; actual remote
//! cross-platform CI results remain unverified.
//! Database creation and search results are in memory; disk persistence is outside scope.

#![deny(unsafe_op_in_unsafe_fn)]

pub mod core;
pub mod database;
pub mod dustmask;
pub mod nucleotide_setup;
pub mod protein_query;
pub mod protein_setup;
pub mod query;
pub mod query_batch;
pub mod report;
pub mod seed;
pub mod seqdb;
pub mod util;

pub use database::{
    DatabaseBuilder, DatabaseError, NucleotideDatabase, NucleotideSequence, ResidueWarning,
};
pub use query::{NucleotideTask, PreparedNucleotideQuery, QueryError};
pub use query_batch::{PreparedNucleotideBatch, QueryBatchError};

/// The upstream release whose behavior the translation targets.
pub const UPSTREAM_BLAST_VERSION: &str = "2.17.0+";

pub use protein_query::{PreparedProteinBatch, PreparedProteinQuery, ProteinQueryError};
pub use protein_setup::{
    ProteinInitialSearch, ProteinPreliminaryFrame, ProteinPreliminarySubject, ProteinSetupError,
};

pub mod nucleotide_search;

mod protein_chunks;
pub mod protein_redo;
pub mod protein_search;

pub use nucleotide_search::{
    search_nucleotide_database, NucleotideSearchError, NucleotideSearchResults,
};
pub use protein_search::{search_protein_database, ProteinSearchError, ProteinSearchResults};
pub mod query_batch_size;
pub mod query_split;
