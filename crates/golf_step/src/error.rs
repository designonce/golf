use crate::value::Id;

/// Why a file couldn't be read at all.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum ParseError {
    #[error("not a STEP physical file: it doesn't start with ISO-10303-21")]
    NotPart21,
    #[error("line {line}: {message}")]
    Syntax { line: usize, message: String },
}

/// An instance that was skipped because it couldn't be parsed.
#[derive(Clone, Debug, PartialEq)]
pub struct ParseIssue {
    /// The instance's id, if that much could be read.
    pub id: Option<Id>,
    pub line: usize,
    pub message: String,
}

/// Why an entity couldn't be read as what it was expected to be.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum DecodeError {
    #[error("#{id} is missing")]
    Missing { id: Id },
    #[error("#{id} is a {found}, not a {expected}")]
    WrongEntity {
        id: Id,
        expected: String,
        found: String,
    },
    #[error("{entity} parameter {index}: expected {expected}, found {found}")]
    WrongValue {
        entity: String,
        index: usize,
        expected: &'static str,
        found: String,
    },
    #[error("{entity} has no parameter {index}")]
    MissingParameter { entity: String, index: usize },
}
