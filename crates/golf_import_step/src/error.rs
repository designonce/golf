use golf_step::DecodeError;
use golf_step::Id;
use golf_step::ParseError;

/// Why a file couldn't be imported at all.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum ImportError {
    #[error(transparent)]
    Parse(#[from] ParseError),
}

/// Something in a file that was skipped or approximated, and why.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[error("#{id}: {problem}")]
pub struct ImportWarning {
    /// The instance concerned.
    pub id: Id,
    pub problem: Problem,
}

/// What went wrong with an instance.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum Problem {
    #[error(transparent)]
    Decode(#[from] DecodeError),
    #[error("{0} isn't supported")]
    Unsupported(String),
    #[error("{0}")]
    Invalid(String),
    #[error("healing left a problem: {0}")]
    Heal(String),
}

/// A problem with instance `id`.
pub(crate) fn problem(id: Id, problem: impl Into<Problem>) -> ImportWarning {
    ImportWarning {
        id,
        problem: problem.into(),
    }
}

pub(crate) fn invalid(id: Id, message: impl Into<String>) -> ImportWarning {
    problem(id, Problem::Invalid(message.into()))
}

pub(crate) fn unsupported(id: Id, what: impl Into<String>) -> ImportWarning {
    problem(id, Problem::Unsupported(what.into()))
}

/// Results of reading one instance: a warning naming it on failure.
pub(crate) type Read<T> = Result<T, ImportWarning>;

/// Attaches `id` to a decoding error.
pub(crate) trait At<T> {
    fn at(self, id: Id) -> Read<T>;
}

impl<T> At<T> for Result<T, DecodeError> {
    fn at(self, id: Id) -> Read<T> {
        self.map_err(|e| problem(id, e))
    }
}
