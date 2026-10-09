use std::collections::HashMap;

use crate::error::DecodeError;
use crate::error::ParseError;
use crate::error::ParseIssue;
use crate::parse::parse;
use crate::value::Entity;
use crate::value::Id;
use crate::value::Record;

/// A parsed STEP file: its header records and its entity instances by id.
#[derive(Clone, Debug, Default)]
pub struct StepData {
    /// `FILE_DESCRIPTION`, `FILE_NAME`, `FILE_SCHEMA` and any others.
    pub header: Vec<Record>,
    pub(crate) entities: HashMap<Id, Entity>,
    /// Instances skipped because they couldn't be parsed.
    pub issues: Vec<ParseIssue>,
}

impl StepData {
    /// Parses a file. Its text is read as UTF-8 where it is, and as Latin-1
    /// where it isn't, so no instance is lost to an encoding.
    pub fn parse(bytes: &[u8]) -> Result<Self, ParseError> {
        parse(bytes)
    }

    /// The instance `#id`.
    pub fn get(&self, id: Id) -> Option<&Entity> {
        self.entities.get(&id)
    }

    /// The instance `#id`, or an error if there's none.
    pub fn entity(&self, id: Id) -> Result<&Entity, DecodeError> {
        self.get(id).ok_or(DecodeError::Missing { id })
    }

    /// The record named `name` of instance `#id`, or an error if there's no
    /// such instance or it isn't one.
    pub fn record(&self, id: Id, name: &str) -> Result<&Record, DecodeError> {
        self.entity(id)?.expect(name)
    }

    /// Every instance, in no particular order.
    pub fn entities(&self) -> impl ExactSizeIterator<Item = &Entity> {
        self.entities.values()
    }

    /// Every instance that is, or has a part, named `name`, by increasing id.
    pub fn all(&self, name: &str) -> Vec<&Entity> {
        let mut found: Vec<&Entity> = self.entities.values().filter(|e| e.is(name)).collect();
        found.sort_by_key(|e| e.id);
        found
    }

    pub fn len(&self) -> usize {
        self.entities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entities.is_empty()
    }

    /// The schema names from `FILE_SCHEMA`, e.g. `AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF`.
    pub fn schemas(&self) -> Vec<&str> {
        self.header
            .iter()
            .find(|r| r.name == "FILE_SCHEMA")
            .and_then(|r| r.params.first())
            .and_then(|v| v.as_list())
            .map(|items| items.iter().filter_map(|v| v.as_str()).collect())
            .unwrap_or_default()
    }
}
