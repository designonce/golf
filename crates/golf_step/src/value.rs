use core::fmt;

use crate::error::DecodeError;

/// An entity instance's id: the `n` of `#n`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Id(pub u64);

impl fmt::Display for Id {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(f)
    }
}

/// A parameter value.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    /// `$`: no value.
    Null,
    /// `*`: derived from elsewhere.
    Derived,
    Integer(i64),
    Real(f64),
    /// A string, with its escapes decoded.
    String(String),
    /// `.NAME.`, without the dots.
    Enum(String),
    /// `"..."`, as written.
    Binary(String),
    Ref(Id),
    List(Vec<Value>),
    /// `NAME(value)`: a value of a defined type, e.g. `LENGTH_MEASURE(25.4)`.
    Typed(String, Box<Value>),
}

/// One keyword and its parameters: a simple entity, or one part of a complex
/// one.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    /// Upper case, as STEP keywords are.
    pub name: String,
    pub params: Vec<Value>,
}

/// An entity instance: one record, or several for a complex instance
/// `(A(...) B(...))`, in the order written.
#[derive(Clone, Debug, PartialEq)]
pub struct Entity {
    pub id: Id,
    pub records: Vec<Record>,
}

impl Entity {
    /// Whether this is a simple instance named `name`, or a complex one with a
    /// part named `name`.
    pub fn is(&self, name: &str) -> bool {
        self.records.iter().any(|r| r.name == name)
    }

    /// The record named `name`.
    pub fn record(&self, name: &str) -> Option<&Record> {
        self.records.iter().find(|r| r.name == name)
    }

    /// The record named `name`, or an error naming what was found.
    pub fn expect(&self, name: &str) -> Result<&Record, DecodeError> {
        self.record(name).ok_or_else(|| DecodeError::WrongEntity {
            id: self.id,
            expected: name.to_string(),
            found: self.name(),
        })
    }

    /// The first record of a simple instance, which names it.
    pub fn simple(&self) -> Option<&Record> {
        match self.records.as_slice() {
            [record] => Some(record),
            _ => None,
        }
    }

    /// Its keyword, or its keywords joined by `+` for a complex instance.
    pub fn name(&self) -> String {
        self.records
            .iter()
            .map(|r| r.name.as_str())
            .collect::<Vec<_>>()
            .join("+")
    }
}

impl Record {
    /// Parameter `index`.
    pub fn param(&self, index: usize) -> Result<&Value, DecodeError> {
        self.params
            .get(index)
            .ok_or_else(|| DecodeError::MissingParameter {
                entity: self.name.clone(),
                index,
            })
    }

    fn wrong(&self, index: usize, expected: &'static str) -> DecodeError {
        DecodeError::WrongValue {
            entity: self.name.clone(),
            index,
            expected,
            found: self
                .params
                .get(index)
                .map_or("nothing".to_string(), Value::describe),
        }
    }

    pub fn real(&self, index: usize) -> Result<f64, DecodeError> {
        self.param(index)?
            .as_real()
            .ok_or_else(|| self.wrong(index, "a number"))
    }

    pub fn integer(&self, index: usize) -> Result<i64, DecodeError> {
        self.param(index)?
            .as_integer()
            .ok_or_else(|| self.wrong(index, "an integer"))
    }

    pub fn reference(&self, index: usize) -> Result<Id, DecodeError> {
        self.param(index)?
            .as_ref_id()
            .ok_or_else(|| self.wrong(index, "a reference"))
    }

    /// A reference, or `None` for `$`.
    pub fn optional_reference(&self, index: usize) -> Result<Option<Id>, DecodeError> {
        match self.param(index)? {
            Value::Null => Ok(None),
            value => value
                .as_ref_id()
                .map(Some)
                .ok_or_else(|| self.wrong(index, "a reference or $")),
        }
    }

    pub fn string(&self, index: usize) -> Result<&str, DecodeError> {
        self.param(index)?
            .as_str()
            .ok_or_else(|| self.wrong(index, "a string"))
    }

    pub fn enumeration(&self, index: usize) -> Result<&str, DecodeError> {
        self.param(index)?
            .as_enum()
            .ok_or_else(|| self.wrong(index, "an enumeration"))
    }

    /// `.T.` or `.F.`.
    pub fn boolean(&self, index: usize) -> Result<bool, DecodeError> {
        self.param(index)?
            .as_bool()
            .ok_or_else(|| self.wrong(index, "a boolean"))
    }

    pub fn list(&self, index: usize) -> Result<&[Value], DecodeError> {
        self.param(index)?
            .as_list()
            .ok_or_else(|| self.wrong(index, "a list"))
    }

    /// A list of references.
    pub fn references(&self, index: usize) -> Result<Vec<Id>, DecodeError> {
        self.list(index)?
            .iter()
            .map(|v| {
                v.as_ref_id()
                    .ok_or_else(|| self.wrong(index, "a list of references"))
            })
            .collect()
    }

    /// A list of numbers.
    pub fn reals(&self, index: usize) -> Result<Vec<f64>, DecodeError> {
        self.list(index)?
            .iter()
            .map(|v| {
                v.as_real()
                    .ok_or_else(|| self.wrong(index, "a list of numbers"))
            })
            .collect()
    }

    /// A list of integers.
    pub fn integers(&self, index: usize) -> Result<Vec<i64>, DecodeError> {
        self.list(index)?
            .iter()
            .map(|v| {
                v.as_integer()
                    .ok_or_else(|| self.wrong(index, "a list of integers"))
            })
            .collect()
    }
}

impl Value {
    /// A real, or an integer as one; looks through a typed value.
    pub fn as_real(&self) -> Option<f64> {
        match self {
            Value::Real(x) => Some(*x),
            Value::Integer(n) => Some(*n as f64),
            Value::Typed(_, inner) => inner.as_real(),
            _ => None,
        }
    }

    pub fn as_integer(&self) -> Option<i64> {
        match self {
            Value::Integer(n) => Some(*n),
            Value::Typed(_, inner) => inner.as_integer(),
            _ => None,
        }
    }

    pub fn as_ref_id(&self) -> Option<Id> {
        match self {
            Value::Ref(id) => Some(*id),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::String(s) => Some(s),
            Value::Typed(_, inner) => inner.as_str(),
            _ => None,
        }
    }

    pub fn as_enum(&self) -> Option<&str> {
        match self {
            Value::Enum(e) => Some(e),
            _ => None,
        }
    }

    pub fn as_bool(&self) -> Option<bool> {
        match self.as_enum()? {
            "T" => Some(true),
            "F" => Some(false),
            _ => None,
        }
    }

    pub fn as_list(&self) -> Option<&[Value]> {
        match self {
            Value::List(items) => Some(items),
            _ => None,
        }
    }

    /// A short description for error messages.
    pub fn describe(&self) -> String {
        match self {
            Value::Null => "$".to_string(),
            Value::Derived => "*".to_string(),
            Value::Integer(n) => n.to_string(),
            Value::Real(x) => format!("{x:?}"),
            Value::String(s) => format!("'{s}'"),
            Value::Enum(e) => format!(".{e}."),
            Value::Binary(b) => format!("\"{b}\""),
            Value::Ref(id) => format!("#{id}"),
            Value::List(items) => format!("a list of {}", items.len()),
            Value::Typed(name, _) => format!("{name}(...)"),
        }
    }
}
