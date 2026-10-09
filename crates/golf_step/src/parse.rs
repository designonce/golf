//! The Part 21 grammar: sections of instances `#n = KEYWORD(params);`.

use crate::data::StepData;
use crate::error::ParseError;
use crate::error::ParseIssue;
use crate::string::decode;
use crate::value::Entity;
use crate::value::Id;
use crate::value::Record;
use crate::value::Value;

pub(crate) fn parse(bytes: &[u8]) -> Result<StepData, ParseError> {
    let mut p = Parser { bytes, pos: 0 };
    p.skip();
    if !p.eat_word("ISO-10303-21") {
        return Err(ParseError::NotPart21);
    }
    p.expect(b';')?;
    let mut data = StepData::default();
    loop {
        p.skip();
        if p.at_end() || p.eat_word("END-ISO-10303-21") {
            break;
        }
        if p.eat_word("HEADER") {
            p.expect(b';')?;
            while !p.section_end() {
                let start = p.pos;
                match p.record().and_then(|r| p.expect(b';').map(|_| r)) {
                    Ok(record) => data.header.push(record),
                    Err(error) => {
                        data.issues.push(p.issue(None, start, error));
                        p.recover();
                    }
                }
            }
        } else if p.eat_word("DATA") {
            // `DATA;`, or `DATA('name', (schema));` in later editions.
            p.skip();
            if p.peek() == Some(b'(') {
                p.value()?;
            }
            p.expect(b';')?;
            while !p.section_end() {
                let start = p.pos;
                match p.instance() {
                    Ok(entity) => {
                        data.entities.insert(entity.id, entity);
                    }
                    Err((id, error)) => {
                        data.issues.push(p.issue(id, start, error));
                        p.recover();
                    }
                }
            }
        } else if p.eat_word("ANCHOR") || p.eat_word("REFERENCE") || p.eat_word("SIGNATURE") {
            // Sections this reader has no use for: skip to their ENDSEC.
            while !p.section_end() {
                p.recover();
            }
        } else {
            return Err(p.error("expected HEADER, DATA or END-ISO-10303-21"));
        }
    }
    data.entities.shrink_to_fit();
    Ok(data)
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl Parser<'_> {
    fn at_end(&self) -> bool {
        self.pos >= self.bytes.len()
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    /// Skips whitespace and `/* comments */`.
    fn skip(&mut self) {
        loop {
            while self.peek().is_some_and(|b| b.is_ascii_whitespace()) {
                self.pos += 1;
            }
            if self.bytes[self.pos..].starts_with(b"/*") {
                match find(&self.bytes[self.pos + 2..], b"*/") {
                    Some(end) => self.pos += 2 + end + 2,
                    None => self.pos = self.bytes.len(),
                }
            } else {
                return;
            }
        }
    }

    fn eat_word(&mut self, word: &str) -> bool {
        self.skip();
        let rest = &self.bytes[self.pos..];
        let after = rest.get(word.len()).copied();
        let ends = !after.is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-');
        if rest.starts_with(word.as_bytes()) && ends {
            self.pos += word.len();
            true
        } else {
            false
        }
    }

    fn expect(&mut self, byte: u8) -> Result<(), ParseError> {
        self.skip();
        if self.peek() == Some(byte) {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.error(&format!("expected '{}'", byte as char)))
        }
    }

    /// At `ENDSEC;` (consumed) or the end of the input.
    fn section_end(&mut self) -> bool {
        self.skip();
        if self.eat_word("ENDSEC") {
            let _ = self.expect(b';');
            return true;
        }
        self.at_end()
    }

    fn line(&self, pos: usize) -> usize {
        1 + self.bytes[..pos.min(self.bytes.len())]
            .iter()
            .filter(|&&b| b == b'\n')
            .count()
    }

    fn error(&self, message: &str) -> ParseError {
        let near: String = self.bytes[self.pos..]
            .iter()
            .take(24)
            .map(|&b| b as char)
            .collect();
        ParseError::Syntax {
            line: self.line(self.pos),
            message: format!("{message}, near {near:?}"),
        }
    }

    fn issue(&self, id: Option<Id>, start: usize, error: ParseError) -> ParseIssue {
        ParseIssue {
            id,
            line: self.line(start),
            message: match error {
                ParseError::Syntax { message, .. } => message,
                other => other.to_string(),
            },
        }
    }

    /// Skips past the next `;` outside strings, comments and binary.
    fn recover(&mut self) {
        while let Some(b) = self.peek() {
            match b {
                b'\'' => {
                    self.pos += 1;
                    let _ = self.string_body();
                }
                b'"' => {
                    self.pos += 1;
                    while self.peek().is_some_and(|b| b != b'"') {
                        self.pos += 1;
                    }
                    self.pos += 1;
                }
                b'/' if self.bytes[self.pos..].starts_with(b"/*") => self.skip(),
                b';' => {
                    self.pos += 1;
                    return;
                }
                _ => self.pos += 1,
            }
        }
    }

    /// `#n = record;` or `#n = (record record ...);`.
    fn instance(&mut self) -> Result<Entity, (Option<Id>, ParseError)> {
        let id = self.reference().map_err(|e| (None, e))?;
        let with_id = |e| (Some(id), e);
        self.expect(b'=').map_err(with_id)?;
        self.skip();
        let records = if self.peek() == Some(b'(') {
            self.pos += 1;
            let mut records = Vec::new();
            loop {
                self.skip();
                if self.peek() == Some(b')') {
                    self.pos += 1;
                    break;
                }
                records.push(self.record().map_err(with_id)?);
            }
            records
        } else {
            vec![self.record().map_err(with_id)?]
        };
        self.expect(b';').map_err(with_id)?;
        Ok(Entity { id, records })
    }

    fn reference(&mut self) -> Result<Id, ParseError> {
        self.expect(b'#')?;
        let start = self.pos;
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            self.pos += 1;
        }
        core::str::from_utf8(&self.bytes[start..self.pos])
            .ok()
            .and_then(|s| s.parse().ok())
            .map(Id)
            .ok_or_else(|| self.error("expected an instance number"))
    }

    fn keyword(&mut self) -> Result<String, ParseError> {
        self.skip();
        let start = self.pos;
        // User-defined keywords start with '!'.
        if self.peek() == Some(b'!') {
            self.pos += 1;
        }
        while self
            .peek()
            .is_some_and(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
        {
            self.pos += 1;
        }
        if self.pos == start {
            return Err(self.error("expected a keyword"));
        }
        Ok(String::from_utf8_lossy(&self.bytes[start..self.pos]).to_ascii_uppercase())
    }

    /// `KEYWORD(params)`.
    fn record(&mut self) -> Result<Record, ParseError> {
        let name = self.keyword()?;
        self.expect(b'(')?;
        let params = self.list_items()?;
        Ok(Record { name, params })
    }

    /// The items of a list whose `(` has been read, and its `)`.
    fn list_items(&mut self) -> Result<Vec<Value>, ParseError> {
        let mut items = Vec::new();
        self.skip();
        if self.peek() == Some(b')') {
            self.pos += 1;
            return Ok(items);
        }
        loop {
            items.push(self.value()?);
            self.skip();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b')') => {
                    self.pos += 1;
                    return Ok(items);
                }
                _ => return Err(self.error("expected ',' or ')'")),
            }
        }
    }

    fn value(&mut self) -> Result<Value, ParseError> {
        self.skip();
        let Some(b) = self.peek() else {
            return Err(self.error("unexpected end of file"));
        };
        match b {
            b'$' => {
                self.pos += 1;
                Ok(Value::Null)
            }
            b'*' => {
                self.pos += 1;
                Ok(Value::Derived)
            }
            b'#' => self.reference().map(Value::Ref),
            b'\'' => {
                self.pos += 1;
                self.string_body().map(Value::String)
            }
            b'"' => {
                self.pos += 1;
                let start = self.pos;
                while self.peek().is_some_and(|b| b != b'"') {
                    self.pos += 1;
                }
                let text = String::from_utf8_lossy(&self.bytes[start..self.pos]).into_owned();
                self.expect(b'"')?;
                Ok(Value::Binary(text))
            }
            b'.' if self
                .bytes
                .get(self.pos + 1)
                .is_some_and(|b| b.is_ascii_alphabetic() || *b == b'_') =>
            {
                self.pos += 1;
                let start = self.pos;
                while self.peek().is_some_and(|b| b != b'.') {
                    self.pos += 1;
                }
                let name =
                    String::from_utf8_lossy(&self.bytes[start..self.pos]).to_ascii_uppercase();
                self.expect(b'.')?;
                Ok(Value::Enum(name))
            }
            b'(' => {
                self.pos += 1;
                self.list_items().map(Value::List)
            }
            b'+' | b'-' | b'.' | b'0'..=b'9' => self.number(),
            b if b.is_ascii_alphabetic() || b == b'!' => {
                let name = self.keyword()?;
                self.expect(b'(')?;
                let mut items = self.list_items()?;
                let inner = match items.len() {
                    1 => items.pop().expect("one item"),
                    _ => Value::List(items),
                };
                Ok(Value::Typed(name, Box::new(inner)))
            }
            _ => Err(self.error("expected a value")),
        }
    }

    /// The rest of a string whose opening quote has been read.
    fn string_body(&mut self) -> Result<String, ParseError> {
        let start = self.pos;
        loop {
            match self.peek() {
                None => return Err(self.error("unterminated string")),
                Some(b'\'') if self.bytes.get(self.pos + 1) == Some(&b'\'') => self.pos += 2,
                Some(b'\'') => break,
                Some(_) => self.pos += 1,
            }
        }
        let text = decode(&self.bytes[start..self.pos]);
        self.pos += 1;
        Ok(text)
    }

    fn number(&mut self) -> Result<Value, ParseError> {
        let start = self.pos;
        if matches!(self.peek(), Some(b'+' | b'-')) {
            self.pos += 1;
        }
        let mut real = false;
        while let Some(b) = self.peek() {
            match b {
                b'0'..=b'9' => {}
                b'.' => real = true,
                b'E' | b'e' => {
                    real = true;
                    if matches!(self.bytes.get(self.pos + 1), Some(b'+' | b'-')) {
                        self.pos += 1;
                    }
                }
                _ => break,
            }
            self.pos += 1;
        }
        let text = core::str::from_utf8(&self.bytes[start..self.pos]).expect("ascii");
        let parsed = match real {
            // Rust doesn't read "1." or "1.E3", which Part 21 writes.
            true => text
                .replace(".E", ".0E")
                .replace(".e", ".0e")
                .trim_end_matches('.')
                .parse()
                .ok()
                .map(Value::Real),
            false => text.parse().ok().map(Value::Integer),
        };
        parsed.ok_or_else(|| self.error("expected a number"))
    }
}

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('a test'),'2;1');
FILE_NAME('cube.step','2026-10-09T12:00:00',('me'),(''),'golf','golf','');
FILE_SCHEMA(('AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF { 1 0 10303 442 1 1 4 }'));
ENDSEC;
DATA;
/* a comment */
#1 = CARTESIAN_POINT('',(0.,1.5,-2.E-3));
#2=DIRECTION('',(0.0,0.0,1.0));
#3 = AXIS2_PLACEMENT_3D('', #1, #2, $);
#10 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#11 = LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#10);
#12 = PRODUCT('it''s a part','name','',(#13));
#20 = BROKEN(#1,;
#21 = VERTEX_POINT('',#1);
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn a_file_parses() {
        let data = StepData::parse(FILE.as_bytes()).unwrap();
        assert_eq!(data.len(), 7);
        assert_eq!(data.header.len(), 3);
        assert_eq!(
            data.schemas()[0].split(' ').next(),
            Some("AP242_MANAGED_MODEL_BASED_3D_ENGINEERING_MIM_LF")
        );

        let point = data.record(Id(1), "CARTESIAN_POINT").unwrap();
        assert_eq!(point.reals(1).unwrap(), vec![0.0, 1.5, -2e-3]);
        let axes = data.record(Id(3), "AXIS2_PLACEMENT_3D").unwrap();
        assert_eq!(axes.reference(1).unwrap(), Id(1));
        assert_eq!(axes.optional_reference(3).unwrap(), None);

        let unit = data.entity(Id(10)).unwrap();
        assert!(unit.is("LENGTH_UNIT") && unit.is("SI_UNIT"));
        assert_eq!(
            unit.expect("SI_UNIT").unwrap().enumeration(0).unwrap(),
            "MILLI"
        );
        assert_eq!(unit.name(), "LENGTH_UNIT+NAMED_UNIT+SI_UNIT");

        let measure = data.record(Id(11), "LENGTH_MEASURE_WITH_UNIT").unwrap();
        assert_eq!(measure.real(0).unwrap(), 25.4);
        assert!(
            matches!(measure.param(0).unwrap(), Value::Typed(name, _) if name == "LENGTH_MEASURE")
        );
        assert_eq!(
            data.record(Id(12), "PRODUCT").unwrap().string(0).unwrap(),
            "it's a part"
        );

        // The broken instance is skipped, and what follows it still read.
        assert_eq!(data.issues.len(), 1);
        assert_eq!(data.issues[0].id, Some(Id(20)));
        assert!(data.get(Id(21)).is_some());
    }

    #[test]
    fn decode_errors_say_what_was_wrong() {
        let data = StepData::parse(FILE.as_bytes()).unwrap();
        let point = data.record(Id(1), "CARTESIAN_POINT").unwrap();
        assert_eq!(
            point.reference(1).unwrap_err().to_string(),
            "CARTESIAN_POINT parameter 1: expected a reference, found a list of 3"
        );
        assert_eq!(
            data.record(Id(2), "VERTEX_POINT").unwrap_err().to_string(),
            "#2 is a DIRECTION, not a VERTEX_POINT"
        );
        assert_eq!(
            data.entity(Id(99)).unwrap_err().to_string(),
            "#99 is missing"
        );
    }

    #[test]
    fn other_text_is_rejected() {
        assert_eq!(
            StepData::parse(b"solid cube").unwrap_err(),
            ParseError::NotPart21
        );
    }
}
