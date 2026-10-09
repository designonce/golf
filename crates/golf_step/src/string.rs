//! Part 21 string escapes.

/// Decodes a string's contents (between its quotes, with `''` still doubled):
/// `''` is a quote, `\\` a backslash, `\X2\…\X0\` UTF-16 and `\X4\…\X0\` UTF-32
/// code units in hex, `\X\hh` one Latin-1 byte, and `\S\c` the character `c`
/// with its top bit set. Code-page directives (`\PA\`) are dropped. Bytes that
/// aren't UTF-8 are read as Latin-1.
pub(crate) fn decode(raw: &[u8]) -> String {
    let mut out = String::with_capacity(raw.len());
    let mut i = 0;
    let latin1 = |b: u8| char::from(b);
    while i < raw.len() {
        let rest = &raw[i..];
        if rest.starts_with(b"''") {
            out.push('\'');
            i += 2;
        } else if rest.starts_with(b"\\\\") {
            out.push('\\');
            i += 2;
        } else if rest.starts_with(b"\\X2\\") || rest.starts_with(b"\\X4\\") {
            let width = if rest[2] == b'2' { 4 } else { 8 };
            let body = &rest[4..];
            let end = find(body, b"\\X0\\").unwrap_or(body.len());
            let units: Vec<u32> = body[..end]
                .chunks(width)
                .filter_map(|c| u32::from_str_radix(core::str::from_utf8(c).ok()?, 16).ok())
                .collect();
            match width {
                4 => {
                    let units: Vec<u16> = units.iter().map(|&u| u as u16).collect();
                    out.extend(char::decode_utf16(units).map(|c| c.unwrap_or('\u{fffd}')));
                }
                _ => out.extend(
                    units
                        .iter()
                        .map(|&u| char::from_u32(u).unwrap_or('\u{fffd}')),
                ),
            }
            i += 4 + (end + 4).min(body.len());
        } else if rest.starts_with(b"\\X\\") && rest.len() >= 5 {
            match core::str::from_utf8(&rest[3..5])
                .ok()
                .and_then(|h| u8::from_str_radix(h, 16).ok())
            {
                Some(b) => out.push(latin1(b)),
                None => out.push('\u{fffd}'),
            }
            i += 5;
        } else if rest.starts_with(b"\\S\\") && rest.len() >= 4 {
            out.push(latin1(rest[3] | 0x80));
            i += 4;
        } else if rest.len() >= 4 && rest[0] == b'\\' && rest[1] == b'P' && rest[3] == b'\\' {
            i += 4;
        } else {
            // Copy up to the next escape or quote, as UTF-8 if it is.
            let len = rest
                .iter()
                .position(|&b| b == b'\\' || b == b'\'')
                .map_or(rest.len(), |p| p.max(1));
            let run = &rest[..len];
            match core::str::from_utf8(run) {
                Ok(s) => out.push_str(s),
                Err(_) => out.extend(run.iter().map(|&b| latin1(b))),
            }
            i += len;
        }
    }
    out
}

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack.windows(needle.len()).position(|w| w == needle)
}

#[cfg(test)]
mod tests {
    use super::decode;

    #[test]
    fn escapes_decode() {
        assert_eq!(decode(b"it''s"), "it's");
        assert_eq!(decode(b"a\\\\b"), "a\\b");
        assert_eq!(decode(b"\\X2\\00E9006C00E8\\X0\\ve"), "élève");
        assert_eq!(decode(b"\\X4\\0001F600\\X0\\"), "😀");
        assert_eq!(decode(b"caf\\X\\E9"), "café");
        assert_eq!(decode(b"\\S\\i"), "é");
        assert_eq!(decode(b"\\PA\\x"), "x");
        assert_eq!(decode("ünï".as_bytes()), "ünï");
        assert_eq!(decode(b"caf\xe9"), "café");
    }
}
