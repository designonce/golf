use core::fmt;
use core::str::FromStr;

/// An sRGB colour, 8 bits a channel, optionally with alpha.
///
/// Alpha is `None` when nothing was said about transparency, which is distinct
/// from an explicit fully opaque `Some(255)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: Option<u8>,
}

impl Color {
    /// An opaque colour.
    pub const fn rgb(r: u8, g: u8, b: u8) -> Self {
        Self { r, g, b, a: None }
    }

    /// A colour with alpha, 0 transparent to 255 opaque.
    pub const fn rgba(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self {
            r,
            g,
            b,
            a: Some(a),
        }
    }

    /// Red, green and blue as fractions in [0, 1].
    pub fn unit_rgb(&self) -> [f64; 3] {
        [self.r, self.g, self.b].map(|c| f64::from(c) / 255.0)
    }

    /// Opacity as a fraction in [0, 1]: 1 if unspecified.
    pub fn opacity(&self) -> f64 {
        self.a.map_or(1.0, |a| f64::from(a) / 255.0)
    }
}

/// Why text isn't a colour.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{0:?} is not a colour; expected #rrggbb or #rrggbbaa")]
pub struct ParseColorError(pub String);

impl FromStr for Color {
    type Err = ParseColorError;

    /// Parses `#rrggbb` or `#rrggbbaa` (the `#` is optional).
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        let error = || ParseColorError(text.to_string());
        let hex = text.strip_prefix('#').unwrap_or(text);
        if !matches!(hex.len(), 6 | 8) || !hex.is_ascii() {
            return Err(error());
        }
        let channel = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).map_err(|_| error());
        let (r, g, b) = (channel(0)?, channel(2)?, channel(4)?);
        Ok(match hex.len() {
            8 => Self::rgba(r, g, b, channel(6)?),
            _ => Self::rgb(r, g, b),
        })
    }
}

impl fmt::Display for Color {
    /// `#rrggbb`, or `#rrggbbaa` with alpha.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "#{:02x}{:02x}{:02x}", self.r, self.g, self.b)?;
        if let Some(a) = self.a {
            write!(f, "{a:02x}")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_round_trips() {
        let c: Color = "#ff8000".parse().unwrap();
        assert_eq!(c, Color::rgb(255, 128, 0));
        assert_eq!(c.to_string(), "#ff8000");
        let c: Color = "20406080".parse().unwrap();
        assert_eq!(c, Color::rgba(0x20, 0x40, 0x60, 0x80));
        assert_eq!(c.to_string(), "#20406080");
        for bad in ["", "#fff", "#gg0000", "#ff00000", "#ffé000"] {
            assert!(bad.parse::<Color>().is_err(), "{bad}");
        }
    }

    #[test]
    fn unit_channels_and_opacity() {
        assert_eq!(Color::rgb(255, 0, 51).unit_rgb(), [1.0, 0.0, 0.2]);
        assert_eq!(Color::rgb(1, 2, 3).opacity(), 1.0);
        assert_eq!(Color::rgba(1, 2, 3, 0).opacity(), 0.0);
    }
}
