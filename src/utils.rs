
const EPSILON: f64 = 1e-6;

#[inline]
pub const fn fuzzy_equals(a: f64, b: f64) -> bool {
    (a - b).abs() < EPSILON
}

pub trait MaxScore {
    fn max_score(self) -> f64;
}

impl MaxScore for &[f64] {
    fn max_score(self) -> f64 {
        self.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
}