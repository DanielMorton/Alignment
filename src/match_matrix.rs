use std::collections::HashMap;

#[derive(Debug, Clone, Default)]
pub struct MatchMatrix(HashMap<char, HashMap<char, f64>>);

impl MatchMatrix {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn set_score(&mut self, a: char, b: char, score: f64) {
        self.0.entry(a).or_default().insert(b, score);
    }

    pub(crate) fn get_score(&self, a: char, b: char) -> f64 {
        self.0
            .get(&a)
            .and_then(|inner| inner.get(&b))
            .copied()
            .unwrap_or(0.0)
    }
}